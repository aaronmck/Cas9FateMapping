import scala.io._
import java.io._
import scala.collection.mutable._
import aaronmck.github.io.seq.crispr._
import scala.util.matching.Regex

val tearsheet = Source.fromFile(args(0)).getLines() // .drop(1)
val output = new PrintWriter(args(1))
val cutSites = args(2)

// constants we use in the computation
val window = 5
val maxCigarLength = 40
val minReadCount = 10

// load up the cutsite information so we can map indels onto it
def loadCutSites(cutSitePositions: String, cutSiteWindow: Int): HashMap[Int,Tuple3[Int,Int,Int]] = {
  val cutPos = new HashMap[Int,Tuple3[Int,Int,Int]]()
  Source.fromFile(cutSitePositions).getLines().drop(1).zipWithIndex.foreach{case(line,index) =>
    cutPos(index) = (line.split("\t")(2).toInt - cutSiteWindow,line.split("\t")(2).toInt,line.split("\t")(2).toInt + cutSiteWindow)
  }
  return(cutPos)
}

// create our cut site map
val cutPositionsArray = loadCutSites(cutSites, window)

// setup the output file
output.write("sample\tUMI\tused\tfail.reason\tfwd.reads\tprimered.fwd.reads\tfiltered.fwd.reads\trev.reads\tprimered.rev.reads\tfiltered.rev.reads\t")
output.write("fwd.error.rate\tfwd.read\trev.error.rate\trev.read\tmerged.cigar")
(1 until cutPositionsArray.size + 1).foreach{index => output.write("\ttarget" + index)}
(1 until cutPositionsArray.size + 1).foreach{index => output.write("\tencoded.target" + index)}
output.write("\n")

// setup a way to encode alleles
val encodings = HashMap[Int,HashMap[String,Int]]()
val encodingsNext = HashMap[Int,Int]()
cutPositionsArray.zipWithIndex.foreach{case(positions,index) => {
  encodings(index) = new HashMap[String,Int]()
  encodingsNext(index) = 1
}}

// now for each file in the tearsheet do the following: load the barcodes from the sam file,
// load and add those to the appropriate row of the stats file
tearsheet.foreach{line => {
  val sp = line.split("\t")
  val outputDir = sp(3)
  val sample = sp(0)

  val statsFile = new File(outputDir + "/" + sample + "/" + sample + ".stats")
  val alignmentFile = new File(outputDir + "/" + sample + "/" + sample + ".sam")

  if (statsFile.exists && alignmentFile.exists) {
    println("processing " + alignmentFile)
    val barcodeToCigar = new HashMap[String,String]()

    var lines = 0
    // read in the sam file
    Source.fromFile(alignmentFile).getLines().foreach{line => {
      if (!(line startsWith "@")) {
        val tokens = line.split("\t")
        barcodeToCigar(tokens(0).split("_")(0)) = tokens(5)
      }
      lines += 1
    }}

    println("barcode size: " + barcodeToCigar.size + " lines " + lines)

    // create some storage for our cut indels
    val cutIndels = IndelMap(cutPositionsArray)
    var index = 0

    val counts = Source.fromFile(statsFile).getLines().drop(1).map{line => line.split("\t")(3).toInt}.toArray
    val mean = counts.sum.toDouble / counts.size
    val stdDevTmp = counts.map{cnt => (cnt - mean) * (cnt - mean) }.toArray
    val stdDev = math.sqrt(stdDevTmp.sum / stdDevTmp.length.toDouble)

    println("Processing " + statsFile + " with mean/std of " + mean + "/" + stdDev)

    var total = 0
    var kept = 0

    val haveWeSeenIt = new HashMap[String,Boolean]()

    Source.fromFile(statsFile).getLines().drop(1).foreach{line => {
      val umi = line.split("\t")(0)
      val fwdCount = line.split("\t")(3).toInt
      val revCount = line.split("\t")(6).toInt
      total += 1
      val cigar = barcodeToCigar.getOrElse(umi,"NA")

      if (!(cigar contains "NA") &&
        cigar.length < maxCigarLength &&
        fwdCount < (mean + (stdDev * 2)) &&
        fwdCount > (mean - (stdDev * 2)) &&
        fwdCount > minReadCount) {



          //if (cigar == "NA") println("failed on " + line)
          val hit = cutIndels.addToIndelMapping(cigar)

          // map the events to integer values
          val events = hit._2.zipWithIndex.map{case(cigarEvent,index) => {
            if (cigarEvent == "") {
              0
            } else {
              if (encodings(index) contains cigarEvent)
                encodings(index)(cigarEvent)
              else {
                encodings(index)(cigarEvent) = encodingsNext(index)
                encodingsNext(index) = encodingsNext(index) + 1
                encodings(index)(cigarEvent)
              }
            }
          }}
          val eventStr = hit._2.map{h => if (h == "") "NONE" else h}.mkString("-")
          if (!(haveWeSeenIt contains eventStr)) {
            kept += 1
            output.write(sample + "\t" + line + "\t" + cigar + "\t")
            output.write(hit._2.map{h => if (h == "") "NONE" else h}.mkString("\t") + "\t")
            output.write(events.mkString("\t") + "\n")
          }
          haveWeSeenIt(eventStr) = true

      } else {
        //println("Dropping " + umi + " fwd count " + fwdCount + " mean " + mean + " stddev " + stdDev)
      }
    }}
    println("total = " + total + " kept = " + kept + " prop = " + (kept.toDouble/total.toDouble))
  }
}}

output.close()
