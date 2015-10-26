import scala.io._
import java.io._
import scala.collection.mutable._
import aaronmck.github.io.seq.crispr._

val tearsheet = Source.fromFile(args(0)).getLines() // .drop(1)
val output = new PrintWriter(args(1))
val cutSites = args(2)
val window = 5

// load up the cutsite information for each of the cutsite maps
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
output.write("sample\tUMI\tall.treads.with.umi\tumi.reads.with.primers\treads.used.in.final.consensus")
output.write("\tmean.proportion.of.final.bases.matching.consensus\treads.retained.primers.to.final\t")
output.write("umi.used.in.analysis\tfinal.consensus\tref\tcigar")
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
        barcodeToCigar(tokens(0).split("_")(3)) = tokens(5)
      }
      lines += 1
    }}

    println("barcode size: " + barcodeToCigar.size + " lines " + lines)

    // create some storage for our cut indels
    val cutIndels = IndelMap(cutPositionsArray)
    var index = 0

    Source.fromFile(statsFile).getLines().drop(1).foreach{line => {
      val umi = line.split("\t")(0)
      val cigar = barcodeToCigar.getOrElse(umi,"NA")

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

      output.write(sample + "\t" + line + "\t" + cigar + "\t")
      output.write(hit._2.map{h => if (h == "") "NONE" else h}.mkString("\t") + "\t")
      output.write(events.mkString("\t") + "\n")

    }}
  }
}}

output.close()
