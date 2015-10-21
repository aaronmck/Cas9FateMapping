import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._

// ------------------------------------------------------------------------------------------
//
// for each read, count the reads by UMI, output only those with an high enough number of
// UMI occurances
//
// ------------------------------------------------------------------------------------------

// store the results of our read merging
case class MergeResult(umi: String, consensusRead: String, qualities: Array[Double], overallQuality: Double)

case class Merger(alignedFile: File) {
  val nameToRead = new HashMap[String,String]()
  var currentReadName: Option[String] = None
  var currentRead = ""
  var readLength = 0

  Source.fromFile(alignedFile).getLines().foreach{line => {
    if (line startsWith ">") {
      if (currentReadName.isDefined) {
        nameToRead(currentReadName.get) = currentRead
        readLength = currentRead.length
        currentRead = ""
        currentReadName = Some(line.stripPrefix(">"))
      } else {
        currentReadName = Some(line.stripPrefix(">"))
      }
    } else {
      currentRead += line
    }
  }}

  if (currentReadName.isDefined)
      nameToRead(currentReadName.get) = currentRead
}


// develop a consensus over the read -- look at each aligned position and find the most common character and it's proportion
def consensus(nameToRead : HashMap[String,String], readLength: Int): Tuple2[String,Array[Double]] = {
    // each string should be the same length coming back from CLUSTAL
    var highestOccuringBase = Array[Char]()
    var highestOccuringBaseProp = Array[Double]()

    (0 until readLength).foreach{index => {
      var bases = new HashMap[Char,Int]()

      nameToRead.foreach{case(name,read) => {bases(read(index)) = bases.getOrElse(read(index),0) + 1 }}
      var maxChar = 'U'
      var total = 0
      bases.foreach{case(base,count) => {
        if (!(bases contains maxChar) || count > bases(maxChar)) maxChar = base
        total += count
      }}
      highestOccuringBase :+= maxChar
      highestOccuringBaseProp :+= (if (bases contains maxChar) bases(maxChar).toDouble / total.toDouble else 0.0)
    }}

    (highestOccuringBase.mkString(""),highestOccuringBaseProp)
  }

// ------------------------------------------------------------------------------------------
// helper functions
// ------------------------------------------------------------------------------------------

// do a first pass consensus, then remove any reads that are X bases (proportionally) different than the concensus read, then make a new consensus
def twoPassConsensus(nameToRead : HashMap[String,String], readLength: Int, differenceThreshold : Double): Tuple4[String,Array[Double],Int,Double] = {
    val firstPass = consensus(nameToRead,readLength)
    val newNameToRead = HashMap[String,String]()
    nameToRead.foreach{case(name,bases) => {
      val differences = firstPass._1.zip(bases).map{case(concBase,readBase) => if (readBase != concBase) 1 else 0}.sum.toDouble / bases.length.toDouble
      if (differences <= differenceThreshold)
        newNameToRead(name) = bases
    }}

    val conc = consensus(newNameToRead,readLength)
    return Tuple4[String,Array[Double], Int, Double](conc._1, conc._2, newNameToRead.size, newNameToRead.size.toDouble / nameToRead.size.toDouble)
}

// do a first pass consensus, then remove any reads that are X bases (proportionally) different than the concensus read, then make a new consensus
def collapseConsensus(consensusValues: Tuple2[String,Array[Double]]): Tuple2[String,Double] = {
    var retString = ""
    var retProps = Array[Double]()
    consensusValues._1.zipWithIndex.foreach{case(base,index) => {
      if (base != '-') {
        retString += base
        retProps :+= consensusValues._2(index)
      }
    }}
    (retString,retProps.sum.toDouble / retProps.length.toDouble)
}

def outputReads(umi: String, reads: Array[String], names: Array[String], outputFASTA: PrintWriter, outputStats: PrintWriter, ref: String, primers: List[String], sample: String) { //: MergeResult = {
  // filter down the reads to those that have the primers on each end
  var newReads = Array[String]()
  var newNames = Array[String]()
  reads.zip(names).foreach{case(read,name) => {
    if ((read contains primers(0)) && (read contains primers(1))) {
      newReads :+= read
      newNames :+= name
    }
  }}

  if (newReads.size < 2) {
    println("WARN: dropping umi " + umi + " since we have less than two surviving reads")
    return
  }

  // setup a CLUSTAL run and farm it out the machine
  val tmp = java.io.File.createTempFile("UMIMerger", ".txt")
  val tmpWriter = new PrintWriter(tmp)
  newReads.zip(newNames).foreach{case(read,name) => tmpWriter.write("> " + name + "\n" + read + "\n")}
  tmpWriter.close()
  val tmpOutput = java.io.File.createTempFile("UMIMerger", ".post.txt")
  val clustoResult = ("clustalo --force --pileup -i " + tmp + " -o " + tmpOutput).!!

  // read in the output from CLUSTAL and see how many reads we have
  val merger = Merger(tmpOutput)

  // merge the reads down to a single consensus read
  val consensusRead = twoPassConsensus(merger.nameToRead,merger.readLength,0.10)
  val collapsedConc = collapseConsensus((consensusRead._1,consensusRead._2))

  // require that the consensus read has kept at least 80% of the original reads
  val usingRead = (consensusRead._4 > 0.8) && (collapsedConc._1.slice(0,primers(0).length + 5) contains primers(0)) && (collapsedConc._1.slice(collapsedConc._1.size - (primers(0).length + 5),collapsedConc._1.size) contains primers(1))
  outputStats.write(umi + "\t" + reads.size + "\t" + newReads.size + "\t" + consensusRead._3 + "\t" + collapsedConc._2 + "\t" + consensusRead._4 + "\t" + usingRead + "\t" + collapsedConc._1 + "\t" + ref + "\n")
  if ( usingRead ) {
    outputFASTA.write(">" + sample + "_umi_" + umi + "_reads_used_" + consensusRead._3 + "_overall_match_prop_" + collapsedConc._2 + "_read_prop_retained_" + consensusRead._4 + "\n" + collapsedConc._1 + "\n")
  }
}

// *********************************************************************************************************
// print the usage parameters and quit
def printUsageAndQuit() {
  println("usage:\n\nscala umi_count.scala <reads.fq> <output_fasta> <output_stats> <umi_cutoff> <primers_file> <reference_file> <sample.name>")
  sys.exit(1)
}


// *********************************************************************************************************
// main script point
// *********************************************************************************************************

// check that the correct number of inputs was provided
if (args.length != 7) printUsageAndQuit()

// process the command line args
val inputFileReads = args(0)
val outputFile = new PrintWriter(args(1))
val outputStats = new PrintWriter(args(2))
val umiLength = args(3).toInt
val primersEachEnd = args(4)
val reference = args(5)
val samplename= args(6)

// get the reference string
var referenceString = ""
Source.fromFile(reference).getLines().foreach{line => if (!line.startsWith(">")) referenceString += line}

// first find out what UMIs to keep
// ------------------------------------------------------------------------------------------
val forwardReads = Source.fromFile(inputFileReads).getLines().grouped(4)

val primers = Source.fromFile(primersEachEnd).getLines().map{line => line}.toList
if (primers.length != 2)
  throw new IllegalStateException("You should only provide a primer file with two primers")

// our counts
val umiReads = new HashMap[String,ArrayBuilder[String]]()
val umiReadNames = new HashMap[String,ArrayBuilder[String]]()

// find the UMIs in the forward reads
var readsProcessed = 0
forwardReads foreach {grp => {
  val umi = grp(1).slice(0,umiLength)
  val readNoUMI = grp(1).slice(umiLength,grp(1).length)
  val readBuilder = umiReads.getOrElse(umi,ArrayBuilder.make[String])
  val nameBuilder = umiReadNames.getOrElse(umi,ArrayBuilder.make[String])
  readBuilder += readNoUMI
  nameBuilder += grp(0)
  umiReads(umi) = readBuilder
  umiReadNames(umi) = nameBuilder
  readsProcessed += 1
}}
//println("Processed " + readsProcessed + " reads....")

outputStats.write("UMI\tall.treads.with.umi\tumi.reads.with.primers\treads.used.in.final.consensus\tmean.proportion.of.final.bases.matching.consensus\treads.retained.primers.to.final\tumi.used.in.analysis\tfinal.consensus\tref\n")

var passingUMI = 0
var hereReads = 0
println("Processing " + umiReads.size + " UMIs")
umiReads.foreach{ case(umi,reads) => {

  val readArray = reads.result()
  val nameArray = umiReadNames(umi).result()

  hereReads += readArray.length
  if (readArray.length > 20) {
    outputReads(umi,readArray,nameArray, outputFile, outputStats, referenceString,primers, samplename)
    passingUMI += 1
  }
}}

outputStats.close()
outputFile.close()
