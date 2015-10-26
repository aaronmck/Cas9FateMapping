import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._

// read in compressed input streams with scala source commands
def gis(s: String) = new GZIPInputStream(new BufferedInputStream(new FileInputStream(s)))

// complement a base as a character
def compBase(b: Char): Char = b match {
  case 'A' => 'T'
  case 'C' => 'G'
  case 'G' => 'C'
  case 'T' => 'A'
  case 'a' => 't'
  case 'c' => 'g'
  case 'g' => 'c'
  case 't' => 'a'
  case _ => 'N'
}

// reverse complement a string of DNA bases
def reverseComplement(str: String) = str.map{t => compBase(t)}.reverse.mkString("")

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

def phredCharToDouble(ch: Char): Double = math.pow(10,((ch.toInt -32) / 10.0) * -1.0)
def phredCharToQscore(ch: Char): Int = (ch.toInt -32)

// quality filter reads.  Find any window of n bases that have a quality score less than X, cut and drop
// the rest of the read.  returns the remaining read string
def qualityControlRead(read: String, qualString: String, windowSize: Int = 5, minWindowQual: Double = 25): String = {
  val cutPos = qualString.toArray.sliding(windowSize).zipWithIndex.map{case(bases,index) => {
    if (bases.map{b => phredCharToQscore(b)}.sum / windowSize.toDouble < minWindowQual)
      index
    else
      0
  }}.filter(x => x != 0).toArray

  if (cutPos.size > 0)
    read.slice(0,cutPos(0))
  else
    read
}

// develop a consensus over the read -- look at each aligned position and find the most common character and it's proportion
// return the consensus read and proportion of bases at each position that match the called base
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
// returns a tuple-4: the consensus, it's average per-base error rate, number of reads used in the final assembly, and the proportion retained of the total
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
// return a tuple-2 of the consensus and it's average base level consensus
def collapseConsensus(consensusValues: Tuple2[String,Array[Double]]): Tuple2[String,Double] = {
    var retString = ""
    var retProps = Array[Double]()
    consensusValues._1.zipWithIndex.foreach{case(base,index) => {
      if (base != '-') {
        retString += base
        retProps :+= consensusValues._2(index)
      }
    }}
    (retString, retProps.sum.toDouble / retProps.length.toDouble)
}

def outputReads(umi: String,
  readsF: Array[String],
  readsR: Array[String],
  namesF: Array[String],
  namesR: Array[String],
  outputFASTAF: PrintWriter,
  outputFASTAR: PrintWriter,
  outputStats: PrintWriter,
  ref: String,
  primers: List[String],
  sample: String): Tuple2[Int,Int] = { //: MergeResult = {


  // filter down the reads to those that have the primers on each end
  var newReadsF = Array[String]()
  var newNamesF = Array[String]()
  var newReadsR = Array[String]()
  var newNamesR = Array[String]()

  val fwdPairs = readsF.zip(namesF)
  val revPairs = readsR.zip(namesR)

  var matchingFor = 0
  var matchingRev = 0
  var total = 0
  var kept = 0

  fwdPairs.zip(revPairs).foreach{case((fwd,fwdName),(rev,revName)) => { //case(read,name) => {
    if (false) {
      println("ERER")
      println(fwd)
      println(rev)
      println(fwdName)
      println(revName)
      println(primers(0))
      println(reverseComplement(primers(1)))
    }

    if (fwd contains primers(0))
      matchingFor += 1
    if (rev contains reverseComplement(primers(1)))
      matchingRev += 1
    total += 1

    //println(reverseComplement(primers(1)))
    if ((fwd contains primers(0)) && (rev contains reverseComplement(primers(1)))) {
      // println("HIT")
      newReadsF :+= fwd
      newNamesF :+= fwdName

      newReadsR :+= rev
      newNamesR :+= revName
      kept += 1
    }
  }}

  if (newReadsF.size < 2) {
    println("WARN: dropping umi " + umi + " since we have less than two surviving reads. FM = " + matchingFor + " RV = " + matchingRev + " total = " + total)
    return (kept,total)
  }

  // setup a CLUSTAL run and farm it out the machine
  val tmpOutputF = clustal_reads(newReadsF, newNamesF)
  val tmpOutputR = clustal_reads(newReadsR, newNamesR)

  // read in the output from CLUSTAL and see how many reads we have
  val mergerF = Merger(tmpOutputF)
  val mergerR = Merger(tmpOutputR)

  // merge the reads down to a single consensus read
  val consensusReadF = twoPassConsensus(mergerF.nameToRead,mergerF.readLength,0.10)
  val consensusReadR = twoPassConsensus(mergerR.nameToRead,mergerR.readLength,0.10)

  val collapsedConcF = collapseConsensus((consensusReadF._1,consensusReadF._2))
  val collapsedConcR = collapseConsensus((consensusReadR._1,consensusReadR._2))

  // require that the consensus read has kept at least 80% of the original reads
  val highConcensusF = (consensusReadF._4 > 0.8)
  val highConcensusR = (consensusReadR._4 > 0.8)
  val forwardPrimer = (collapsedConcF._1.slice(0,primers(0).length + 5) contains primers(0))
  val reversePrimer = (collapsedConcR._1.slice(0,primers(1).length + 5) contains reverseComplement(primers(1)))
  val usingRead =  highConcensusF && highConcensusR && forwardPrimer && reversePrimer

  var failureReason = ""
  if (!highConcensusF) failureReason += "lowConsensusF;"
  if (!highConcensusR) failureReason += "lowConsensusR;"
  if (!forwardPrimer) failureReason += "forwardPrimerMissing;"
  if (!reversePrimer) failureReason += "reversePrimerMissing;"
  if (failureReason == "") failureReason = "PASS"

  outputStats.write(umi + "\t" + usingRead + "\t" + failureReason + "\t") // write the UMI out and if we're using the read (and if it was rev or fwd reads that ruined it)
  outputStats.write(readsF.size + "\t" + newReadsF.size + "\t" + consensusReadF._3 + "\t" + + readsR.size + "\t" + newReadsR.size + "\t" + consensusReadR._3 + "\t") // output the before and after filtering counts of reads (fwd and rev)

  // twoPassConsensus: returns a tuple-4: the consensus, it's average per-base error rate, number of reads used in the final assembly, and the proportion retained of the total
  outputStats.write(collapsedConcF._2 + "\t" + collapsedConcF._1 + "\t")
  outputStats.write(collapsedConcR._2 + "\t" + collapsedConcR._1 + "\n")



  if ( usingRead ) {
    outputFASTAF.write(">" + sample + "_umi_" + umi + "_reads_used_" + consensusReadF._3 + "_overall_match_prop_" + collapsedConcF._2 + "_read_prop_retained_" + consensusReadF._4 + "\n" + collapsedConcF._1 + "\n")
    outputFASTAR.write(">" + sample + "_umi_" + umi + "_reads_used_" + consensusReadR._3 + "_overall_match_prop_" + collapsedConcR._2 + "_read_prop_retained_" + consensusReadR._4 + "\n" + collapsedConcR._1 + "\n")
  }
  return (kept,total)
}

def clustal_reads(newReads: Array[String], newNames:Array[String] ): File = {
  // setup a CLUSTAL run and farm it out the machine
  val tmp = java.io.File.createTempFile("UMIMerger", ".txt")
  val tmpWriter = new PrintWriter(tmp)
  newReads.zip(newNames).foreach{case(read,name) => tmpWriter.write("> " + name + "\n" + read + "\n")}
  tmpWriter.close()
  val tmpOutput = java.io.File.createTempFile("UMIMerger", ".post.txt")
  val clustoResult = ("clustalo --force --pileup -i " + tmp + " -o " + tmpOutput).!!
  return tmpOutput
}

def discover_cluster(newReads: Array[String], newNames:Array[String] ): File = {
  // setup a CLUSTAL run and farm it out the machine
  val tmp = java.io.File.createTempFile("UMIMerger", ".txt")
  val tmpWriter = new PrintWriter(tmp)
  newReads.zip(newNames).foreach{case(read,name) => tmpWriter.write("> " + name + "\n" + read + "\n")}
  tmpWriter.close()
  val tmpOutput = java.io.File.createTempFile("UMIMerger", ".post.txt")
  val clustoResult = ("clustalo --force --pileup -i " + tmp + " -o " + tmpOutput).!!
  return tmpOutput
}


// *********************************************************************************************************
// print the usage parameters and quit
def printUsageAndQuit() {
  println("usage:\n\nscala umi_count.scala <reads1.fq> <reads2.fq> <output_fasta_forward_reads> <output_fasta_reverse_reads> <output_stats> <umi_cutoff> <primers_file> <reference_file> <sample.name>")
  sys.exit(1)
}


// *********************************************************************************************************
//
//
// main script point
//
//
// *********************************************************************************************************

// check that the correct number of inputs was provided
if (args.length != 9) printUsageAndQuit()
val inputFileReads1 = args(0)
val inputFileReads2 = args(1)
val outputFileF = new PrintWriter(args(2))
val outputFileR = new PrintWriter(args(3))
val outputStats = new PrintWriter(args(4))
val umiLength = args(5).toInt
val primersEachEnd = args(6)
val reference = args(7)
val samplename= args(8)

// get the reference string
var referenceString = ""
Source.fromFile(reference).getLines().foreach{line => if (!line.startsWith(">")) referenceString += line}

// first find out what UMIs to keep
// ------------------------------------------------------------------------------------------
val forwardReads = Source.fromInputStream(gis(inputFileReads1)).getLines().grouped(4)
val reverseReads = Source.fromInputStream(gis(inputFileReads2)).getLines().grouped(4)

val primers = Source.fromFile(primersEachEnd).getLines().map{line => line}.toList
if (primers.length != 2)
  throw new IllegalStateException("You should only provide a primer file with two primers")

// our containers for forward and reverse reads
val umiReadsFWD = new HashMap[String,ArrayBuilder[String]]()
val umiReadNamesFWD = new HashMap[String,ArrayBuilder[String]]()
val umiReadsRVS = new HashMap[String,ArrayBuilder[String]]()
val umiReadNamesRVS = new HashMap[String,ArrayBuilder[String]]()


// find the UMIs in the forward reads
var readsProcessed = 0
forwardReads foreach {fGroup => {
  val rGroup = reverseReads.next()

  val umi = fGroup(1).slice(0,umiLength)
  val readNoUMI = fGroup(1).slice(umiLength,fGroup(1).length)
  val qualNoUMI = fGroup(3).slice(umiLength,fGroup(3).length)

  val readBuilderF = umiReadsFWD.getOrElse(umi,ArrayBuilder.make[String])
  val nameBuilderF = umiReadNamesFWD.getOrElse(umi,ArrayBuilder.make[String])

  val readBuilderR = umiReadsRVS.getOrElse(umi,ArrayBuilder.make[String])
  val nameBuilderR = umiReadNamesRVS.getOrElse(umi,ArrayBuilder.make[String])

  readBuilderF += qualityControlRead(readNoUMI,qualNoUMI)
  nameBuilderF += fGroup(0)
  umiReadsFWD(umi) = readBuilderF
  umiReadNamesFWD(umi) = nameBuilderF

  readBuilderR += qualityControlRead(rGroup(1),rGroup(3))
  nameBuilderR += rGroup(0)
  umiReadsRVS(umi) = readBuilderR
  umiReadNamesRVS(umi) = nameBuilderR

  readsProcessed += 1
}}

outputStats.write("UMI\tused\tffail.reason\tfwd.reads\tprimered.fwd.reads\tfiltered.fwd.reads\trev.reads\tprimered.rev.reads\tfiltered.rev.reads\t")
outputStats.write("fwd.error.rate\tfwd.read\trev.error.rate\trev.read\n")

var passingUMI = 0
var totalWithUMI = 0
var usedReads = 0

println("Processing " + umiReadsFWD.size + " UMIs")
var index = 0
umiReadsFWD.foreach{ case(umi,reads) => {
  val reverseReads = umiReadsRVS(umi).result()
  val forwardReads = reads.result()

  val nameArrayF = umiReadNamesFWD(umi).result()
  val nameArrayR = umiReadNamesRVS(umi).result()

  if (reverseReads.length > 30) {
    val res = outputReads(umi,forwardReads, reverseReads, nameArrayF, nameArrayR, outputFileF, outputFileR, outputStats, referenceString,primers, samplename)
    totalWithUMI += res._2
    usedReads += res._1
    passingUMI += 1
  }
  index += 1
  if (index % 100 == 0) {
    println("INFO: Processed " + index + " umis so far")
  }
}}

println("Summary:\ntotal reads:\t"+ readsProcessed + "\ntotal reads withing a valid UMI:\t" + totalWithUMI + "\ntotal reads used for consensus reads:\t" + usedReads + "\ntotal passing UMIs:\t" + passingUMI)

outputStats.close()
outputFileF.close()
outputFileR.close()
