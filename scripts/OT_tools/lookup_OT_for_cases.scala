import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._
import scala.util.matching.Regex

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
def reverseComplement(str: String) = str.map { t => compBase(t) }.reverse.mkString("")

// *********************************************************************************************************
// print the usage parameters and quit
def printUsageAndQuit() {
  println("usage:\n\nscala umi_count.scala <reads1.fq> <reads2.fq> <output_fasta_forward_reads> <output_fasta_reverse_reads> <output_stats> <umi_cutoff> <primers_file> <reference_file> <sample.name>")
  sys.exit(1)
}

def readsToUMICounts(lines: Iterator[String], umiStart: Int, umiStop: Int, refString: String, editStart: Int, editStop: Int): HashMap[String, Tuple2[Int,Int]] =  {
  val umiToCount = new HashMap[String, Tuple2[Int,Int]]()
  var cnt = 0
  lines.foreach { line =>
    {
      if (!line.startsWith("@")) {
        val sp = line.split("\t")
        val bases = sp(9)
        val cigar = sp(5)

        val zippedToRef = zipCigar(bases, refString, cigar)

        val edited = if (zippedToRef.map{tk => if ((tk._4 >= editStart && tk._4 <= editStop) && (tk._5 == 'I' || tk._5 == 'D')) 1 else 0}.sum > 0) 1 else 0
        val umi = zippedToRef.filter(zp => zp._4 >= umiStart && zp._4 < umiStop).map{zp => zp._1}.mkString("")
        var oldValue = umiToCount.getOrElse(umi,(0,0))
        umiToCount(umi) = (oldValue._1 + 1,oldValue._2 + edited)

        cnt += 1
        if (cnt % 10000 == 0) {
          println("Count = " + cnt)
        }
      }
    }
  }
  umiToCount
}

  // zip by the cigar string, returning a tuple of:
  // 1: read base
  // 2: reference base
  // 3: do they match? no if an indel
  // 4: position in the reference
  // 5: cigar character
  def zipCigar(baseString: String, referenceString: String, cigarString: String): Array[Tuple5[Char, Char, Boolean, Int, Char]] = {
    var returnTuples = Array[Tuple5[Char, Char, Boolean, Int, Char]]()
    var baseIndex = 0
    var refIndex = 0


    (new Regex("(\\d+)([MIDSNHP=X])") findAllIn cigarString).matchData.foreach { matching => {
      //println("HERE " + matching.subgroups(1) + " for " + (matching.subgroups(0).toInt))
      try {
      (0 until matching.subgroups(0).toInt).foreach(ignored => {
        matching.subgroups(1) match {
          case "M" => {
            try {
              returnTuples :+=(baseString(baseIndex), referenceString(refIndex), baseString(baseIndex) == referenceString(refIndex), refIndex, 'M')
            } catch {
              case e: Exception => {
                println("base index = " + baseIndex + " refIndex = " + " " + cigarString); throw e
              }
            }
            //println("M + " + ignored + " = " + baseIndex + " refindex = " + refIndex)
            baseIndex += 1
            refIndex += 1
          }
          case "I" => {
            returnTuples :+=(baseString(baseIndex), '-', false, refIndex, 'I')
            //println("I + " + ignored + " = " + baseIndex + " refindex = " + refIndex)
            baseIndex += 1
          }
          case "D" => {
            returnTuples :+=('-', referenceString(refIndex), false, refIndex, 'D')
            //println("D + " + ignored + " = " + baseIndex + " refindex = " + refIndex)
            refIndex += 1
          }
          case "S" => {
            returnTuples :+=('+', '+', false, refIndex, 'S')
            //println("S + " + ignored + " = " + baseIndex + " refindex = " + refIndex)
            baseIndex += 1
          }
          case _ => {
            println("Uknown cigar")
          }
        }
      })
    } catch {
      case e: java.lang.StringIndexOutOfBoundsException => {
        println("DROPPED: failed on cigar " + cigarString + " and base string " + baseString)
        return returnTuples
      }
    }
    }
    }
    returnTuples
  }


def ployTrack(str: String): Int = {
    var lastLetter = 'N'
    var currentRun = 0
    var largestRun = 0
    str.foreach{letter => {
      if (letter == lastLetter && letter != "-") {
        currentRun += 1
        if (currentRun > largestRun) largestRun = currentRun
      } else {
        lastLetter = letter
        currentRun = 1
      }
    }}
    largestRun
}
// *********************************************************************************************************
//
//
// main script point
//
//
// *********************************************************************************************************

// a mapping from umi to insert / degenerate site
val umiToReference = new HashMap[String,ArrayBuilder[String]]

val reference = "GATCACTCTCGGCATGGACGAGCTGTACAAGTCCGGACTCAGATCTCGAGCTCAAGCTTCGAATTCTCGACCTCGAGACAAATGCGAGCTCTGACTACATCGCTACATTCAGCTCGATCACAGCTTACATCGCTTCTGATCACGCGATATCAGTCGGATGCGATGCTGATACTGTATNNNNNNNNNNNNCTGTTAGCCATAGCAGTCGTGCGATGTAGCAGCCTGAGACAAATGGCAGTATTCATCCAC"

val insertStart = 134
val insertStop =  154
val umiStart =    177
val umiStop =     189

// get the mapping from experimental condition to the umi -> target sequence for those that pass the filter
val umiToSeq = new HashMap[Int,HashMap[String,String]]()
var dropped = 0
var collisions = 0
var same = 0

Source.fromFile("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_31_OT_Deep_Sequencing/data/umi_to_reference.txt").getLines().drop(1).foreach{str => {
  val sp = str.split("\t")
  val sample = sp(0)
  val umi = sp(1)
  val consensus = sp(2)
  val cnt = sp(4).toInt

  val day = sample.slice(6,7).toInt
  val experiment = sample.slice(16,17).toInt

  val umiToSeqExp = umiToSeq.getOrElse(experiment,HashMap[String,String]())


  if (umiToSeqExp contains umi) {
    //println("Collision detected with UMI " + umi + " old seq " + umiToSeqExp(umi) + " new seq " + consensus + " same? " + (umiToSeqExp(umi) == consensus))
    if (umiToSeqExp(umi) == consensus) same += 1
    collisions += 1
  } else if (!(sp(0) contains "-") && cnt >= 5)
    umiToSeqExp(umi) = consensus
  else
    dropped += 1

  umiToSeq(experiment) =  umiToSeqExp
}}
println("Dropped " + dropped + " control samples due to low read backing or insertions")
println("Collisions " + collisions + " in UMI assignments, same sequence count = " + same)

// setup the output files
val output = new PrintWriter("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_31_OT_Deep_Sequencing/data/edit_sites.txt")
output.write("file\tumi\tinsert\tcigar\tedited\tcontrol.seq\tcontrol.edited\tmatched\tcontains.ns\tcontains.poly\n")




// *********************************************************************************************************
// read in the Cas9 exposed samples, find if they have a matched control, figure out how probable it is that
// the reads map back to the given reference sequence
// *********************************************************************************************************
Source.fromFile("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_31_OT_Deep_Sequencing/data/case_files_sam.txt").getLines().foreach{caseFile => {
  println("Case file " + caseFile)
  val name = (new File(caseFile)).getName.stripSuffix(".fixed.sam")

  var included = 0
  var dropped = 0

  val experiment =
    if(name startsWith "OT_day")
      name.slice(6,7).toInt
    else
      name.slice(11,12).toInt

  val comparisonUMIs = umiToSeq(experiment)

  Source.fromFile(caseFile).getLines.foreach{line => {
    if (!(line startsWith "@")) {
      val sp = line.split("\t")
      val cigar = sp(5)
      val read = sp(9)

      // 1: read base
      // 2: reference base
      // 3: do they match? no if an indel
      // 4: position in the reference
      // 5: cigar character
      val zippedCigar = zipCigar(read, reference, cigar)
      var umi = ""
      var insert = ""

      zippedCigar.foreach{case(readBase,refBase,matched,refPos,cigar) => {
        if (refPos >= umiStart && refPos < umiStop)
          umi += readBase.toString
        if (refPos >= insertStart && refPos < insertStop)
          insert += readBase.toString
      }}

      if (comparisonUMIs contains umi) {
        // check a bunch of conditions for each read site / UMI combination.  Make the output easy to grep from (don't resuse TRUE / FALSE, etc)
        val edited = insert.slice(10,18).contains("-") || insert.size > 20 // do we have either an insertion or a deletion?
        val editedWT = if (comparisonUMIs(umi).slice(10,18).contains("-") || comparisonUMIs(umi).size > 20) "REFEDIT" else "WT" // what about the wild type? !! We're not catching insertions here as we align to the reference, catch them below
        val matched = if (insert.slice(0,20) == comparisonUMIs(umi).slice(0,20)) "MATCHED" else "DIFF" // do the first 20 bases of the sequences match?
        val containsNs = if ((insert.slice(0,20) contains "N") || (comparisonUMIs(umi).slice(0,20) contains "N")) "YES" else "NO" // do either contain Ns
        val containsPoly = if (ployTrack(insert) > 4 || ployTrack(comparisonUMIs(umi).slice(0,20)) > 4) "POLY" else "NP"
        output.write(name + "\t" + umi + "\t" + insert + "\t" + cigar + "\t" + edited + "\t" + comparisonUMIs(umi) + "\t" + editedWT + "\t" + matched + "\t" + containsNs + "\t" + containsPoly + "\n")
        included += 1
      } else {
        dropped += 1
      }
    }
  }}

  println("Included " + included + " dropped " + dropped)
}}

output.close()
