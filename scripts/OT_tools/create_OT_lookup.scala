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
// *********************************************************************************************************
//
//
// main script point
//
//
// *********************************************************************************************************

// a mapping from umi to insert / degenerate site
val reference = "GATCACTCTCGGCATGGACGAGCTGTACAAGTCCGGACTCAGATCTCGAGCTCAAGCTTCGAATTCTCGACCTCGAGACAAATGCGAGCTCTGACTACATCGCTACATTCAGCTCGATCACAGCTTACATCGCTTCTGATCACGCGATATCAGTCGGATGCGATGCTGATACTGTATNNNNNNNNNNNNCTGTTAGCCATAGCAGTCGTGCGATGTAGCAGCCTGAGACAAATGGCAGTATTCATCCAC"

val insertStart = 134
val insertStop =  154
val umiStart =    177
val umiStop =     189

val output = new PrintWriter("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_31_OT_Deep_Sequencing/data/umi_to_reference.txt")
output.write("sample\tumi\tconsensus\tmatching\treads\n")

// for each input file, zip each read to reference_file
Source.fromFile("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_31_OT_Deep_Sequencing/data/control_files_sam.txt").getLines().foreach{controlFile => {
  println("Control file " + controlFile)
  val name = (new File(controlFile)).getName.stripSuffix(".fixed.sam")
  val umiToReference = new HashMap[String,ArrayBuilder[String]]()

  Source.fromFile(controlFile).getLines.foreach{line => {
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

      umiToReference(umi) = umiToReference.getOrElse(umi,ArrayBuilder.make[String]()) += insert
    }
  }}

  umiToReference.foreach{case(umi,inserts) => {
    val insertSize = inserts.result.map{ins => ins.length}.toArray.min
    var consensus = ""
    var basesProb = 0.0
    for (i <- 0 until insertSize) {
      val conc = new HashMap[Char,Int]()
      inserts.result.foreach{ins => conc(ins(i)) = conc.getOrElse(ins(i),0) + 1}
      val maxValue = conc.valuesIterator.max.toDouble
      val sumValue = conc.valuesIterator.sum.toDouble
      val matching = maxValue / sumValue
      basesProb += matching

      if (matching > 0.90)
        consensus += (conc.retain{case(a,b) => b == maxValue}).head._1
      else
        consensus += "N"
    }
    // println(insertSize + " " + consensus)
    output.write(name + "\t" + umi + "\t" + consensus + "\t" + (basesProb / insertSize.toDouble) + "\t" + inserts.result.size + "\n")
  }}
}}
  output.close()
