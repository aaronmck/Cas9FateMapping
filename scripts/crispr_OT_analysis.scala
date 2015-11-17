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

// add a read's cigar string, adding any indels we see, returning
// a tuple, the first value is true if there were any events and an
// array of events at each position
// --------------------------------------------------------------------------------
def mapToSite(cigar: String, cutsite: Tuple3[Int, Int, Int]): Tuple2[Boolean, Array[String]] = {
  var currentPos = 0

  var ret = Array[String]()

  (new Regex("(\\d+)([MIDSNHP=X])") findAllIn cigar).matchData.foreach { matching =>
    {
      matching.subgroups(1) match {

        // matching -- just update our read position
        case "M" => { currentPos += matching.subgroups(0).toInt }

        // insertion -- check that we are close enough to the cutsite
        case "I" => {
          val len = matching.subgroups(0).toInt
          if ((currentPos < cutsite._2 && len + currentPos > cutsite._1) || (currentPos >= cutsite._2 && currentPos < cutsite._3)) {
            ret :+= len + "I-" + currentPos
          }
          currentPos += len
        }

        case "D" => {
          val len = matching.subgroups(0).toInt
          if ((currentPos < cutsite._2 && len + currentPos > cutsite._1) || (currentPos >= cutsite._2 && currentPos < cutsite._3)) {
            ret :+= len + "D-" + currentPos
          }
          currentPos += len
        }
        case "S" => {
          // do nothing
        }
        case "H" => {
          // do nothing
        }
      }
    }
  }

  return Tuple2[Boolean, Array[String]](ret.mkString("") != "", ret)
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

val umiStart = 177
val umiStop = 189

val editStart = 142
val editEnd = 153
// command line args:
// 1) control SAM file
// 2) edit SAM file
// 3) reference fasta file
// 4) output file

// check that the correct number of inputs was provided
val referenceString = Source.fromFile(args(2)).getLines().drop(1).mkString("")


println("loading reads for control UMI counts...")
val umiToCountControl = readsToUMICounts(Source.fromFile(args(0)).getLines(),umiStart,umiStop,referenceString,editStart,editEnd)
println("loading reads for case UMI counts...")
val umiToCountEdit = readsToUMICounts(Source.fromFile(args(1)).getLines(),umiStart,umiStop,referenceString,editStart,editEnd)
val umiOutput = new PrintWriter(args(3))
println("setting up output...")
/*
val ref = "GATCACTCTCGGCATGGACGAGCTGTACAAGTCCGGACTCAGATCTCGAGCTCAAGCTTCGAATTCTCGACCTCGAGACAAATGCGAGCTCTGACTACATCGCTACATTCAGCTCGATCACAGCTTACATCGCTTCTGATCANNNNNNATCAGTCGGATGCGATGCTGATACTGTATNNNNNNNNNNNNCTGTTAGCCATAGCAGTCGTGCGATGTAGCAGCCTGAGACAAATGGCAGTATTCATCCAC"
val bses = "GATCACTCTCGGCATGGACGAGCTGTACAAGTCCGGACTCAGATCTCGAGCTCAAGCTTCGAATTCTCGACCTCGAGACAAATGCGAGCTCTGACTACATCGCTACATTCAGCTCGATCACAGCTTACATCGCTGCTGGTCGCGCGAATATAAGTCGGATGCGATGCTGATACTGTATCCATTTAGTGCACTGTTAGCCATAGCAGTCGTGCGATGTAGCAGCCTGAGACAAATGGCAGTATTCATCCAC"
println(zipCigar(bses, ref, "145M1I104M").map{tk => tk._4 + "-" + tk._5 }.mkString(","))
println(bses)
println(zipCigar(bses, ref, "145M1I104M").map{tk => if ((tk._4 > editStart && tk._4 < editEnd) && (tk._5 == 'D' || tk._5 == 'I')) '_' else tk._5}.mkString)
*/

umiOutput.write("UMI\tnumberInControl\teditedInControls\tnumberInCase\teditedInCase\t" + ((1 until 13).map{i => "umiBase" + i}.mkString("\t")) + "\n")
(umiToCountControl.keySet ++ umiToCountEdit.keySet).foreach{umi => {
  val controlCount = umiToCountControl.getOrElse(umi,(0,0))
  val editCount = umiToCountEdit.getOrElse(umi,(0,0))
  umiOutput.write(umi + "\t" + controlCount._1 + "\t" + controlCount._2 + "\t" + editCount._1 + "\t" + editCount._2 + "\t" + umi.map{t => t}.slice(0,12).mkString("\t") + "\n")
}}
umiOutput.close()
