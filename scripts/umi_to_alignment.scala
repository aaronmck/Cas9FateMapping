// given a stats file and a UMI, split out the read pair in question, align with MAFFT,
// and pretty print the output

import scala.io._
import scala.sys.process._
import java.io._

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

// target 3
//val refSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTAATGATACGGCGACCACCGAGATCTACACtctaacgtctTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNATCTCGAGCTCAAGCTTCGGATACGATACGCGCACGCTATGGAGTCGACACGACTCGCGCATACGATGGAGTCGATAGTATGCGTATACGCTATGGAGTCGATATGCATAGCGCATGCTATGGAGTCGAGTCGAGACGCTGACGATATGGAGTCGCTACGATACACTCTGACTATGGAGTCGCGACTGTACGCACACGCGATGGAGTCGATACGTAGCACGCAGACTATGGAGTCGACACAGTACTCTCACTCTATGGAGTCGATATGAGACTCGCATGTGATGGGAATTCTCGACCTCGAGACAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACAACGACGTCCATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAA"

// target 1
val refSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTAATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNTTCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNATCTCGAGCTCAAGCTTCGGACAGCAGTATCATCGACTATGGAGTCGAGAGCGCGCTCGTCGACTATGGAGTCGTCAGCAGTACTACTGACGATGGAGTCGACAGCAGTGTGTGAGTCTATGGAGTCGAGAGCATAGACATCGAGTATGGAGTCGACTACAGTCGCTACGACTATGGAGTCGACAGAGATATCATGCAGTATGGAGTCGACAGCAGTATCTGCTGTCATGGAGTCGACTGCACGACAGTCGACTATGGAGTCGCGAGCGCTATGAGCGACTATGGGAATTCTCGACCTCGAGACAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNNNATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAAAAAAAAAAAA"

val umi = args(1)
val lines = Source.fromFile(args(0)).getLines()
val header = lines.next.split("\t").zipWithIndex.map{case(tk,ind) => (tk,ind)}.toMap
val statsLine = lines.filter{line => line startsWith umi}.toArray

if (statsLine.size != 1) {
  println("Failed to find one and only one umi line for " + umi + " found " + statsLine.size)
  sys.exit(1)
}

val read1 = statsLine(0).split("\t")(header("readF"))
val read2 = reverseComplement(statsLine(0).split("\t")(header("readR")))

val outRef = "ref.fasta"
val outRefF = new PrintWriter(outRef)
outRefF.write(">ref\n" + refSeq.toUpperCase + "\n")
outRefF.close()

val outRead1 = "read1.fasta"
val outFRead1File = new PrintWriter(outRead1)
outFRead1File.write(">read1\n" + read1.toUpperCase + "\n")
outFRead1File.close()

val outRead2 = "read2.fasta"
val outFRead2File = new PrintWriter(outRead2)
outFRead2File.write(">read2\n" + read2.toUpperCase + "\n")
outFRead2File.close()

val outRead12 = "read12.fasta"
val outFRead12File = new PrintWriter(outRead12)
outFRead12File.write(">read1\n" + read1.toUpperCase + "\n")
outFRead12File.write(">read2\n" + read2.toUpperCase + "\n")
outFRead12File.close()

val outFl = "temp.fasta"
val output = new PrintWriter(outFl)
output.write(">ref\n" + refSeq.toUpperCase + "\n")
output.write(">read1\n" + read1.toUpperCase + "\n")
output.write(">read2\n" + read2.toUpperCase + "\n")
output.close()

val outFlF = "tempF.fasta"
val outputF = new PrintWriter(outFlF)
outputF.write(">ref\n" + refSeq.toUpperCase + "\n")
outputF.write(">read1\n" + read1.toUpperCase + "\n")
outputF.close()

val outFlR = "tempR.fasta"
val outputR = new PrintWriter(outFlR)
outputR.write(">ref\n" + refSeq.toUpperCase + "\n")
outputR.write(">read2\n" + read2.toUpperCase + "\n")
outputR.close()

def align(readFile: String, refFile:String, aligner: String): Array[String] = {

  val out = new StringBuilder
  val err = new StringBuilder
  val logger = ProcessLogger((o: String) => out.append(o + "\n"),(e: String) => err.append(e + "\n"))

  aligner match {
    case "mafft" => ("mafft --maxiterate 1000 --genafpair " + readFile) ! logger
    case "clustal" => ("clustalo -i " + readFile) ! logger
    case "water" => ("water -datafile /net/shendure/vol10/projects/CRISPR.lineage/nobackup/reference/EDNAFULL -asequence " + readFile + " -bsequence " + refFile + " -gapopen 10 -gapextend 0.1 -stdout -auto") ! logger
  }

//println(("water -asequence " + readFile + " -bsequence " + refFile + " -gapopen 10 -gapextend 0.5 -stdout"))
  var readNames = Array[String]()
  var readStrings = Array[String]()
  var currentRead = ""

  out.toString().split("\n") foreach { line =>
        if (line startsWith ">") {
          if (currentRead != "")
            readStrings :+= currentRead
          currentRead = ""
          readNames :+= line.stripPrefix(">")
        } else {
          currentRead += line.toUpperCase
        }
  }
  if (currentRead != "")
    readStrings :+= currentRead

  var refF = ""
  var read1F = ""
  var read2F = ""

  //println("OUT " + out)
  //println("ERR " + err)
  //println("READ NAMES " + readNames.mkString("-"))
  (0 until readStrings(0).size).foreach{i =>
    if (readStrings(0)(i) == readStrings(1)(i) && readStrings(0)(i) == readStrings(2)(i)) {
      refF += readStrings(0)(i).toString.toUpperCase
    }
    else {
      refF += readStrings(0)(i).toString.toLowerCase
    }

    if (readStrings(0)(i) == readStrings(1)(i)) {
      read1F += readStrings(1)(i).toString.toUpperCase
    } else {
      read1F += readStrings(1)(i).toString.toLowerCase
    }

    if (readStrings(0)(i) == readStrings(2)(i)) {
      read2F += readStrings(2)(i).toString.toUpperCase
    } else {
      read2F += readStrings(2)(i).toString.toLowerCase
    }
  }
  return Array[String](refF,read1F,read2F)
}


def alignSingle(readFile: String, refFile:String, aligner: String): Array[String] = {

  val out = new StringBuilder
  val err = new StringBuilder
  val logger = ProcessLogger((o: String) => out.append(o + "\n"),(e: String) => err.append(e + "\n"))

  aligner match {
    case "mafft" => ("mafft --maxiterate 1000 --genafpair " + readFile) ! logger
    case "clustal" => ("clustalo -i " + readFile) ! logger
    case "water" => ("water -datafile /net/shendure/vol10/projects/CRISPR.lineage/nobackup/reference/EDNAFULL -asequence " + readFile + " -bsequence " + refFile + " -gapopen 10 -gapextend 0.1 -stdout -auto") ! logger
  }

  var readNames = Array[String]()
  var readStrings = Array[String]()
  var currentRead = ""

  out.toString().split("\n") foreach { line =>
        if (line startsWith ">") {
          if (currentRead != "")
            readStrings :+= currentRead
          currentRead = ""
          readNames :+= line.stripPrefix(">")
        } else {
          currentRead += line.toUpperCase
        }
  }
  if (currentRead != "")
    readStrings :+= currentRead

  var refF = ""
  var read1F = ""
  var read2F = ""

  (0 until readStrings(0).size).foreach{i =>
    if (readStrings(0)(i) == readStrings(1)(i)) {
      refF += readStrings(0)(i).toString.toUpperCase
      read1F += readStrings(1)(i).toString.toUpperCase
    }
    else {
      refF += readStrings(0)(i).toString.toLowerCase
      read1F += readStrings(1)(i).toString.toLowerCase
    }
  }
  return Array[String](refF,read1F,read2F)
}

println("MAFFT\n" + align(outFl,outFl,"mafft").mkString("\n") + "\nFWD:\n" + alignSingle(outFlF,outFlF,"mafft").mkString("\n") + "\nREV:\n" + alignSingle(outFlR,outFlR,"mafft").mkString("\n"))

println("CLUSTALO\n" + align(outFl,outFl,"clustal").mkString("\n") + "\nFWD:\n" + alignSingle(outFlF,outFlF,"clustal").mkString("\n") + "\nREV:\n" + alignSingle(outFlR,outFlR,"clustal").mkString("\n"))

//println("WATER\n" + align(outRead12,outRef,"water").mkString("\n") + "\nFWD:\n" + alignSingle(outRead1,outRef,"water").mkString("\n") + "\nREV:\n" + alignSingle(outRead2,outRef,"water").mkString("\n"))
