package aligner
package aligner

import java.io._

import _root_.aligner.MAFFT._
import main.scala.utils.Utils

import scala.collection.mutable._
import scala.io._
import scala.main.{ReadDirection, SequencingRead}
import scala.sys.process._

/**
 * Created by aaronmck on 10/22/15.
 */
// ------------------------------------------------------------------------------------------
//
// process reads to and from the clustalo algorithm
//
// ------------------------------------------------------------------------------------------
object Waterman extends Aligner {
  val referenceName = "waterRef"


  override def alignTo(reads: Array[SequencingRead], ref: Option[String], debug: Boolean): Array[SequencingRead] = {
    alignTo(reads, ref, debug, 10.0, 0.5)
  }

  def alignTo(reads: Array[SequencingRead], ref: Option[String], debug: Boolean, gapOpen: Double, gapExtend: Double): Array[SequencingRead] = {
    if (reads.size != 1 && ref.isDefined)
      throw new IllegalStateException("Unable to run waterman with " + reads.size + " reads with a reference")
    if (reads.size != 2 && !ref.isDefined)
      throw new IllegalStateException("Unable to run waterman with " + reads.size + " reads without a reference")

    val tmp = java.io.File.createTempFile("UMIMergerBWA", ".fasta")
    val tmp2 = java.io.File.createTempFile("UMIMergerBWA", ".fasta")
    //println(tmp + " " + tmp2)
    val tmpWriter = new PrintWriter(tmp)
    val tmpWriter2 = new PrintWriter(tmp2)

    if (ref.isDefined) {
      tmpWriter.write(">" + referenceName + "\n" + ref.get + "\n")
      tmpWriter.close()

      if (reads(0).reverseCompAlign)
        tmpWriter2.write(">" + reads(0).name + "\n" + Utils.reverseComplement(reads(0).bases.filter { bs => bs != '-' }.mkString("")) + "\n")
      else
        tmpWriter2.write(">" + reads(0).name + "\n" + reads(0).bases.filter { bs => bs != '-' }.mkString("") + "\n")
      tmpWriter2.close()
    } else {
      //println(reads(0).reverseCompAlign)
      //println(reads(1).reverseCompAlign)

      if (reads(0).reverseCompAlign)
        tmpWriter.write(">" + reads(0).name + "\n" + Utils.reverseComplement(reads(0).bases.filter { bs => bs != '-' }.mkString("")) + "\n")
      else
        tmpWriter.write(">" + reads(0).name + "\n" + reads(0).bases.filter { bs => bs != '-' }.mkString("") + "\n")
      tmpWriter.close()


      if (reads(1).reverseCompAlign)
        tmpWriter2.write(">" + reads(1).name + "\n" + Utils.reverseComplement(reads(1).bases.filter { bs => bs != '-' }.mkString("")) + "\n")
      else
        tmpWriter2.write(">" + reads(1).name + "\n" + reads(1).bases.filter { bs => bs != '-' }.mkString("") + "\n")
      tmpWriter2.close()

    }

    /**
     * Run needlall and capture the output
     */
    val out = new StringBuilder
    val err = new StringBuilder
    val logger = ProcessLogger(
      (o: String) => out.append(o + "\n"),
      (e: String) => err.append(e + "\n"))

    val call = "needleall -datafile /net/shendure/vol10/projects/CRISPR.lineage/nobackup/reference_data/EDNAFULL -snucleotide1 -snucleotide2 -gapopen " +
      gapOpen + " -gapextend " + gapExtend + " -asequence " + tmp + " -bsequence " +
      tmp2 + " -aformat3 fasta -auto -stdout"

    (call) ! logger

    var readNames = Array[String]()
    var readStrings = Array[String]()
    var currentRead = ""


    out.toString().split("\n") foreach { line =>
      //println(line)
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

    // make an array of SequenceReads to store the result in
    var ret = Array[SequencingRead]()

    if (readStrings.size != 2)
      return Array[SequencingRead]()

    if (ref.isDefined) {
      ret :+= Aligner.SequencingReadFromNameBases(referenceName, readStrings(0))
      ret :+= Aligner.SequencingReadFromNameBases(readNames(1), readStrings(1))
      //println("reference: " + readStrings(0))
      //println("read     : " + readStrings(1))
    } else {
      ret :+= Aligner.SequencingReadFromNameBases(readNames(0), readStrings(0))
      ret :+= Aligner.SequencingReadFromNameBases(readNames(1), readStrings(1))
      //println("read     : " + readStrings(0))
    }

    tmp.delete()
    tmp2.delete()
    return ret
  }
}
