package aligner

import java.io._

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
object Clustalo extends Aligner {
  /**
   * align two sequences
   * @param reads a list of reads, which we pull sequence out of
   * @param ref the reference sequence
   * @param debug dump a lot of info
   * @return a array of aligned sequencing reads, the reference will be the first
   */
  def alignTo(reads: Array[SequencingRead], ref: Option[String], debug: Boolean = false): Array[SequencingRead] = {
    val tmp = java.io.File.createTempFile("UMIMergerBWA", ".txt")
    val tmpWriter = new PrintWriter(tmp)
    var readDirections = Array[ReadDirection]()

    // write the reads / reference to the input file
    if (ref.isDefined)
      tmpWriter.write(">reference\n" + ref.get + "\n")
    reads.foreach {
      rd =>
        if (rd.reverseCompAlign)
          tmpWriter.write(">" + rd.name + "\n" + Utils.reverseComplement(rd.bases.filter {
            bs => bs != '-'
          }.mkString("")) + "\n")
        else
          tmpWriter.write(">" + rd.name + "\n" + rd.bases.filter {
            bs => bs != '-'
          }.mkString("") + "\n")
    }
    tmpWriter.close()

    // make an array of SequenceReads to store the result in
    var ret = Array[SequencingRead]()

    /**
     * Run MAFFT and capture the output
     */
    val out = new StringBuilder
    val err = new StringBuilder
    val logger = ProcessLogger(
      (o: String) => out.append(o + "\n"),
      (e: String) => err.append(e + "\n"))

    if (debug)
      println("clustalo --dealign --force --pileup -i " + tmp)
    ("clustalo --dealign --force --pileup -i " + tmp) ! logger

    /**
     * now readback each read figure out the alignment
     */
    var readNames = Array[String]()
    var readStrings = Array[String]()
    var currentRead = ""

    out.toString().split("\n") foreach {
      line =>
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

    if ((ref.isDefined && readStrings.length != reads.size + 1) || (!ref.isDefined && readStrings.length != reads.size)) {
      throw new IllegalStateException("DIDNT get all our reads back: " + readStrings.length + " instead of " + reads.size)
    }
    for (i <- 0 until readStrings.length) {
      if (debug)
        println(readStrings(i))
      ret :+= Aligner.SequencingReadFromNameBases(readNames(i), readStrings(i))
    }

    //tmp.delete()
    return ret
  }

}
