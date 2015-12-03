package aligner
import java.io.{PrintWriter, File}

import main.scala.utils.Utils

import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.main.{ReadDirection, ReverseReadOrientation, ForwardReadOrientation, SequencingRead}
import scala.sys.process._

/**
 * Created by aaronmck on 11/17/15.
 */
case class SmithWaterman(reads: Array[SequencingRead], referenceFile: File) { //extends Aligner {

  // ----------------------------------------------------------------------
  // setup a BWA run and farm it out the machine
  // ----------------------------------------------------------------------
  val tmp = java.io.File.createTempFile("UMIMergerBWA", ".txt")
  val tmpWriter = new PrintWriter(tmp)
  var readDirections = Array[ReadDirection]()

  println("read len " + reads.size)
  reads.foreach {
    case (read) => {
      readDirections :+= read.readOrientation
      println(read.readOrientation)
      val readBases = read.bases.filter { bs => bs != '-' }
      val readQuals = read.bases.zip(read.quals).filter { case(bs,ql) => bs != '-' }.map{case(bs,ql) => ql}.mkString("")
      tmpWriter.write("@" + read.name.replace(" ","-") + "\n" + readBases + "\n+\n" + readQuals + "\n")
      //tmpWriter.write(read.name + "\n" + readBases + "\n+\n" + readQuals + "\n")
    }
  }
  tmpWriter.close()


  /**
   * we use this code to capture the output from the BWA run
   */
  val out = new StringBuilder
  val err = new StringBuilder
  val logger = ProcessLogger(
    (o: String) => out.append(o + "\n"),
    (e: String) => err.append(e + "\n"))

  // /net/gs/vol1/home/aaronmck/tools/smith_waterman/src/./ssw_test -r -m 6 -x 4 -o 6 -e 1 -c -s


  println("/net/gs/vol1/home/aaronmck/tools/smith_waterman/src/./ssw_test -r -m 10 -x 90 -o 60 -e 1 -c -s " + referenceFile.getAbsolutePath + " " + tmp)
  ("/net/gs/vol1/home/aaronmck/tools/smith_waterman/src/./ssw_test -r -m 10 -x 90 -o 60 -e 1 -c -s " + referenceFile.getAbsolutePath + " " + tmp) ! logger

  var returnReads = Array[SequencingRead]()

  var readID = 0
  out.toString().split("\n")foreach{ line => {
    if (!(line startsWith "@")) {
      println("LINE:" + line + " " + readID )

      val splitRead = line.split("\t")
      //val readNameSplit = splitRead(0).split("&&")
      // val (name, umi) = (readNameSplit(0), readNameSplit(1))
      val name = splitRead(0)
      returnReads :+= SequencingRead(name, splitRead(9), splitRead(10),readDirections(readID) , "UNKNOWN", Some(splitRead(5)), splitRead(4).toInt)
      readID += 1
    }
  }}

  def alignedReads(): Array[SequencingRead] = returnReads
}
