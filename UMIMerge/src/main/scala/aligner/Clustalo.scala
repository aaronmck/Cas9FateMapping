package aligner

import java.io._

import scala.collection.mutable._
import scala.io._
import scala.main.SequencingRead
import scala.sys.process._

/**
 * Created by aaronmck on 10/22/15.
 */
// ------------------------------------------------------------------------------------------
//
// process reads to and from the clustalo algorithm
//
// ------------------------------------------------------------------------------------------
case class Clustalo(reads: Array[SequencingRead]) extends Aligner {

  var currentReadName: Option[String] = None
  var currentRead = ""
  var readLength = 0


  // map the input read name to the read, so we can reset the UMI, direction, etc
  // ------------------------------------------------------------------------------
  val mappedReads = new HashMap[String,SequencingRead]()
  reads.foreach{rd => mappedReads(rd.name) = rd}


  // our temp file we've written to
  // ------------------------------------------------------------------------------
  val tmpFile = clustal_reads()
  val postMergedReads = fetchReads(tmpFile)


  // map our input reads from their name
  // ------------------------------------------------------------------------------
  def fetchReads(alignedFile: File): Array[SequencingRead] = {
    var returnReads = Array[SequencingRead]()
    Source.fromFile(alignedFile).getLines().foreach { line => {
      if (line startsWith ">") {
        if (currentReadName.isDefined) {
          if (!(mappedReads contains currentReadName.get))
            throw new IllegalStateException("Unable to find the read we send out: " + currentReadName.get + " in our input reads")
          val oldRead = mappedReads(currentReadName.get)
          returnReads :+= SequencingRead(currentReadName.get, currentRead, {"H" * currentRead.length}.mkString(""), oldRead.readOrientation, oldRead.umi)

          readLength = currentRead.length
          currentRead = ""
          currentReadName = Some(line.stripPrefix(">"))
        } else {
          currentReadName = Some(line.stripPrefix(">"))
        }
      } else {
        currentRead += line
      }
    }
    }

    // if we have a read, save it off
    // ------------------------------------------------------------------------------
    if (currentReadName.isDefined) {
      if (!(mappedReads contains currentReadName.get))
        throw new IllegalStateException("Unable to find the read we send out: " + currentReadName.get + " in our input reads")
      val oldRead = mappedReads(currentReadName.get)

      returnReads :+= SequencingRead(currentReadName.get, currentRead, {
        "H" * currentRead.length
      }.mkString(""), oldRead.readOrientation, oldRead.umi)
    }
    return(returnReads)
  }

  /**
   * do the actual clustal run with our reads
   * @return the file containing the alignment
   */
  def clustal_reads(): File = {
    // setup a CLUSTAL run and farm it out the machine
    val tmp = java.io.File.createTempFile("UMIMerger", ".txt")
    val tmpWriter = new PrintWriter(tmp)
    reads.foreach { case (read) => tmpWriter.write("> " + read.name + "\n" + read.bases.filter{bs => bs != '-'}.mkString("") + "\n") }
    tmpWriter.close()
    val tmpOutput = java.io.File.createTempFile("UMIMerger", ".post.txt")
    val clustoResult = ("clustalo --dealign --force --pileup -i " + tmp + " -o " + tmpOutput).!!
    //println("Clustal: Input " + tmp + " and output " + tmpOutput)
    return tmpOutput
  }


}
