package aligner

import java.io.{PrintWriter, File}

import main.scala.utils.Utils
import utils.CutSites

import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.main.{ReadDirection, ReverseReadOrientation, ForwardReadOrientation, SequencingRead}
import scala.sys.process._

/**
 * Created by aaronmck on 11/17/15.
 */
case class MAFFT(forwardRead: SequencingRead, reverseRead: SequencingRead, referenceFile: String, cuts: CutSites) {

  // ----------------------------------------------------------------------
  // setup a BWA run and farm it out the machine
  // ----------------------------------------------------------------------
  val tmp = java.io.File.createTempFile("UMIMergerBWA", ".txt")
  val tmpWriter = new PrintWriter(tmp)
  var readDirections = Array[ReadDirection]()

  // write the reads / reference to the input file
  tmpWriter.write(">reference\n" + referenceFile + "\n")
  tmpWriter.write(">forward\n" + forwardRead.bases.filter { bs => bs != '-' }.mkString("") + "\n")
  tmpWriter.write(">reverse\n" + forwardRead.bases.filter { bs => bs != '-' }.mkString("") + "\n")
  tmpWriter.close()


  /**
   * Run MAFFT and capture the output
   */
  val out = new StringBuilder
  val err = new StringBuilder
  val logger = ProcessLogger(
    (o: String) => out.append(o + "\n"),
    (e: String) => err.append(e + "\n"))

  println("mafft --maxiterate 1000 --localpair " + tmp)
  ("mafft --maxiterate 1000 --localpair " + tmp) ! logger

  /**
   * now readback each read figure out the alignment
   */
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
      currentRead += line
    }
  }
  if (currentRead != "")
    readStrings :+= currentRead

  if (readNames.length != 3 || readStrings.length != 3)
    throw new IllegalStateException("Unable to parse out 3 reads")

  if (readStrings(1).length != readStrings(2).length || readStrings(1).length != readStrings(3).length)
    throw new IllegalStateException("The return strings are not the right size!")

  // find the first non-dash in the reverse, and the last non-dash in the forward
  val refFirstNonDash     = readStrings(0).zipWithIndex.filter(st => st._1 != '-').take(1)(0)._2
  val refLastNonDash      = readStrings(0).length - readStrings(1).reverse.zipWithIndex.filter(st => st._1 != '-').take(1)(0)._2

  val reverseFirstNonDash = readStrings(2).zipWithIndex.filter(st => st._1 != '-').take(1)(0)._2
  val forwardLastNonDash  = readStrings(1).length - readStrings(1).reverse.zipWithIndex.filter(st => st._1 != '-').take(1)(0)._2

  //println("refFirstNonDash " + refFirstNonDash + " reverseFirstNonDash " + reverseFirstNonDash + " forwardFirstNonDash " + forwardFirstNonDash)

  val callFromFwd = (refFirstNonDash, reverseFirstNonDash)
  val callFromRev = (refFirstNonDash, reverseFirstNonDash)

  var basePosition = 0
  var eventType = ""
  var eventLenth = 0
  var accumulatedEvents = new Array[String](10)

  // for each position in the reference
  for (pos <- refFirstNonDash until refLastNonDash) {

    // do a match on the three bases (ref, fwd, rev)
    (readStrings(0)(pos),readStrings(1)(pos),readStrings(2)(pos)) match {

      case('-','-','-') => {/* do nothing */}

      // we're in between the reads, flush the events and record any overlaps
      case(a ,b, c) if pos < reverseFirstNonDash && pos > forwardLastNonDash =>
      {
        // flush everything, we're in no-mans land
        if (eventLenth != 0)
          accumulatedEvents = MAFFT.updateHits(MAFFT.checkOverlap(basePosition - eventLenth,basePosition,eventType,cuts ), accumulatedEvents)
        eventType = ""
        eventLenth = 0
        basePosition += 1
      }

      // we have at least one non-dash in the reads
      case('-',b  ,c  ) if b != '-' || c != '-' =>
      {
        // flush everything, we're in no-mans land
        if (eventLenth != 0 && eventType != "I")
          accumulatedEvents = MAFFT.updateHits(MAFFT.checkOverlap(basePosition - eventLenth,basePosition,eventType,cuts ), accumulatedEvents)
        eventType = "I"
        eventLenth = 1
      }

      //
      case(a , b, c) if b == '-' && c == '-' =>
      {
        if (eventLenth != 0 && eventType != "D")
          accumulatedEvents = MAFFT.updateHits(MAFFT.checkOverlap(basePosition - eventLenth,basePosition,eventType,cuts ), accumulatedEvents)
        eventType = "D"
        eventLenth = 1
        basePosition += 1
      }

      case(a, b, c) => {
        // flush everything, we're in no-mans land
        if (eventLenth != 0)
          accumulatedEvents = MAFFT.updateHits(MAFFT.checkOverlap(basePosition - eventLenth,basePosition,eventType,cuts ), accumulatedEvents)
        eventType = ""
        eventLenth = 1
        basePosition += 1
      }
    }
  }
  if (eventLenth != 0)
    accumulatedEvents = MAFFT.updateHits(MAFFT.checkOverlap(basePosition - eventLenth,basePosition,eventType,cuts ), accumulatedEvents)

  print(accumulatedEvents.mkString("*"))
}

object MAFFT {
  def checkOverlap(start: Int, stop: Int, eventType: String, cutSites: CutSites): Option[Tuple2[Int,String]] = {
    cutSites.windows.zipWithIndex.foreach{case((cutfront,cutmiddle,cutrear),index) =>
      if ((start > cutfront && start < cutrear) || (stop > cutfront && stop < cutrear))
        return Some(index,(stop-start) + eventType + "-" + start)
    }
    return None
  }

  def updateHits(maybeHit: Option[Tuple2[Int,String]], eventSlots: Array[String]) : Array[String] = {
    if (!maybeHit.isDefined )
      return eventSlots
    val pos = maybeHit.get._1
    val event = maybeHit.get._2

    (eventSlots.slice(0,pos) ++ Array[String](event) ++ eventSlots.slice(pos+1, eventSlots.length))
  }
}