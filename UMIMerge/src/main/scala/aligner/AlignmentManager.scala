package aligner

import java.io.{PrintWriter, File}

import _root_.aligner.aligner.Waterman
import main.scala.Consensus
import main.scala.utils.Utils
import utils.CutSites

import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.main._
import scala.sys.process._


/**
 * a simple case class to hold alignments -- results we get back from parsing reads aligned with MAFFTv7
 * @param refPos the reference position for the start of the event
 * @param refBase the reference bases over the event
 * @param readBase the read bases over the event
 * @param cigarCharacter the cigar character for the event -- I, D, and M are valid
 */
case class Alignment(val refPos: Int, refBase: String, readBase: String, cigarCharacter: String) {
  def combine(next: Alignment): Array[Alignment] =
    if (next.cigarCharacter != cigarCharacter)
      return Array[Alignment](this, next)
    else
      return Array[Alignment](Alignment(this.refPos, this.refBase + next.refBase, this.readBase + next.readBase, cigarCharacter))

  def prettyPrint: String = refPos + ":" + refBase + ":" + readBase + ":" + cigarCharacter

  def toEditString: String = readBase.length + "" + cigarCharacter + "+" + refPos + {
    if (cigarCharacter == "I") {
      "+" + readBase
    } else {
      ""
    }
  }
}

object AlignmentManager {

  /**
   * call edits over a matched reference and read string
   * @param reference the reference string
   * @param read the read string, THE SAME LENGTH as the reference,  i.e. out of an MSA program
   * @param minMatchOnEnd the minimum number of matches to end our calling, if we don't see an event like this we backtrack Is and Ds until we find an M of this size
   * @param debugInfo should we dump a ton of debug info
   * @return a list of alignments over the read/ref combo
   */
  def callEdits(reference: String, read: String, minMatchOnEnd: Int, debugInfo: Boolean = false): List[Alignment] = {
    var referencePos = 0
    var inRef = false

    var refToEvent = List[Alignment]()

    reference.zip(read).foreach { case (refBase: Char, readBase: Char) =>
      //if (debugInfo)
      //  println("BASES: " + refBase + "," + readBase + " ")
      (refBase, readBase) match {
        case ('-', readB) if !inRef => {
          /* we might be in the situation where we haven't started the real alignment, take the offset */
        }
        case ('-', readB) if inRef => {
          // insertion

          if (refToEvent.isEmpty) refToEvent :+= Alignment(referencePos, refBase.toString, readBase.toString, "I")
          else refToEvent = refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos, refBase.toString, readBase.toString, "I"))

          if (debugInfo)
            println("1: " + refToEvent.map { st => st.prettyPrint }.mkString("<>") + " " + refToEvent.size)
        }
        case (refB, '-') if !inRef => {
          // deletion before read starts -- we haven't aligned yet
          referencePos += 1
        }
        case (refB, '-') => {
          // deletion
          inRef = true
          if (refToEvent.isEmpty) refToEvent :+= Alignment(referencePos, refBase.toString, readBase.toString, "D")
          else refToEvent = refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos, refBase.toString, readBase.toString, "D"))
          referencePos += 1
          if (debugInfo)
            println("2: " + refToEvent.map { st => st.prettyPrint }.mkString("<>") + " " + refToEvent.size)
        }
        case (refB, readB) => {
          // match / mismatch
          inRef = true
          if (refToEvent.isEmpty) refToEvent :+= Alignment(referencePos, refBase.toString, readBase.toString, "M")
          else refToEvent = refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos, refBase.toString, readBase.toString, "M"))
          referencePos += 1
          if (debugInfo)
            println("3: " + refToEvent.map { st => st.prettyPrint }.mkString("<>") + " " + refToEvent.size)
        }


      }
    }
    // get a bit aggressive here -- start from both ends -- strip off insertions and deletions until we hit a match or mismatch of at least 10 bases
    return filterEnds(refToEvent, minMatchOnEnd)
  }

  /**
   * filter the alignments, if we have poor matches on the ends that allow us to accept and indel, roll it back until we've matched enough bases
   * @param eventList the list of alignments over the read
   * @param minMatch the minimum number of match bases to 'anchor' the ends, otherwise strip the trash off
   * @return a filtered alignment set
   */
  def filterEnds(eventList: List[Alignment], minMatch: Int): List[Alignment] = {

    // make this as clear as possible
    var firstIndex = -1
    for (i <- 0 until eventList.size)
      if (firstIndex < 0 && eventList(i).cigarCharacter == "M" && eventList(i).readBase.length >= minMatch)
        firstIndex = i

    var lastIndex = -1
    for (i <- (eventList.size - 1).until(-1, -1))
      if (lastIndex < 0 && eventList(i).cigarCharacter == "M" && eventList(i).readBase.length >= minMatch)
        lastIndex = i

    //println(firstIndex+ " " + (lastIndex+1) + " " + "PRE: " + eventList.mkString("-") + " POST: " + eventList.slice(firstIndex,lastIndex+1).mkString("-"))
    return (eventList.slice(firstIndex, lastIndex + 1))
  }


  /**
   * given a read and reference, align and call events at the cut-sites
   * @param reference ref string
   * @param fwdRead read string
   * @param revRead reverse read string
   * @param cutSites the cutsutes to consider
   * @param minMatchOnEnd the minimum number of matches on the ends to keep from peeling crappy indels off
   * @param debug should we dump a lot of debug info
   * @return the rate of matching for cigar "M" bases for both reads and the array of events over cutsites
   */
  def cutSiteEvents(umi: String,
                    reference: String,
                    fwdRead: SequencingRead,
                    revRead: SequencingRead,
                    cutSites: CutSites,
                    minMatchOnEnd: Int,
                    debug: Boolean = false): Tuple9[Double, Double, Array[String],List[Alignment],List[Alignment], String, String, String, String] = {

    fwdRead.reverseCompAlign = false
    revRead.reverseCompAlign = true

    val alignmentsF = Waterman.alignTo(Array[SequencingRead](fwdRead), Some(reference), debug)
    val alignmentsR = Waterman.alignTo(Array[SequencingRead](revRead), Some(reference), debug)

    val events1 = AlignmentManager.callEdits(alignmentsF(0).bases, alignmentsF(1).bases, minMatchOnEnd, false)
    val events2 = AlignmentManager.callEdits(alignmentsR(0).bases, alignmentsR(1).bases, minMatchOnEnd, false)

    val combined = editsToCutSiteCalls(List[List[Alignment]](events1, events2), cutSites, debug)

    val matchRate1 = percentMatch(alignmentsF(0).bases, alignmentsF(1).bases)
    val matchRate2 = percentMatch(alignmentsR(0).bases, alignmentsR(1).bases)


    if (debug) {
      println("Events1 : \n" + events1.mkString("\n"))
      println("Events2 : \n" + events2.mkString(","))
    }

    return (matchRate1, matchRate2, combined._2,events1,events2,alignmentsF(0).bases,alignmentsF(1).bases,alignmentsR(0).bases,alignmentsR(1).bases)
  }

  /**
   * given a read and reference, align and call events at the cut-sites
   * @param reference ref string
   * @param mergedRead read string
   * @param cutSites the cutsutes to consider
   * @param minMatchOnEnd the minimum number of matches on the ends to keep from peeling crappy indels off
   * @param debug should we dump a lot of debug info
   * @return the rate of matching for cigar "M" bases for both reads and the array of events over cutsites
   */
  def cutSiteEvent(umi: String,
                    reference: String,
                    mergedRead: SequencingRead,
                    cutSites: CutSites,
                    minMatchOnEnd: Int,
                    debug: Boolean = false): Tuple6[Double, Double, Array[String],List[Alignment], String, String] = {

    mergedRead.reverseCompAlign = false

    val alignmentsMerged = Waterman.alignTo(Array[SequencingRead](mergedRead), Some(reference), debug)

    val events1 = AlignmentManager.callEdits(alignmentsMerged(0).bases, alignmentsMerged(1).bases, minMatchOnEnd, false)

    val combined = editsToCutSiteCalls(List[List[Alignment]](events1), cutSites, debug)

    val matchRate1 = percentMatch(alignmentsMerged(0).bases, alignmentsMerged(1).bases)


    if (debug) {
      println("Events1 : \n" + events1.mkString("\n"))
    }

    return (matchRate1, 0.0, combined._2,events1, alignmentsMerged(0).bases, alignmentsMerged(1).bases)
  }

  /**
   * do two intervals overlap
   * @param pos1Start
   * @param pos1End
   * @param pos2Start
   * @param pos2End
   * @return
   */
  def overlap(pos1Start: Int, pos1End: Int, pos2Start: Int, pos2End: Int): Boolean = (pos1Start, pos1End, pos2Start, pos2End) match {
    case (pos1S, pos1E, pos2S, pos2E) if pos1S <= pos2S && pos1E >= pos2E => true
    case (pos1S, pos1E, pos2S, pos2E) if pos2S <= pos1S && pos2E >= pos1E => true
    case (pos1S, pos1E, pos2S, pos2E) if pos1S <= pos2E && pos1E >= pos2S => true
    case (pos1S, pos1E, pos2S, pos2E) if pos2S <= pos1E && pos2E >= pos1S => true
    case _ => false
  }

  /**
   * combine the edits over two reads
   * @param edits the set of edits over the reads
   * @param cutSites the cutsites we consider
   * @param debug should we dump a lot of debugging info
   * @return a tuple of: an indicator if there was a collision between edits, and an array of events over the target cut sites
   */
  def editsToCutSiteCalls(edits: List[List[Alignment]], cutSites: CutSites, debug: Boolean = false): Tuple2[Boolean, Array[String]] = {
    var ret = Array[String]()
    var retCovered = Array[Boolean]()
    var collision = false

    cutSites.windows.foreach { case (start, cut, end) => {
      var candidates = Array[Alignment]()


      edits.foreach { editList =>
        editList.foreach { edit => {
          if ((edit.cigarCharacter == "D" && overlap(start, end, edit.refPos, edit.refPos + edit.refBase.length)) ||
            edit.cigarCharacter == "I" && overlap(start, end, edit.refPos, edit.refPos + 1))
            candidates :+= edit
        }
        }
      }

      if (debug)
        println("Site: " + start + "-" + end + ": " + candidates.mkString("\t") + "<<<")

      if (candidates.size == 0)
        ret :+= "NONE"
      else if (candidates.size == 1)
        ret :+= candidates(0).toEditString
      else if (candidates.size == 2) {
        // Do collision detection here -- do we have the same event?
        if (candidates(0).toEditString == candidates(1).toEditString)
          ret :+= candidates(0).toEditString
        else {
          ret :+= candidates(0).toEditString + "&" + candidates(1).toEditString
          collision = true
        }
      } else {
        ret :+= candidates.map { t => t.toEditString }.mkString("&")
        collision = true
      }
    }
    }

    return (collision, ret)
  }

  /**
   * for non gap bases, what is our matching proportion?
   * @param ref the reference string
   * @param read the read string of the same length as the reference string
   * @return a proportion of bases that match
   */
  def percentMatch(ref: String, read: String, minimumAlignedBases: Int = 50): Double = {
    var bases = 0
    var matches = 0
    ref.zip(read).foreach { case (refBase, readBase) => {
      if (refBase != '-' && readBase != '-') {
        if (refBase == readBase)
          matches += 1
        bases += 1
      }
    }
    }

    if (bases < minimumAlignedBases)
      return -1.0
    matches.toDouble / bases.toDouble
  }





}