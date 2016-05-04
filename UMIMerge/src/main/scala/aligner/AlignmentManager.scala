package aligner

import java.io.{PrintWriter, File}

import _root_.aligner.aligner.Waterman
import main.scala.Consensus
import main.scala.utils.{RefReadPair, Utils}
import utils.CutSites

import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.main._
import scala.sys.process._


/**
 * a simple case class to hold alignments -- results we get back from parsing reads aligned with MAFFTv7 or other aligners
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

  def pctBasesMatching(): Double = refBase.toUpperCase.zip(readBase.toUpperCase).map{case(b1,b2) => if (b1 == b2) 1 else 0}.sum / refBase.length.toDouble
}

object AlignmentManager {

  /**
   * call edits over a matched reference and read string
   * @param reference the reference string
   * @param read the read string, THE SAME LENGTH as the reference,  i.e. out of an MSA program
   * @param minMatchOnEnd the minimum number of matches to end our calling, if we don't see an event like this we backtrack Is and Ds until we find an M of this size
   * @param debugInfo should we dump a ton of debug info
   * @return a list of alignments over the read/ref combo, and a list of sequences over the target
   */
  def callEdits(reference: String, read: String, minMatchOnEnd: Int, cutSites: CutSites, debugInfo: Boolean = false): Tuple2[List[Alignment], List[String]] = {

    if (reference.length != read.length)
      throw new IllegalStateException("The read and reference lengths are unequal!")

    var referencePos = 0
    var inRef = false

    var refToEvent = List[Alignment]()
    var targetIndexToSequence = List[StringBuilder]()
    cutSites.fullSites.foreach { fullCutSite => targetIndexToSequence :+= new StringBuilder() }

    if (debugInfo) {
      println(reference)
      println(read)
    }
    reference.zip(read).foreach { case (refBase: Char, readBase: Char) => {

      // first add to our target sequence object if we overlap target
      cutSites.fullSites.zipWithIndex.foreach { case (fullCut, index) =>
        if (fullCut._2 <= referencePos && fullCut._3 >= referencePos) {
          if (debugInfo) {
            println("adding reference " + referencePos + " to window " + fullCut._2 + "," + fullCut._3)
          }
          targetIndexToSequence(index) += readBase
        }
      }

      // now we match the ref and read bases to see if we have an indel
      (refBase, readBase) match {
        case ('-', readB) if !inRef => {
          /* the situation where we haven't started the real alignment, the read has aligned off the end of the reference */
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
    }}

    // return the filtered list of the edits, plus the list of sequences over each target region
    // for filtering get a bit aggressive here -- start from both ends -- strip off insertions and deletions until we hit a match or mismatch of at least 10 bases
    return (filterEnds(refToEvent, minMatchOnEnd, debugInfo),targetIndexToSequence.map{bld => bld.result()}.toList)
  }

  /**
   * filter the alignments, when we have poor matches on the ends that allow us to accept and indel, peel back alignments until we've matched enough bases
   * @param eventList the list of alignments over the read
   * @param minMatch the minimum number of match bases to 'anchor' the ends, otherwise strip the trash off
   * @return a filtered alignment set
   */
  def filterEnds(eventList: List[Alignment], minMatch: Int, debugInfo: Boolean = false, matchingBasePCT: Double = 0.60): List[Alignment] = {

    // find the first and last alignments that have at least X bases matched that aren't Ns
    var firstIndex = -1
    for (i <- 0 until eventList.size)
      if (firstIndex < 0 && eventList(i).cigarCharacter == "M" &&  // have to be a match
        eventList(i).pctBasesMatching > matchingBasePCT &&  // and have over 60% of the bases match
        (eventList(i).readBase.length - eventList(i).refBase.count(p => p == 'N')) >= minMatch) // and meet the minimum length, not counting Ns
        firstIndex = i

    var lastIndex = -1
    for (i <- (eventList.size - 1).until(-1, -1))
      if (lastIndex < 0 && eventList(i).cigarCharacter == "M" && // have to be a match
        eventList(i).pctBasesMatching > matchingBasePCT && // and have over 60% of the bases match
        (eventList(i).readBase.length  - eventList(i).refBase.count(p => p == 'N')) >= minMatch) // and meet the minimum length, not counting Ns
        lastIndex = i

    if (debugInfo) {
      println(firstIndex + " " + (lastIndex + 1) + " " + "PRE: " + eventList.mkString("-") + " POST: " + eventList.slice(firstIndex, lastIndex + 1).mkString("-"))
    }
    return (eventList.slice(firstIndex, lastIndex + 1))
  }


  /**
   * given a read and reference, align and call events at the cut-sites
   * @param fwdRead read string
   * @param revRead reverse read string
   * @param cutSites the cutsutes to consider
   * @param minMatchOnEnd the minimum number of matches on the ends to keep from peeling crappy indels off
   * @param debug should we dump a lot of debug info
   * @return the rate of matching for cigar "M" bases for both reads and the array of events over cutsites
   */
  def cutSiteEventsPair(fwdRead: RefReadPair,
                    revRead: RefReadPair,
                    cutSites: CutSites,
                    minMatchOnEnd: Int = 8,
                    debug: Boolean = false): PairedReadCutSiteEvent = {

    fwdRead.read.reverseCompAlign = false
    revRead.read.reverseCompAlign = true

    val alignmentsF = Array[SequencingRead](fwdRead.reference,fwdRead.read)
    val alignmentsR = Array[SequencingRead](revRead.reference,revRead.read)

    val events1 = AlignmentManager.callEdits(alignmentsF(0).bases, alignmentsF(1).bases, minMatchOnEnd, cutSites)
    val events2 = AlignmentManager.callEdits(alignmentsR(0).bases, alignmentsR(1).bases, minMatchOnEnd, cutSites)

    val combined = editsToCutSiteCalls(
      List[List[Alignment]](events1._1, events2._1),
      List[List[String]](events1._2, events2._2), cutSites, debug)

    val matchRate1 = percentMatch(alignmentsF(0).bases, alignmentsF(1).bases)
    val matchRate2 = percentMatch(alignmentsR(0).bases, alignmentsR(1).bases)

    return PairedReadCutSiteEvent(matchRate1._1,
      matchRate1._2,
      matchRate2._1,
      matchRate2._2,
      combined._2,
      combined._3,
      combined._1,
      alignmentsF(1).bases,
      alignmentsF(0).bases,
      alignmentsR(1).bases,
      alignmentsR(0).bases
    )
  }

  // an inline case class to make the return of a cutsite call more readable
  case class PairedReadCutSiteEvent(matchingRate1: Double,
                                    matchingBaseCount1: Int,
                                    matchingRate2: Double,
                                    matchingBaseCount2: Int,
                                    alignments: Array[String],
                                    basesOverTargets: Array[String],
                                    collision: Boolean,
                                    read1: String,
                                    read1Ref: String,
                                    read2: String,
                                    read2Ref: String)

  /**
   * given a read and reference, align and call events at the cut-sites
   * @param mergedRead read object
   * @param cutSites the cutsutes to consider
   * @param minMatchOnEnd the minimum number of matches on the ends to keep from peeling crappy indels off
   * @param debug should we dump a lot of debug info
   * @return a cutsite event object
   */
  def cutSiteEvent(mergedRead: RefReadPair,
                   cutSites: CutSites,
                   minMatchOnEnd: Int = 8,
                   debug: Boolean = false): SingleReadCutSiteEvent = {

    mergedRead.read.reverseCompAlign = false

    val alignmentsMerged = Array[SequencingRead](mergedRead.reference, mergedRead.read)

    val events1 = AlignmentManager.callEdits(alignmentsMerged(0).bases, alignmentsMerged(1).bases, minMatchOnEnd, cutSites)

    val combined = editsToCutSiteCalls(List[List[Alignment]](events1._1), List[List[String]](events1._2),cutSites, debug)

    val matchRate1 = percentMatch(alignmentsMerged(0).bases, alignmentsMerged(1).bases)

    return SingleReadCutSiteEvent(matchRate1._1, matchRate1._2, combined._2, combined._3, combined._1)
  }

  // an inline case class to make the return of a cutsite call more readable
  case class SingleReadCutSiteEvent(matchingRate: Double, matchingBaseCount: Int, alignments: Array[String],basesOverTargets: Array[String], collision: Boolean)


  def overlap(pos1Start: Int, pos1End: Int, pos2Start: Int, pos2End: Int): Boolean = (pos1Start, pos1End, pos2Start, pos2End) match {
    case (pos1S, pos1E, pos2S, pos2E) if pos1S <= pos2S && pos1E >= pos2E => true
    case (pos1S, pos1E, pos2S, pos2E) if pos2S <= pos1S && pos2E >= pos1E => true
    case (pos1S, pos1E, pos2S, pos2E) if pos1S <= pos2E && pos1E >= pos2S => true
    case (pos1S, pos1E, pos2S, pos2E) if pos2S <= pos1E && pos2E >= pos1S => true
    case _ => false
  }

  // does position span at least the whole length of pos 2? THIS IS DIRECTIONAL, watch out!
  def span(pos1Start: Int, pos1End: Int, pos2Start: Int, pos2End: Int): Boolean = {
    if (pos1Start < pos2Start && pos1End > pos2End)
      true
    else
      false
  }

  /**
   * combine the edits over two reads -- this is bit complicated, as we have to check for collisions between aligned reads
   * @param readAlignments the set of edits aggregated by aligned read
   * @param readSequences the read sequences -- we have to choose the right one in the two read case
   * @param cutSites the cutsites over the reference we consider
   * @param debug should we dump a lot of debugging info
   * @return a tuple3 of: an indicator if there was a collision between edits, an array of events over the target cut sites, and an array of reference sequences
   */
  def editsToCutSiteCalls(readAlignments: List[List[Alignment]],
                          readSequences: List[List[String]],
                          cutSites: CutSites,
                          debug: Boolean = false): Tuple3[Boolean, Array[String], Array[String]] = {

    if (readAlignments.size != readSequences.size) {
      throw new IllegalStateException("Unable to call edits from read alignments of length " + readAlignments.size + " and read sequences of length " + readSequences.size)
    }

    var ret = Array[String]()
    var retTargetSeq = Array[String]()
    var collision = false

    cutSites.windows.zipWithIndex.foreach { case ((start, cut, end), cutSiteIndex) => {
      var candidates = Array[String]()
      var cigarMatchOverlapsEdit = false
      var nonWildType = Array[String]()
      var referenceSeq = Array[String]()

      readAlignments.zipWithIndex.foreach { case(singleReadEdits,readIndex) => {
        var singleSampleEvent = Array[String]()
        singleReadEdits.foreach { edit => {

          referenceSeq :+= readSequences(readIndex)(cutSiteIndex)

          if ((edit.cigarCharacter == "D" && overlap(start, end, edit.refPos, edit.refPos + edit.refBase.length)) ||
            edit.cigarCharacter == "I" && overlap(start, end, edit.refPos, edit.refPos)) {
            // check that we haven't already added this exact edit to the list -- this will happen in paired reads where the edits agree

            if (!(nonWildType contains readSequences(readIndex)(cutSiteIndex)))
              nonWildType :+= readSequences(readIndex)(cutSiteIndex)

            singleSampleEvent :+= edit.toEditString

          } else if (edit.cigarCharacter == "M" && span(edit.refPos, edit.refPos + edit.refBase.length, start, end)) {
            cigarMatchOverlapsEdit = true
          }
        }}
        if (singleSampleEvent.size > 0 && !(candidates contains singleSampleEvent.mkString("&")))
          candidates :+= singleSampleEvent.mkString("&")
        }
      }

      if (debug) {
        println("Site: " + start + "-" + end + ": " + candidates.mkString("\t") + "<<<")
      }

      // create target strings
      nonWildType.size match {
        case 0 => retTargetSeq :+= {
          referenceSeq.size match {
            case 0 => "BLANK"
            case 1 => referenceSeq(0)
            case _ => {
              // this is an edge case we should do better on: both reads 'overlap', we choose the one with the reference sequence
              val nonWTbases = referenceSeq.map{sq => sq.map{base => if (base == '-') 0 else 1}.sum}.toArray
              referenceSeq(nonWTbases.zipWithIndex.maxBy(_._1)._2)
            }
          }
        }
        case 1 => retTargetSeq :+= nonWildType(0)
        case _ => retTargetSeq :+= nonWildType.mkString(",")
      }

      // look at the number of candidate events,
      // if there is WT sequence covering the locus,
      // figure out what we should do with the edit
      (candidates.size, cigarMatchOverlapsEdit) match {
        case (0, true) => {
          ret :+= "NONE"
        }
        case (0, false) => {
          ret :+= "UNKNOWN"
        }
        case (1, true) => {
          ret :+= "WT_" + candidates(0)
        }
        case (1, false) => {
          ret :+= candidates(0)
        }
        case (2, true) if (candidates(0)== candidates(1)) => {
          ret :+= "WT_" + candidates(0)
        }
        case (2, false) if (candidates(0)== candidates(1)) => {
          ret :+= candidates(0)
        }
        case (2, true) => {
          collision = true
          ret :+= "WT_" + candidates(0)+ "&" + candidates(1)
        }
        case (2, false) => {
          collision = true
          ret :+= candidates(0)+ "&" + candidates(1)
        }
        case _ => {
          ret :+= candidates.mkString("&")
          collision = true
        }
      }
    }
    }

    return (collision, ret, retTargetSeq)
  }

  /**
   * for non gap bases, what is our matching proportion?
   * @param ref the reference string
   * @param read the read string of the same length as the reference string
   * @return a proportion of bases that match, and the count of non-gap bases
   */
  def percentMatch(ref: String, read: String, minimumAlignedBases: Int = 50): Tuple2[Double, Int] = {
    var bases = 0
    var matches = 0
    var nonGap = 0
    ref.zip(read).foreach { case (refBase, readBase) => {
      if (refBase != '-' && readBase != '-') {
        if (refBase == readBase) {
          matches += 1
          nonGap += 1
        }
        bases += 1
      }
    }
    }

    if (bases < minimumAlignedBases)
      (-1.0, 0)
    else
      (matches.toDouble / bases.toDouble, nonGap)
  }
}