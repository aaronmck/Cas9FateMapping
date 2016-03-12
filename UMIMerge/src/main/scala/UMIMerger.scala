package main.scala

import _root_.aligner.{AlignmentManager, MAFFT}
import _root_.utils.CutSites
import aligner._
import main.scala.stats.{StatsContainer, StatsOutput}

import java.io._
import scala.main.SequencingRead
import main.scala.utils._

/**
 * handles the somewhat complex process of merging down read with matching UMIs
 */
object UMIMerger {

  def mergeTogether(umi: String,
                    readsF: Array[SequencingRead],
                    readsR: Array[SequencingRead],
                    readsFCount: Int,
                    readsRCount: Int,
                    cutSites: CutSites,
                    outputFastq1: PrintWriter,
                    outputFastq2: PrintWriter,
                    outputStats: StatsOutput,
                    ref: String,
                    refFile: File,
                    primers: List[String],
                    sample: String,
                    minSurvivingReads: Int,
                    cutSiteInfo: CutSites
                     ): Int = {


    // some constants we should probably bubble-up
    val minReadLength = 30
    val minMeanQualScore = 30.0
    val debug = false

    // use MSA to align all the reads
    val preparedFWD = Consensus.prepareConsensus(readsF, minReadLength, minMeanQualScore)
    val preparedREV = Consensus.prepareConsensus(readsR, minReadLength, minMeanQualScore)

    if (preparedFWD.size > 1 && preparedREV.size > 1) {

      val mergedF = MAFFT.alignTo(preparedFWD, None)
      val mergedR = MAFFT.alignTo(preparedREV, None)

      // remove the reads that are a really poor match
      val fwdCleanedUp = Consensus.removeMismatchedReads(mergedF)
      val revCleanedUp = Consensus.removeMismatchedReads(mergedR)

      if (fwdCleanedUp.size > 1 && revCleanedUp.size > 1) {

        // make a consensus from the remaining 'good' reads
        val fwdConsensus = SequencingRead.stripDownToJustBases(Consensus.consensus(fwdCleanedUp, "ConsensusFWD"))
        val revConsensus = SequencingRead.stripDownToJustBases(Consensus.consensus(revCleanedUp, "ConsensusREV"))
        revConsensus.reverseCompAlign = true

        val forwardPrimer = fwdConsensus.bases contains primers(0)
        val reversePrimer = Utils.reverseComplement(revConsensus.bases) contains primers(1)
        val readsKept = (mergedF.size + mergedR.size).toDouble / (readsF.size + readsR.size).toDouble

        if (debug) {
          println("UMI::::: " + umi)
        }

        val merged = ReadMerger.mergeRead(fwdConsensus, revConsensus)
        if (debug && merged.isDefined) {
          println("Merged candidate ")
          println(merged.get.matches)
          println(merged.get.overlap)
          println(merged.get.orientationOK)

        }

        if (merged.isDefined && merged.get.matches > 20 && merged.get.matches.toDouble / merged.get.overlap.toDouble > 0.95 && merged.get.orientationOK) {
          // get the overlap of cutsite events
          try {
            //if (debug) {
              println("MERGED")
            //}
            val mergedRefPair = RefReadPair(SequencingRead.readFromNameAndSeq("ref",ref),merged.get.read)
            val cutEvents = AlignmentManager.cutSiteEvent(mergedRefPair, cutSites,debug=debug)

            var failureReason: String = getFailureStatus(readsKept, cutEvents.matchingRate, None, forwardPrimer, reversePrimer)

            if (failureReason == "PASS" && !cutEvents.collision) {
              outputFastq1.write(fwdConsensus.toFastqString(umi + "FWD", true) + "\n")
              outputFastq2.write(revConsensus.toFastqString(umi + "REV", false) + "\n")
            }

            outputStats.outputStatEntry(StatsContainer(umi, failureReason == "PASS" && !cutEvents.collision, forwardPrimer, reversePrimer,
              true, true, readsF.size, readsR.size, mergedF.size, mergedR.size, cutEvents.matchingRate, -1, cutEvents.matchingBaseCount, -1,
              cutEvents.alignments, cutEvents.basesOverTargets, None, None, Some(mergedRefPair.read.bases), None, None, Some(mergedRefPair.reference.bases)))

          } catch {
            case e: Exception => {
              println("ERROR: Unable to process UMI " + umi + " tossing out!!!" + e.getMessage)
              return 0
            }
            // case e: Exception => throw e
          }
        } else {
          // get the overlap of cutsite events
          try {
            if (debug) {
              println("SINGLE READS")
            }
            val fwdPair = RefReadPair(SequencingRead.readFromNameAndSeq("ref",ref),fwdConsensus)
            val revPair = RefReadPair(SequencingRead.readFromNameAndSeq("ref",ref),revConsensus)
            val cutEvents = AlignmentManager.cutSiteEventsPair(fwdPair, revPair, cutSites,debug=debug)

            var failureReason: String = getFailureStatus(readsKept, cutEvents.matchingRate1, Some(cutEvents.matchingRate2), forwardPrimer, reversePrimer)

            if (failureReason == "PASS" && !cutEvents.collision) {
              outputFastq1.write(fwdConsensus.toFastqString(umi + "FWD", true) + "\n")
              outputFastq2.write(revConsensus.toFastqString(umi + "REV", false) + "\n")
            }

            outputStats.outputStatEntry(StatsContainer(umi, failureReason == "PASS" && !cutEvents.collision, forwardPrimer, reversePrimer,
              true, false, readsF.size, readsR.size, mergedF.size, mergedR.size, cutEvents.matchingRate1, cutEvents.matchingRate2, cutEvents.matchingBaseCount1, cutEvents.matchingBaseCount2,
              cutEvents.alignments, cutEvents.basesOverTargets, Some(cutEvents.read1), Some(cutEvents.read2), None, Some(cutEvents.read1Ref), Some(cutEvents.read1Ref), None))

          } catch {
            case e: Exception => {
              println("ERROR: Unable to process UMI " + umi + " tossing out!!!" + e.getMessage)
              return 0
            }
            // case e: Exception => throw e
          }
        }


      }
    }
    return 1

  }

  def getFailureStatus(readsKept: Double, fwdMatchBase: Double, revMatchBase: Option[Double], fwdPrimer: Boolean, revPrimer: Boolean): String = {
    var failureReason = ""
    if (readsKept < 0.5) failureReason += "notEnoughReadsRemaining;"
    if (fwdMatchBase < 0.85) failureReason += "tooManyForwardMismatches;"
    if (revMatchBase.getOrElse(1.0) < 0.85) failureReason += "tooManyForwardMismatches;"
    if (!fwdPrimer) failureReason += "forwardPrimerMissing;"
    if (!revPrimer) failureReason += "reversePrimerMissing;"

    if (failureReason == "")
      failureReason = "PASS"
    failureReason
  }

}
