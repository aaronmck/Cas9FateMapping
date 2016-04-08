package main.scala

import _root_.aligner.{AlignmentManager, MAFFT}
import _root_.utils.CutSites
import aligner._
import main.scala.stats.{StatsContainer, StatsOutput}

import java.io._
import scala.main.SequencingRead
import main.scala.utils._

/**
  * handles the process of merging down read with matching UMIs
  */
object UMIMerger {

  def mergeTogether(umi: String,
                    readsF: Array[SequencingRead],
                    readsR: Array[SequencingRead],
                    readsFCount: Int,
                    readsRCount: Int,
                    outputFastq1: PrintWriter,
                    outputFastq2: PrintWriter,
                    primers: List[String],
                    sample: String,
                    minSurvivingReads: Int,
                    index: Int): Int = {


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

        outputFastq1.write(fwdConsensus.toFastqString(umi + "FWD" + "_" + fwdCleanedUp.size + "_" + readsFCount, false, index, 0) + "\n")
        outputFastq2.write(revConsensus.toFastqString(umi + "REV" + "_" + revCleanedUp.size + "_" + readsRCount, false, index, 0) + "\n")
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
