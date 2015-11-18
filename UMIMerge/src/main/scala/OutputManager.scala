package main.scala

import _root_.utils.{CutSites, CutSite, Clustering}
import aligner.{SmithWaterman, PairedReadsToEdits, MAFFT, Clustalo}
import aligner.aligner.Waterman

import scala.io._
import java.io._
import scala.collection.mutable._
import scala.main.SequencingRead
import scala.sys.process._
import java.util.zip._
import main.scala.utils._
import main.scala._

/**
 * Created by aaronmck on 10/22/15.
 */
object OutputManager {
  /*
    def outputReads(umi: String,
                    readsF: Array[SequencingRead],
                    readsR: Array[SequencingRead],
                    outputFASTAF: PrintWriter,
                    outputFASTAR: PrintWriter,
                    outputStats: PrintWriter,
                    ref: String,
                    primers: List[String],
                    sample: String,
                    minSurvivingReads: Int
    ): Tuple2[Int, Int] = {
      //: MergeResult = {


      // filter down the reads to those that have the primers on each end
      var newReadsF = Array[SequencingRead]()
      var newReadsR = Array[SequencingRead]()

      var matchingFor = 0
      var matchingRev = 0
      var total = 0
      var kept = 0

      readsF.zip(readsR).foreach { case(fwd,rev) => {
        val containsForward = fwd.startsWithPrimer(primers(0))
        val containsReverse = rev.startsWithPrimer(Utils.reverseComplement(primers(1)))

        if (containsForward) matchingFor += 1
        if (containsReverse) matchingRev += 1
        total += 1

        //println(reverseComplement(primers(1)))
        if (containsForward && containsReverse) {
          newReadsF :+= fwd
          newReadsR :+= rev
          kept += 1
        }
      }}

      val warn = true
      if (newReadsF.size < minSurvivingReads && warn) {
        println("WARN: dropping umi " + umi + " since we have less than two surviving reads. FM = " + matchingFor + " RV = " + matchingRev + " total = " + total)
        return (kept, total)
      } else {
        println("INFO: keeping umi " + umi + " since we have more than two surviving reads. FM = " + matchingFor + " RV = " + matchingRev + " total = " + total)

      }

      // read in the output from CLUSTAL and see how many reads we have
      val mergerF = Merger(newReadsF)
      val mergerR = Merger(newReadsR)

      // merge the reads down to a single consensus read
      val consensusReadF = Consensus.twoPassConsensus( mergerF.postMergedReads, mergerF.readLength, 0.10)
      val consensusReadR = Consensus.twoPassConsensus( mergerR.postMergedReads, mergerR.readLength, 0.10)

      val collapsedConcF = Consensus.collapseConsensus((consensusReadF._1, consensusReadF._2))
      val collapsedConcR = Consensus.collapseConsensus((consensusReadR._1, consensusReadR._2))

      println(collapsedConcF._1)
      println(collapsedConcR._1)
      // require that the consensus read has kept at least 80% of the original reads
      val highConcensusF = (consensusReadF._4 > 0.4)
      val highConcensusR = (consensusReadR._4 > 0.4)
      val enoughReadsKept = (newReadsF.size >= minSurvivingReads)

      val forwardPrimer = collapsedConcF._1.startsWithPrimer(primers(0))
      val reversePrimer = collapsedConcR._1.startsWithPrimer(Utils.reverseComplement(primers(1)))
      val usingRead = highConcensusF && highConcensusR && forwardPrimer && reversePrimer && enoughReadsKept

      var failureReason = ""
      if (!highConcensusF) failureReason += "lowConsensusF;"
      if (!highConcensusR) failureReason += "lowConsensusR;"
      if (!forwardPrimer) failureReason += "forwardPrimerMissing;"
      if (!reversePrimer) failureReason += "reversePrimerMissing;"
      if (!enoughReadsKept) failureReason += "notEnoughReadsRemaining;"
      if (failureReason == "") failureReason = "PASS"

      outputStats.write(umi + "\t" + usingRead + "\t" + failureReason + "\t") // write the UMI out and if we're using the read (and if it was rev or fwd reads that ruined it)
      outputStats.write(readsF.size + "\t" + newReadsF.size + "\t" + consensusReadF._3 + "\t" + +readsR.size + "\t" + newReadsR.size + "\t" + consensusReadR._3 + "\t") // output the before and after filtering counts of reads (fwd and rev)

      outputStats.write(collapsedConcF._2 + "\t" + collapsedConcF._1.bases + "\t")
      outputStats.write(collapsedConcR._2 + "\t" + collapsedConcR._1.bases + "\n")

      if (usingRead) {
        outputFASTAF.write(collapsedConcF._1.toFastqString(umi) + "\n")
        outputFASTAR.write(collapsedConcR._1.toFastqString(umi) + "\n")
      }
      return (kept, total)
    }*/

  def mergeTogether(umi: String,
                    readsF: Array[SequencingRead],
                    readsR: Array[SequencingRead],
                    cutSites: Array[CutSite],
                    outputStats: PrintWriter,
                    ref: String,
                    refFile: File,
                    primers: List[String],
                    sample: String,
                    minSurvivingReads: Int,
                    cutSiteInfo: CutSites
                     ): Tuple2[Int, Int] = {


    // get containers for the forward and reverse reads
    var newReadsF = Array[SequencingRead]()
    var newReadsR = Array[SequencingRead]()

    var matchingFor = 0
    var matchingRev = 0
    var total = 0
    var kept = 0

    readsF.zip(readsR).foreach { case (fwd, rev) => {
      val containsForward = fwd.startsWithPrimer(primers(0))
      val containsReverse = rev.startsWithPrimer(Utils.reverseComplement(primers(1)))

      if (containsForward) matchingFor += 1
      if (containsReverse) matchingRev += 1
      total += 1

      //println(reverseComplement(primers(1)))
      if (containsForward && containsReverse) {
        newReadsF :+= fwd
        newReadsR :+= rev
        kept += 1
      }
    }
    }

    val warn = true
    if (newReadsF.size < minSurvivingReads && warn) {
      println("WARN: dropping umi " + umi + " since we have less than two surviving reads. FM = " + matchingFor + " RV = " + matchingRev + " total = " + total)
      return (kept, total)
    } else {
      println("INFO: keeping umi " + umi + " since we have more than two surviving reads. FM = " + matchingFor + " RV = " + matchingRev + " total = " + total)

    }

    // some constants we should probably bubble-up
    val minReadLength = 30
    val minMeanQualScore = 30.0

    try {

      // use MSA to align all the reads
      val mergedF = Clustalo(newReadsF).postMergedReads
      val mergedR = Clustalo(newReadsR).postMergedReads

      // remove the reads that are a really poor match
      val fwdCleanedUp = Consensus.removeMismatchedReads(mergedF)
      val revCleanedUp = Consensus.removeMismatchedReads(mergedR)

      // make a consensus from the remaining 'good' reads
      val fwdConsensus = Array[SequencingRead](Consensus.consensus(fwdCleanedUp))
      val revConsensus = Array[SequencingRead](Consensus.consensus(revCleanedUp))

      // now align this to the reference using BWA
      val fwdAligned = SmithWaterman(fwdConsensus, refFile).alignedReads()
      val revAligned = SmithWaterman(revConsensus, refFile).alignedReads()

      println(fwdAligned.size)
      println(revAligned.size)
      // now get the cigar events for both strings
      val pairedEdits = PairedReadsToEdits(fwdAligned(0), revAligned(0),ref,cutSiteInfo)
      /*
      // remove the adapters, or anything that aligns past the ends of the reference on either side
      val adapterStripped = Consensus.referenceStrip(mergedF.postMergedReads ++ mergedR.postMergedReads)

      // remove bad reads -- reads that don't align to our target well
      //val failedReadsRemoved = Consensus.removeMismatchedReads(adapterStripped)

      // 2nd round: merge the passing reads with the reference and realign
      val screenReads2 = Consensus.dualConsensus(
        Consensus.prepareConsensus(adapterStripped, minReadLength, minMeanQualScore)
        ,minReadLength, minMeanQualScore)

      val secondRoundMerged = BWAMem(screenReads2, refFile).alignedReads()
      println(secondRoundMerged.size)

      // we have a foward and reverse primer in the reads
      val forwardPrimer = secondRoundMerged(0).bases contains primers(0)
      val reversePrimer = secondRoundMerged(1).bases contains primers(1) // (Utils.reverseComplement(primers(1)))
      val usingRead = forwardPrimer && reversePrimer

      // now get the cigar events for both strings
      val pairedEdits = PairedReadsToEdits(secondRoundMerged(0), secondRoundMerged(1),ref,cutSiteInfo)

      val readsKept = failedReadsRemoved.size.toDouble / (mergedF.postMergedReads ++ mergedR.postMergedReads).size.toDouble

      var failureReason = ""
      if (readsKept < 0.8) failureReason += "notEnoughReadsRemaining;"
      if (!forwardPrimer) failureReason += "forwardPrimerMissing;"
      if (!reversePrimer) failureReason += "reversePrimerMissing;"
      if (failureReason == "") failureReason = "PASS"

      // write out the stats file information
      outputStats.write(umi + "\t" + usingRead + "\t" + failureReason + "\t")
      outputStats.write(readsF.size + "\t" + newReadsF.size + "\t")
      outputStats.write(readsR.size + "\t" + newReadsR.size + "\t")
      outputStats.write(secondRoundMerged(0).cigar.get + "\t" + secondRoundMerged(1).cigar.get + "\t")
      outputStats.write(pairedEdits.mergedEvents.mkString("\t") + "\n")
      */
      return (kept, total)
    } catch {
      case e: Exception => {
        println("Unable to process UMI " + umi)
        e.printStackTrace()
        return (0, 0)
      }
    }
  }
}
