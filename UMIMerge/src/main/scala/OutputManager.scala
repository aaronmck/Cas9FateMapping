package main.scala

import _root_.utils.{CutSites, Clustering}
import aligner._
import _root_.aligner.aligner.Waterman

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
                    cutSites: CutSites,
                    outputFastq1: PrintWriter,
                    outputFastq2: PrintWriter,
                    outputStats: PrintWriter,
                    ref: String,
                    refFile: File,
                    primers: List[String],
                    sample: String,
                    minSurvivingReads: Int,
                    cutSiteInfo: CutSites
                     ): Int = {


    // get containers for the forward and reverse reads
    var newReadsF = Array[SequencingRead]()
    var newReadsR = Array[SequencingRead]()

    var matchingFor = 0
    var matchingRev = 0

    readsF.zip(readsR).foreach { case (fwd, rev) => {
      val containsForward = fwd.startsWithPrimer(primers(0))
      val containsReverse = rev.startsWithPrimer(Utils.reverseComplement(primers(1)))

      if (containsForward) matchingFor += 1
      if (containsReverse) matchingRev += 1

      //println(reverseComplement(primers(1)))
      if (containsForward && containsReverse) {
        newReadsF :+= fwd
        newReadsR :+= rev
      }
    }
    }

    val warn = true
    if (matchingFor < minSurvivingReads || matchingRev < minSurvivingReads) {
      ///println("WARN: dropping umi " + umi + " since we have less than two surviving reads. FM = " + matchingFor + " RV = " + matchingRev)
      return 0
    } else {
      //println("INFO: keeping umi " + umi + " since we have more than two surviving reads. FM = " + matchingFor + " RV = " + matchingRev)
    }

    // some constants we should probably bubble-up
    val minReadLength = 30
    val minMeanQualScore = 30.0

    try {

      // use MSA to align all the reads
      val preparedFWD = Consensus.prepareConsensus(newReadsF,minReadLength,minMeanQualScore)
      val preparedREV = Consensus.prepareConsensus(newReadsR,minReadLength,minMeanQualScore)

      //println("Pre forward: " + newReadsF.size + " post: " + preparedFWD.size + ", Pre reverse: " + newReadsR.size + " post: " + preparedREV.size)

      if (preparedFWD.size > 1 && preparedREV.size > 1) {

        val mergedF = MAFFT2.alignTo(preparedFWD,None)
        val mergedR = MAFFT2.alignTo(preparedREV,None)

        // remove the reads that are a really poor match
        val fwdCleanedUp = Consensus.removeMismatchedReads(mergedF)
        val revCleanedUp = Consensus.removeMismatchedReads(mergedR)

        if (fwdCleanedUp.size > 1 && revCleanedUp.size > 1) {

          // make a consensus from the remaining 'good' reads
          val fwdConsensus = SequencingRead.stripDownToJustBases(Consensus.consensus(fwdCleanedUp))
          val revConsensus = SequencingRead.stripDownToJustBases(Consensus.consensus(revCleanedUp))

          //val mafftAlign = MAFFT2(fwdConsensus(0), revConsensus(0), ref, cutSiteInfo, primers)

          val forwardPrimer = fwdConsensus.bases contains primers(0)
          val reversePrimer = Utils.reverseComplement(revConsensus.bases) contains primers(1)
          val readsKept = (mergedF.size + mergedR.size).toDouble / (readsF.size + readsR.size).toDouble

          var failureReason = ""
          if (readsKept < 0.5) failureReason += "notEnoughReadsRemaining;"
          if (!forwardPrimer) failureReason += "forwardPrimerMissing;"
          if (!reversePrimer) failureReason += "reversePrimerMissing;"
          //if (!enoughReadsKept) failureReason += "notEnoughReadsRemaining;"

          if (failureReason == "") {
            failureReason = "PASS"
            outputFastq1.write(fwdConsensus.toFastqString(umi + "FWD", true) + "\n")
            outputFastq2.write(revConsensus.toFastqString(umi + "REV", false) + "\n")
          }

          // get the overlap of cutsite events
          val cutEvents = MAFFT2.cutSiteEvents(ref,fwdConsensus,revConsensus,cutSites,6)

          //write out the stats file information
          outputStats.write(umi + "\t" + readsKept + "\t" + failureReason + "\t")
          outputStats.write(readsF.size + "\t" + newReadsF.size + "\t" + mergedF.size + "\t")
          outputStats.write(readsR.size + "\t" + newReadsR.size + "\t" + mergedR.size + "\t")
          outputStats.write(fwdConsensus.bases + "\t" + revConsensus.bases + "\t" + cutEvents._1 + "\t" + cutEvents._2 + "\t" + cutEvents._3.mkString("\t") + "\n")
          //outputStats.write(pairedEdits.mergedEvents.mkString("\t") + "\n") */




        }
      }
      return 1
    } catch {
      case e: Exception => {
        println("Unable to process UMI " + umi)
        throw e
        //return (0, 0)
      }
    }
  }
}
