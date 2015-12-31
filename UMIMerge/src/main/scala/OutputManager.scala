package main.scala

import _root_.aligner.{AlignmentManager, MAFFT, AlignmentManager$}
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

  def mergeTogether(umi: String,
                    readsF: Array[SequencingRead],
                    readsR: Array[SequencingRead],
                    readsFCount: Int,
                    readsRCount: Int,
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


    // some constants we should probably bubble-up
    val minReadLength = 30
    val minMeanQualScore = 30.0

    try {

      // use MSA to align all the reads
      val preparedFWD = Consensus.prepareConsensus(readsF,minReadLength,minMeanQualScore)
      val preparedREV = Consensus.prepareConsensus(readsR,minReadLength,minMeanQualScore)

      //println("Pre forward: " + newReadsF.size + " post: " + preparedFWD.size + ", Pre reverse: " + newReadsR.size + " post: " + preparedREV.size)

      if (preparedFWD.size > 1 && preparedREV.size > 1) {

        val mergedF = MAFFT.alignTo(preparedFWD,None)
        val mergedR = MAFFT.alignTo(preparedREV,None)

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



          val merged = ReadMerger.mergeRead(fwdConsensus,revConsensus)
          if (merged.isDefined && merged.get._2 > 20 && merged.get._2.toDouble / merged.get._3.toDouble > 0.90) {
            // get the overlap of cutsite events
            try {
              val cutEvents = AlignmentManager.cutSiteEvent(umi, ref, merged.get._1, cutSites, 6)

              var failureReason: String = getFailureStatus(readsKept,cutEvents._1,cutEvents._2,forwardPrimer,reversePrimer, true)

              if (failureReason == "PASS") {
                outputFastq1.write(fwdConsensus.toFastqString(umi + "FWD", true) + "\n")
                outputFastq2.write(revConsensus.toFastqString(umi + "REV", false) + "\n")
              }

              val edit = cutEvents._4.map{al => al.toEditString}.mkString("-")
              //write out the stats file information
              outputStats.write(umi + "\t" + readsKept + "\t" + failureReason + "\t")
              outputStats.write(readsFCount + "\t" + readsF.size + "\t" + mergedF.size + "\t")
              outputStats.write(readsRCount + "\t" + readsR.size + "\t" + mergedR.size + "\t")
              outputStats.write(cutEvents._1 + "\t" + cutEvents._2 + "\t" + cutEvents._3.mkString("\t") + "\t")
              outputStats.write((if (edit == "") "UNALIGNED" else edit) + "\tmerged\t")
              outputStats.write(cutEvents._5 + "\t" + cutEvents._6 + "\tmerged\tmerged\n")
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
              val cutEvents = AlignmentManager.cutSiteEvents(umi, ref, fwdConsensus, revConsensus, cutSites, 6)

              var failureReason: String = getFailureStatus(readsKept,cutEvents._1,cutEvents._2,forwardPrimer,reversePrimer, false)

              if (failureReason == "PASS") {
                outputFastq1.write(fwdConsensus.toFastqString(umi + "FWD", true) + "\n")
                outputFastq2.write(revConsensus.toFastqString(umi + "REV", false) + "\n")
              }
              //write out the stats file information
              outputStats.write(umi + "\t" + readsKept + "\t" + failureReason + "\t")
              outputStats.write(readsFCount + "\t" + readsF.size + "\t" + mergedF.size + "\t")
              outputStats.write(readsRCount + "\t" + readsR.size + "\t" + mergedR.size + "\t")
              outputStats.write(cutEvents._1 + "\t" + cutEvents._2 + "\t" + cutEvents._3.mkString("\t") + "\t")
              outputStats.write(cutEvents._4.map{al => al.toEditString}.mkString("-") + "\t" + cutEvents._5.map{al => al.toEditString}.mkString("-") + "\t")
              outputStats.write(cutEvents._6 + "\t" + cutEvents._7 + "\t" + cutEvents._8 + "\t" + cutEvents._9 + "\n")
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
    } catch {
      case e: Exception => {
        println("Unable to process UMI " + umi)
        throw e
        //return (0, 0)
      }
    }
  }

  def getFailureStatus(readsKept: Double, fwdMatchBase: Double, revMatchBase: Double, fwdPrimer: Boolean, revPrimer: Boolean, merged: Boolean): String = {
    var failureReason = ""
    if (readsKept < 0.5) failureReason += "notEnoughReadsRemaining;"
    if (fwdMatchBase < 0.9) failureReason += "tooManyForwardMismatches;"
    if (!merged)
      if (revMatchBase < 0.9) failureReason += "tooManyForwardMismatches;"
    if (!fwdPrimer) failureReason += "forwardPrimerMissing;"
    if (!revPrimer) failureReason += "reversePrimerMissing;"

    if (failureReason == "")
      failureReason = "PASS"
    failureReason
  }


  /**
   * given a single read, align it to the reference
   */
  def processIndividualRead(fwdRead: SequencingRead,
                            revRead: SequencingRead,
                            cutSites: CutSites,
                            outputStats: PrintWriter,
                            ref: String,
                            refFile: File,
                            primers: List[String],
                            sample: String): Int = {

    val forwardPrimer = fwdRead.bases contains primers(0)
    val reversePrimer = Utils.reverseComplement(revRead.bases) contains primers(1)

    var failureReason = ""
    if (!forwardPrimer) failureReason += "forwardPrimerMissing;"
    if (!reversePrimer) failureReason += "reversePrimerMissing;"

    if (failureReason == "") failureReason = "PASS"

    val merged = ReadMerger.mergeRead(fwdRead,revRead)
    if (merged.isDefined && merged.get._2 > 20 && merged.get._2.toDouble / merged.get._3.toDouble > 0.90) {
      // get the overlap of cutsite events
      try {
        val cutEvents = AlignmentManager.cutSiteEvent(fwdRead.name, ref, merged.get._1, cutSites, 6)

        val edit = cutEvents._4.map{al => al.toEditString}.mkString("-")

        //write out the stats file information
        outputStats.write(fwdRead.name + "\t" + failureReason + "\t")
        outputStats.write(cutEvents._1 + "\t" + cutEvents._2 + "\t" + cutEvents._3.mkString("\t") + "\t")
        outputStats.write((if (edit == "") "UNALIGNED" else edit) + "\tmerged\t")
        outputStats.write(cutEvents._5 + "\t" + cutEvents._6 + "\tmerged\tmerged\n")
      } catch {
        case e: Exception => {
          println("ERROR: Unable to process UMI " + fwdRead.name + " tossing out!!!" + e.getMessage)
          return 0
        }
        // case e: Exception => throw e
      }
    } else {
      // get the overlap of cutsite events
      try {
        val cutEvents = AlignmentManager.cutSiteEvents(fwdRead.name, ref, fwdRead, revRead, cutSites, 6)

        //write out the stats file information
        outputStats.write(fwdRead.name + "\t" + failureReason + "\t")
        outputStats.write(cutEvents._1 + "\t" + cutEvents._2 + "\t" + cutEvents._3.mkString("\t") + "\t")
        outputStats.write(cutEvents._4.map{al => al.toEditString}.mkString("-") + "\t" + cutEvents._5.map{al => al.toEditString}.mkString("-") + "\t")
        outputStats.write(cutEvents._6 + "\t" + cutEvents._7 + "\t" + cutEvents._8 + "\t" + cutEvents._9 + "\n")
      } catch {
        case e: Exception => {
          println("ERROR: Unable to process UMI " + fwdRead.name + " tossing out!!!" + e.getMessage)
          return 0
        }
        // case e: Exception => throw e
      }
    }
    return 1
  }
}
