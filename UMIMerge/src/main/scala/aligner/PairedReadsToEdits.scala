package aligner

import java.io.File

import utils.CutSites

import scala.main.SequencingRead
import scala.util.matching.Regex

/**
 * Given two aligned reads, find their cigar strings compared to the reference
 */
case class PairedReadsToEdits(read1: SequencingRead, read2: SequencingRead, ref: String, cuts: CutSites) {

  // get a list of target sites
  val events = new Array[String](cuts.cutSites.size)

  // find the forward and reverse targets
  val fwdTargets = cigarString(read1,cuts)
  val revTargets = cigarString(read2,cuts)
  println("FWD:" + fwdTargets.mkString(":"))
  println("REV:" + revTargets.mkString(":"))


  // check that the forward and reverse events are compatible
  val mergedEvents = fwdTargets.zip(revTargets).zipWithIndex.map{case((fwdEvt, revEvt), index) => {
    if ((fwdEvt != "" || revEvt != "") && revEvt != fwdEvt)
      println("Warning We have a conflict at site " + index + " with events -" + fwdEvt + "-,-" + revEvt + "-")
    if (fwdEvt == "")
      revEvt
    else if (revEvt == "")
      fwdEvt
    else
      revEvt + "," + fwdEvt
  }}


  /**
   * given a read and a set of targets, make an array of events over the targets
   * @param read the read to look at -- it must have position and cigar information
   * @param targetSites the cut site object
   * @return an array of strings over the cutsites
   */
  def cigarString(read: SequencingRead, targetSites: CutSites): Array[String] = {

    val currentEvents = (0 until targetSites.size).map{st => ""}.toArray
    var currentPos = read.position
    // march along the cigar elements
    (new Regex("(\\d+)([MIDSNHP=X])") findAllIn read.cigar.get).matchData.foreach{matching => {
      matching.subgroups(1) match {

        // matching -- just update our read position
        case "M" => { currentPos += matching.subgroups(0).toInt }

        // insertion -- check that we are close enough to the cutsite
        case "I" => {
          val len = matching.subgroups(0).toInt
          targetSites.windows.zipWithIndex.foreach{case(pos,index) =>
            if ((currentPos < pos._2 && len + currentPos > pos._1) || (currentPos >= pos._2 && currentPos < pos._3)) {
              currentEvents(index) += len + "I-" + currentPos+ "_"
            }
          }
          //currentPos += len
        }

        case "D" => {
          val len = matching.subgroups(0).toInt
          targetSites.windows.zipWithIndex.foreach{case(pos,index) =>
            if ((currentPos < pos._2 && len + currentPos > pos._1) || (currentPos >= pos._2 && currentPos < pos._3)) {
              currentEvents(index) += len + "D-" + currentPos + "_"
            }
          }
          currentPos += len
        }
        case "S" => {
          // do nothing
        }
        case "H" => {
          // do nothing
        }}
    }}
    currentEvents
  }
}
