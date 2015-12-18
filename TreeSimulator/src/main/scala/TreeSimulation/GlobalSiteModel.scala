package main.scala.TreeSimulation

import java.io.File

import cc.mallet.types.{Alphabet, Dirichlet}
import main.scala.Event

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.util.Random

/**
 * The goal here is to draw events for positions in our simulated trees, given the background probability of events at each site.
 */
case class GlobalSiteModel(leaves: Array[Event], eventSize: Int, alreadySpiltOn: Array[String]) {

  // for each of the target sites, develop a distribution of events within that site.
  // then, aggregate the chance of events happening at each site to a master distribution
  var eventDistributions = Array[NormalizedCounts[String]]()

  // how many events are != NONE
  var eventCounts = Array[Double]()

  // the name of each of our event sites
  var eventSites = Array[String]()

  var rand = new Random()

  // for each site in our target array, find the proportion of events at each position
  (0 until eventSize).foreach { evtPosition => {
    val eventToCount = new mutable.HashMap[String, Int]()

    // count the occurrences of each event type at this position
    // if we have multisite deletions only sample them at the rate of 1/(event target length)
    leaves.foreach { leaf => {
      if (!(alreadySpiltOn contains leaf.eventStrings(evtPosition)) && // remove sites we've already split on from the list
        leaf.eventStrings(evtPosition) != Event.NONE && // remove NONEs as split options
        rand.nextDouble < (1.0 / leaf.eventMap(leaf.eventStrings(evtPosition)).length)) // if we see an event that spans multiple sites, sample at 1 over the total
          eventToCount(leaf.eventStrings(evtPosition)) = eventToCount.getOrElse(leaf.eventStrings(evtPosition), 0) + 1
    }
    }

    val eventList = eventToCount.toList

    // now make an array of normalized values
    val total = eventList.map { case (key, value) => value }.sum
    val normalizedValues = eventList.map { case (key, count) => count.toDouble / total.toDouble }.toArray[Double]
    val alphabet = eventList.map { case (key, value) => key }.toArray

    eventDistributions :+= NormalizedCounts[String](normalizedValues, alphabet)

    eventCounts :+= eventList.map { case (key, count) => if (key == Event.NONE) 0.0 else count.toDouble }.sum
    eventSites :+= "Site_" + evtPosition

  }}

  // --------------------------------------------------
  // now create a distribution over all sites
  // --------------------------------------------------
  val siteLevelDistTmp = eventCounts.map { case (count) => count.toDouble / eventCounts.sum.toDouble }.toArray[Double]
  val siteLevelDist = NormalizedCounts[Int](siteLevelDistTmp,(0 until eventSize).toArray)
  val remainingSamples = eventCounts.sum > 0

}

/**
 * our static functions for the TreeModel class
 */
object GlobalSiteModel {

  /**
   * load events from the stats file and create a TreeModel
   *
   * @param statsFile the stats file
   */
  def createModel(statsFile: File, targetSiteCount: Int): Tuple2[GlobalSiteModel,mutable.HashMap[String, String]] = {
    val events = new ArrayBuffer[Event]()
    val uniqueEvents = new mutable.HashMap[String, String]()

    Source.fromFile(statsFile).getLines().drop(1).foreach { line => {
      val eventStrings = line.split("\t").slice(line.split("\t").size - (targetSiteCount * 2), line.split("\t").size - targetSiteCount)
      val eventAsString = eventStrings.mkString("_")
      if (!(uniqueEvents contains eventAsString)) {

        val isAllNone = eventStrings.map { evt => if (evt == "NONE") 0 else 1 }.sum

        if (isAllNone > 0) {
          val sp = line.split("\t")
          events += Event(sp(1), sp(0), sp(4).toInt, eventStrings, sp(0))
        }
      }
      uniqueEvents(eventAsString) = "Sample_" + line.split("\t")(0)
    }
    }

    println("Loaded " + events.size + " events from stat file " + statsFile)
    (GlobalSiteModel(events.toArray, targetSiteCount, Array[String]()),uniqueEvents)
  }

  /**
   * load events from the stats file and create a TreeModel
   *
   * @param statsFile the stats file
   */
  def createModelReturnEvents(statsFile: File, targetSiteCount: Int): Tuple2[GlobalSiteModel,Array[Event]] = {
    val events = new ArrayBuffer[Event]()
    val uniqueEvents = new mutable.HashMap[String, String]()

    Source.fromFile(statsFile).getLines().drop(1).foreach { line => {
      val eventStrings = line.split("\t").slice(line.split("\t").size - (targetSiteCount * 2), line.split("\t").size - targetSiteCount)
      val eventAsString = eventStrings.mkString("_")
      if (!(uniqueEvents contains eventAsString)) {

        val isAllNone = eventStrings.map { evt => if (evt == "NONE") 0 else 1 }.sum

        if (isAllNone > 0) {
          val sp = line.split("\t")
          events += Event(sp(1), sp(0), sp(4).toInt, eventStrings, sp(0))
        }
      }
      uniqueEvents(eventAsString) = "Sample_" + line.split("\t")(0)
    }
    }

    println("Loaded " + events.size + " events from stat file " + statsFile)
    (GlobalSiteModel(events.toArray, targetSiteCount, Array[String]()),events.toArray)
  }
}