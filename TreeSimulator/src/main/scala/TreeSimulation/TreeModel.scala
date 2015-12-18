package main.scala.TreeSimulation

import java.io.File

import cc.mallet.types.{Alphabet, Dirichlet}
import main.scala.Event

import scala.collection.mutable
import scala.io.Source
import scala.util.Random

/**
 * The goal here is to draw random trees, given the background probability of events at each site.
 */
case class TreeModel(leaves: Array[Event], eventSize: Int) {

  // for each of the Event spaces, develop a distribution of events within that site.
  // then, aggregate the chance of events happening at each site to a master distribution
  var eventDistributions = Array[NumericalDistribution[String]]()

  var eventCounts = Array[Double]()
  var eventSites = Array[String]()

  var rand = new Random()


  // for each position in our event array, find the proportion
  // of events at each position
  (0 until eventSize).foreach{evtPosition => {
    val eventToCount = new mutable.HashMap[String,Int]()

    // count the occurrences of each event type at this position
    // if we have multisite deletions only sample them at the rate of 1/(site length)
    leaves.foreach{leaf => {
      // if we see an event that spans multiple sites, sample at 1 over the total
      if (leaf.eventStrings(evtPosition) == Event.NONE || rand.nextDouble < (1.0 / leaf.eventMap(leaf.eventStrings(evtPosition)).length))
        eventToCount(leaf.eventStrings(evtPosition)) = eventToCount.getOrElse(leaf.eventStrings(evtPosition),0) + 1
    }}

    val eventList = eventToCount.toList

    // now make an array of normalized values
    val total = eventList.map{case(key,value) => value}.sum
    val normalizedValues = eventList.map{case(key,count) => count.toDouble / total.toDouble}.toArray[Double]
    val alphabet = eventList.map{case(key,value) => key}.toArray

    eventDistributions :+= NumericalDistribution[String](normalizedValues,alphabet)

    eventCounts :+= eventList.map{case(key,count) => if (key == Event.NONE) 0.0 else count.toDouble}.sum
    eventSites :+= "Site_" + evtPosition
  }}

  // now create a distribution over all sites
  val siteLevelDist = eventCounts.map{case(count) => count.toDouble / eventCounts.sum.toDouble}.toArray[Double]

  /**
   * draw the next event:
   * 1) choose the target site from the cassette
   * 2) choose the event from the known events at that target site
   * all dependent on their relative probabilities for both choices
   *
   */
  def splitOnEvent() : TreeNode = {

    // draw the next event
    val site = NumericalDistribution.drawRandomIndex(siteLevelDist)
    val eventIndex = eventDistributions(site).drawFromDistribution()
    val event = eventDistributions(site).drawFromDistribution()

    // now split the tree into events that contain the event at the specific site (left) and those that don't (right)
    val leftEvents = TreeModel(leaves.filter{leaf => leaf.eventMap contains event},eventSize)
    val rightEvents = TreeModel(leaves.filter{leaf => !(leaf.eventMap contains event)},eventSize)

    return(TreeNode(leftEvents,rightEvents))
  }
}


object TreeModel {

  /**
   * load events from the stats file
   * @param statsFile the stats file
   */
  def createTree(statsFile: File, targetSiteCount: Int) = {
    var events = Array[Event]()

    // ------------------------------------------------------------------------------------
    // load events from the stats file
    // ------------------------------------------------------------------------------------
    Source.fromFile(statsFile).getLines().drop(1).map { line => {
      val eventStrings = line.split("\t").slice(line.split("\t").size - (targetSiteCount * 2), line.split("\t").size - targetSiteCount)
      val isAllNone = eventStrings.map { evt => if (evt == "NONE") 0 else 1 }.sum

      if (isAllNone > 0) {
        val sp = line.split("\t")
        //eventStrings.zipWithIndex.foreach { case (edit, index) => positionToCounts(index)(edit) = positionToCounts(index).getOrElse(edit, 0) + 1 }
        val id = eventStrings.mkString("-")

        // name: String, sample: String, numberOfReads: Int, eventStrings: Array[String], distanceScore: Int = 1
        // if (sample.toInt < 13) 0 else math.floor((sample.toInt - 11) / 2).toInt
        val evt = Event(sp(1), sp(0), sp(4).toInt, eventStrings, sp(0))

      }
    }
    }
    println("loaded " + totalRows + " events")

  }

}