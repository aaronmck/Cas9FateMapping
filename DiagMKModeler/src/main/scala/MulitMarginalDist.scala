package main.scala

import java.io.File

import main.scala.TreeSimulation.{GlobalSiteModel, NormalizedCounts}

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.util.Random

/**
 * Created by aaronmck on 11/6/15.
 */
case class MulitMarginalDist(leaves: Array[Event]) {
  val randGenerator = new Random()

  // for each of the target sites, develop a distribution of events within that site.
  // then, aggregate the chance of events happening at each site to a master distribution
  var eventNames = Array[Array[String]]()
  var eventProbs = Array[Array[Double]]()

  // the name of each of our event sites
  var eventSites = Array[String]()

  if (leaves.size < 1)
    throw new IllegalStateException("We require at least 1 event (you should have a lot more really)")

  // for each site in our target array, find the proportion of events at each position
  (0 until leaves(0).size).foreach { evtPosition => {
    val eventToCount = new mutable.HashMap[String, Int]()

    // count the occurrences of each event type at this position
    // if we have multisite deletions only sample them at the rate of 1/(event target length)
    leaves.foreach { leaf => {
        eventToCount(leaf.eventStrings(evtPosition)) = eventToCount.getOrElse(leaf.eventStrings(evtPosition), 0) + 1
    }}

    val eventList = eventToCount.toList

    // now make an array of normalized values
    val total = eventList.map { case (key, value) => value }.sum
    eventProbs :+= eventList.map { case (key, count) => count.toDouble / total.toDouble }.toArray[Double]
    eventNames :+= eventList.map { case (key, value) => key }.toArray

    eventSites :+= "Site_" + evtPosition

  }}

  def sampleRandomEvent(notThesePositions: Array[Int], notTheseEvents: Array[String]): SampledEvent = {
    var newSum = 0.0
    eventProbs.zipWithIndex.foreach { case (column, index) =>
      if (!(notThesePositions contains column))
        newSum += column.zipWithIndex.map { case (prob, ind2) => if (!(notTheseEvents contains eventNames(index)(ind2))) prob else 0.0 }.sum
    }

    // ok, now that we have a marginal sum, find a random
    val randvalue = randGenerator.nextDouble()
    var runningTotal = 0.0

    // in case we run off the end
    var runOffEndCol = 0
    var eventNameEnd = ""
    var probEnd = 0.0

    eventProbs.zipWithIndex.foreach { case (column, index) =>
      if (!(notThesePositions contains column)) {
        runOffEndCol = index
        column.zipWithIndex.map { case (prob, ind2) =>
          if (!(notTheseEvents contains eventNames(index)(ind2))) {
            eventNameEnd = eventNames(index)(ind2)
            probEnd = prob / newSum
            if ((runningTotal + prob) / newSum > randvalue)
              return (SampledEvent(runOffEndCol, eventNameEnd, probEnd))
            runningTotal += prob
        }}
      }
    }
    return (SampledEvent(runOffEndCol, eventNameEnd, probEnd))

  }

  def revertRandomEvent(notThesePositions: Array[Int], notTheseEvents: Array[String]): SampledEvent = {
    var newSum = 0.0
    eventProbs.zipWithIndex.foreach { case (column, index) =>
      if (!(notThesePositions contains column))
        newSum += column.zipWithIndex.map { case (prob, ind2) => if (!(notTheseEvents contains eventNames(index)(ind2))) prob else 0.0 }.sum
    }

    // ok, now that we have a marginal sum, find a random
    val randvalue = randGenerator.nextDouble()
    var runningTotal = 0.0

    // in case we run off the end
    var runOffEndCol = 0
    var eventNameEnd = ""
    var probEnd = 0.0

    eventProbs.zipWithIndex.foreach { case (column, index) =>
      if (!(notThesePositions contains column)) {
        runOffEndCol = index
        column.zipWithIndex.map { case (prob, ind2) =>
          if (!(notTheseEvents contains eventNames(index)(ind2))) {
            eventNameEnd = eventNames(index)(ind2)
            probEnd = prob / newSum
            if ((runningTotal + prob) / newSum > randvalue)
              return (SampledEvent(runOffEndCol, eventNameEnd, probEnd))
            runningTotal += prob
          }}
      }
    }
    return (SampledEvent(runOffEndCol, eventNameEnd, probEnd))

  }
}

case class SampledEvent(position: Int, eventString: String, conditionalProb: Double)

/**
 * our static functions for the TreeModel class
 */
object MulitMarginalDist {

  /**
   * load events from the stats file and create a TreeModel
   *
   * @param statsFile the stats file
   */
  def createMode(statsFile: File, targetSiteCount: Int): Tuple2[MulitMarginalDist, Array[Event]] = {
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
    (MulitMarginalDist(events.toArray), events.toArray)
  }
}