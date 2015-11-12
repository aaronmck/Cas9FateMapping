package main.scala

import main.scala.TreeSimulation.GlobalSiteModel
import org.apache.commons.math3.distribution.PoissonDistribution

import scala.collection.mutable

/**
 * Created by aaronmck on 11/6/15.
 */
case class RootPath(event: Event, timing: PoissonDistribution, model: MulitMarginalDist) {

  var eventsToRevert = 0
  var reverted: Array[Boolean] = event.eventStrings.map{evtStr => if (evtStr == Event.NONE) true else {eventsToRevert += 1; false}}
  val probabilities: Array[Double] = event.eventStrings.map{evtStr =>1.0}

  var time = 1
  while (eventsToRevert < event.size) {
    // draw an event count from our timing distribution
    val revertEvents = timing.sample()

    (0 until revertEvents).foreach { evt => {
      val draw = model.sampleRandomEvent(sampledPositions.keySet.toArray, sampledEvents)
      sampledPositions(draw.position) = true
      eventPairs :+= EventTimePair(draw.position, time, draw.conditionalProb)
    }}
    time += 1
  }

  var probability: Option[Double] = None

  def calculateProbability(avgEvents: Double): Double = probability.getOrElse({
      val ret = eventPairs.map{evt => if (evt == Event.NONE) avgEvents * evt.probability else (1.0 - avgEvents) * evt.probability}.product
      probability = Some(ret)
      ret
    })

  def toEventString(): String = eventPairs.map{str => str.toCustomString()}.mkString(".")
}

case class RevertedSite(site: Int, timePoint: Int, probability: Double) {
  def toCustomString(): String = "site:" + site + ",time:" + timePoint + ",prob:" + probability
}