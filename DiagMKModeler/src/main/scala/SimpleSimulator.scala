package main.scala

import org.apache.commons.math3.distribution.PoissonDistribution
import scala.collection.mutable.ArrayBuffer
import akka.actor.Actor
import akka.event.Logging

/**
 * Created by aaronmck on 11/11/15.
 */
class SimpleSimulator(eventsRatios: Array[Array[Double]],
                      events: Array[Event],
                      estimatedGeneration: Int,
                      estimateRate: Double) extends Actor {

  // setup an array of events at the estimated event size we'll have; we can expand it
  // later if need-be, but start with this
  val nodes = new Array[Node](math.pow(2,estimatedGenerations).toInt)

  // where are we right now
  var currentGeneration = 0

  // a logger
  val log = Logging(context.system, this)

  // get a Poisson distribution
  val rateDist = new PoissonDistribution(estimateRate)

  def receive = {
    case "simulate" => {
      log.info("received test")

      //

    }
    case _ => log.info("received unknown message")
  }
}

case class Node(parent: Int,
                left: Option[Node] = None,
                right: Option[Node] = None,
                events:Array[Int] = new Array[Int](10)) {

  def setEvent(position: Int, event: Int): Unit = {
    if (position < 0 || position >= events.length)
      throw new IllegalStateException("Unable to set position " + position)
    events(position) = event
  }
}