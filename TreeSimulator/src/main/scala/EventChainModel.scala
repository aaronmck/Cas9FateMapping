package main.scala

import org.apache.commons.math3.distribution.PoissonDistribution

/**
 * Created by aaronmck on 11/6/15.
 */
class EventChainModel(events: Array[Event], mutationMu: Double) {
  // our timing distribution
  val timing = new PoissonDistribution(mutationMu)

  // for each event, draw out a mutation chain given the mu
  events.map{event => {
    // draw from our timing distribution the number of events to revert
    val revertEvents = timing.sample()
  }}
}
