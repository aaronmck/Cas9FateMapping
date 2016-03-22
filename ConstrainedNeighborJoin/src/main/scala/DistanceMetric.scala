package main.scala

import scala.collection.mutable

/**
 * Created by aaronmck on 11/14/15.
 */
abstract class DistanceMetric[T <: IndexedNode] {
  def distance(event1: T, event2: T): Double
}

/**
 * find the log distance of two events
 * @param logValue the log base to take, assumes 2, but set it to whatever you'd like. Don't be a dick and set it to one
 */
case class SumLogDistance(mutationToCount: EventCounter, noneScore: Double, equalScore: Double, logValue: Double) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = {
    event1.getEventStrings.zip(event2.getEventStrings).map { case (evt1, evt2) => {
      if (evt1 == "NONE" && evt2 == "NONE")
        safeLog2(noneScore)
      else if (evt1 == "NONE" && evt2 != "NONE")
        safeLog2(mutationToCount.getEventCount(event2.getSample,evt2))
      else if (evt1 != "NONE" && evt2 == "NONE")
        safeLog2(mutationToCount.getEventCount(event1.getSample,evt1))
      else if (evt1 != evt2)
        safeLog2(mutationToCount.getEventCount(event1.getSample,evt1) + mutationToCount.getEventCount(event2.getSample,evt2))
      else
        safeLog2(equalScore)
    }
    }.sum / event1.getEventStrings.length.toDouble
  }

  // add one to the value and log-X it
  def safeLog2(count: Double): Double = math.log10(count + 1.0) / math.log10(logValue)
}

/**
 * find the log distance of two events
 * @param logValue the log base to take, assumes 2, but set it to whatever you'd like. Don't be a dick and set it to one
 */
case class AvgLogDistance(mutationToCount: EventCounter, noneScore: Double, equalScore: Double, logValue: Double) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = {
    event1.getEventStrings.zip(event2.getEventStrings).map { case (evt1, evt2) => {
      if (evt1 == "NONE" && evt2 == "NONE")
        safeLog2(noneScore)
      else if (evt1 == "NONE" && evt2 != "NONE")
        safeLog2(mutationToCount.getEventCount(event2.getSample,evt2))
      else if (evt1 != "NONE" && evt2 == "NONE")
        safeLog2(mutationToCount.getEventCount(event1.getSample,evt1))
      else if (evt1 != evt2)
        safeLog2(mutationToCount.getEventCount(event1.getSample,evt1) + mutationToCount.getEventCount(event2.getSample,evt2))
      else
        safeLog2(equalScore)
    }
    }.sum / event1.getEventStrings.length.toDouble
  }

  def safeLog2(count: Double): Double = math.log10(count + 1) / math.log10(logValue)
}

/**
 * find the raw proportion-based distance of two events
 */
case class Hamming(mutationToCount: EventCounter) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = {
    event1.getEventStrings.zip(event2.getEventStrings).map { case (evt1, evt2) => {
      if (evt1 == evt2)
        0
      else
        1
    }
    }.sum / event1.getEventStrings.length.toDouble
  }
}

/**
 * find the raw proportion-based distance of two events
 */
case class SimpleAdditiveDistance(mutationToCount: EventCounter) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = {
    event1.getEventStrings.zip(event2.getEventStrings).map { case (evt1, evt2) => {
      if (evt1 == "NONE" && evt2 == "NONE")
        1
      else if (evt1 == "NONE" && evt2 != "NONE")
        5
      else if (evt1 != "NONE" && evt2 == "NONE")
        5
      else if (evt1 != evt2)
        10
      else
        0
    }
    }.sum
  }
}

/**
  * Cavender-Farris-Neyman corrected distance metric between sequences
  *
  * see page 64 of http://tandy.cs.illinois.edu/textbook.pdf
  *
  * ----DONT USE------ implemented before it was well thought out, the correct isn't appropriate for our data
  *
  * @param mutationToCount just how many mutations we have at each point in the sequence
  */
@deprecated
case class CFNDistance(mutationToCount: EventCounter) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = {
    val hamming = event1.getEventStrings.zip(event2.getEventStrings).map { case (evt1, evt2) => {
      if (evt1 == "NONE" && evt2 == "NONE")
        0
      else if (evt1 == "NONE" && evt2 != "NONE")
        1
      else if (evt1 != "NONE" && evt2 == "NONE")
        1
      else if (evt1 != evt2)
        1
      else
        0
    }
    }

    -0.5 * Math.log(1.0 - (2.0 * (hamming.sum.toDouble / hamming.size.toDouble)))
  }
}