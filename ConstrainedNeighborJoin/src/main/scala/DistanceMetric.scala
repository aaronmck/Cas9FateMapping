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
case class SumLogDistance(mutationToCount: EventCounter, noneScore: Double, logValue: Double) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = {
    event1.getEventStrings.zip(event2.getEventStrings).map { case (evt1, evt2) => {
      if (evt1 == "NONE" && evt2 == "NONE")
        noneScore
      else if (evt1 == "NONE" && evt2 != "NONE")
        safeLog2(mutationToCount.getEventCount(event2.getSample,evt2))
      else if (evt1 != "NONE" && evt2 == "NONE")
        safeLog2(mutationToCount.getEventCount(event1.getSample,evt1))
      else if (evt1 != evt2)
        safeLog2(mutationToCount.getEventCount(event1.getSample,evt1) + mutationToCount.getEventCount(event2.getSample,evt2))
      else
        0
    }
    }.sum / event1.getEventStrings.length.toDouble
  }

  def safeLog2(count: Double): Double = math.log10(count + 1) / math.log10(logValue)
}

/**
 * find the raw proportion-based distance of two events
 */
case class RawDistance(mutationToCount: EventCounter, noneScore: Double) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = {
    event1.getEventStrings.zip(event2.getEventStrings).map { case (evt1, evt2) => {
      if (evt1 == "NONE" && evt2 == "NONE")
        noneScore
      else if (evt1 == "NONE" && evt2 != "NONE")
        mutationToCount.getEventCount(event2.getSample,evt2)
      else if (evt1 != "NONE" && evt2 == "NONE")
        mutationToCount.getEventCount(event1.getSample,evt1)
      else if (evt1 != evt2)
        mutationToCount.getEventCount(event1.getSample,evt1) + mutationToCount.getEventCount(event2.getSample,evt2)
      else
        0
    }
    }.sum / event1.getEventStrings.length.toDouble
  }
}