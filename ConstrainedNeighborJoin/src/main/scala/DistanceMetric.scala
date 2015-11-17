package main.scala

import scala.collection.mutable

/**
 * Created by aaronmck on 11/14/15.
 */
abstract class DistanceMetric(mutationToCount: mutable.HashMap[String,Int], noneScore: Double) {
  def distance(event1: Event, event2: Event): Double
}

/**
 * find the log distance of two events
 * @param logValue the log base to take, assumes 2, but set it to whatever you'd like. Don't be a dick and set it to one
 */
case class SumLogDistance(mutationToCount: mutable.HashMap[String,Int], noneScore: Double, logValue: Double) extends DistanceMetric(mutationToCount,noneScore) {

  def distance(event1: Event, event2: Event): Double = {
    event1.eventStrings.zip(event2.eventStrings).map { case (evt1, evt2) => {
      if (evt1 == "NONE" && evt2 == "NONE")
        noneScore
      else if (evt1 == "NONE" && evt2 != "NONE")
        safeLog2(mutationToCount(evt2))
      else if (evt1 != "NONE" && evt2 == "NONE")
        safeLog2(mutationToCount(evt1))
      else if (evt1 != evt2)
        safeLog2(mutationToCount(evt1) + mutationToCount(evt2))
      else
        0
    }
    }.sum / event1.eventStrings.length.toDouble
  }

  def safeLog2(count: Int): Double = math.log10(count + 1) / math.log10(logValue)
}