package main.scala.TreeSimulation

import scala.util.Random

/**
 * Created by aaronmck on 10/31/15.
 */
case class NormalizedCounts[T](eventDistributions : Array[Double], eventNames : Array[T]) {
  val rand = new Random()

  def drawFromDistribution(): Tuple2[T,Double] = {
    val index = NormalizedCounts.drawRandomIndex(eventDistributions)
    //val index = NumericalDistribution.drawMaxIndex(eventDistributions)
    return (eventNames(index),eventDistributions(index))
  }
}

object NormalizedCounts {
  /**
   *
   * @param eventDistributions
   * @return
   */
  def drawDistEvent(eventDistributions : Array[Double]): Int = {
    if (eventDistributions.length == 0)
      throw new IllegalAccessError("Unable to sample from array of size 0")

    if (eventDistributions.length == 1)
      return 0

    val randInst = (new Random()).nextDouble()
    var sumList = 0.0
    eventDistributions.zipWithIndex.foreach{case(lm,index) => {
      sumList += lm
      if (sumList > randInst)
        return(index)
    }}
    eventDistributions.length - 1
  }

  /**
   *
   * @param eventDistributions
   * @return
   */
  def drawRandomIndex(eventDistributions : Array[Double]): Int = {
    if (eventDistributions.length == 0)
      throw new IllegalAccessError("Unable to sample from array of size 0")

    if (eventDistributions.length == 1)
      return 0

    val randInst = (new Random()).nextDouble()
    var sumList = 0.0
    eventDistributions.zipWithIndex.foreach{case(lm,index) => {
      sumList += lm
      if (sumList > randInst)
        return(index)
    }}
    eventDistributions.length - 1
  }

  /**
   *
   * @param eventDistributions
   * @return
   */
  def drawMaxIndex(eventDistributions : Array[Double]): Int = {
    if (eventDistributions.length == 0)
      throw new IllegalAccessError("Unable to sample from array of size 0")

    if (eventDistributions.length == 1)
      return 0

    var maxValue = 0.0
    var ret = 0
    eventDistributions.zipWithIndex.foreach{case(lm,index) => {
      if (lm > maxValue) {
        maxValue = lm
        ret = index
      }
    }}
    return ret
  }
}

