package main.scala.TreeSimulation

import scala.util.Random

/**
 * Created by aaronmck on 10/31/15.
 */
case class NumericalDistribution[T](eventDistributions : Array[Double], eventNames : Array[T]) {
  val rand = new Random()

  def drawFromDistribution(): T = {
    return(eventNames(NumericalDistribution.drawRandomIndex(eventDistributions)))
  }
}

object NumericalDistribution {
  def drawRandomIndex(eventDistributions : Array[Double]): Int = {
    val randInst = (new Random()).nextDouble()
    var sumList = 0.0
    eventDistributions.zipWithIndex.foreach{case(lm,index) => {
      sumList += lm
      if (sumList > randInst)
        return(index)
    }}
    eventDistributions.length - 1
  }
}

