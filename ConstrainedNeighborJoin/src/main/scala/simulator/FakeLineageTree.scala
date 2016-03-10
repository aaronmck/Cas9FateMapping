package main.scala.simulator

import org.apache.commons.math3.distribution.ExponentialDistribution

import scala.collection.mutable
import scala.util.Random

/**
 * Created by aaronmck on 12/14/15.
 */
case class Node(events: Array[Int], depth: Int, eventDist: Double) {
  val expo = new ExponentialDistribution(eventDist)

  def splitNode(targetDepth: Int, mutationRate: Double,nextMutNumber: Int): Tuple3[Node,Node,Int] = {
    var newLastNumber = nextMutNumber
    val rando = new Random()

    val childLeftEvents = events
    val childRightEvents = events

    // for each event site, figure out if we should mutate it
    events.zipWithIndex.foreach{case(evt,index) => {

      // sample from the distribution: does the event happen?
      if (rando.nextDouble() <= mutationRate) {
        val length = math.max(1,expo.sample(1)(0).ceil.toInt)
        (0 until length).foreach{ind =>
          if (ind + index < childLeftEvents.length && childLeftEvents(ind + index) == 0) childLeftEvents(ind + index) = newLastNumber
        }
        newLastNumber += 1
      } else {
        evt
      }
    }}

    events.zipWithIndex.foreach{case(evt,index) => {
      if (rando.nextDouble() <= mutationRate) {
        val length = math.max(1,expo.sample(1)(0).ceil.toInt)
        (0 until length).foreach{ind =>
          if (ind + index < childRightEvents.length && childLeftEvents(ind + index) == 0) childRightEvents(ind + index) = newLastNumber}
        newLastNumber += 1
      } else {
        evt
      }
    }}

    val childLeft  = Node(childLeftEvents,  depth + 1, eventDist)
    val childRight = Node(childRightEvents, depth + 1, eventDist)
    return (childLeft,childRight,newLastNumber)
  }
}



class FakeLineageTree(depth: Int, siteCount: Int, mutationRate: Double, eventDist: Double) {

  val root = Node(new Array[Int](siteCount),0, eventDist)
  var currentGeneration = new mutable.ArrayBuffer[Node]()
  currentGeneration += root
  var currentMutationIndex = 1

  (0 until depth).foreach{dpt => {
    var newGeneration = new mutable.ArrayBuffer[Node]()
    currentGeneration.toArray.foreach{node => {
      val children = node.splitNode(depth,mutationRate,currentMutationIndex)
      newGeneration += children._1
      newGeneration += children._2
      currentMutationIndex = children._3
    }}
    currentGeneration = newGeneration
  }}
  //println(currentGeneration.toArray.size)

  // now find out how many nodes in the final generation are edited
  var editCount = 0
  var sitesCount = 0
  currentGeneration.foreach{node => {
    node.events.foreach{evt => {
      sitesCount += 1
      if (evt != 0) editCount += 1
    }}
  }}

  val editRate = editCount.toDouble / sitesCount.toDouble
}
