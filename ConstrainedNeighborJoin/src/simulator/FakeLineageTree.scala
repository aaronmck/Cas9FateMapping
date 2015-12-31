package simulator

import org.apache.commons.math3.distribution.ExponentialDistribution

import scala.collection.mutable
import scala.util.Random

/**
 * Created by aaronmck on 12/14/15.
 */
case class Node(parent: Option[Node], events: Array[Int], depth: Int) {
  var childLeft: Option[Node] = None
  var childRight: Option[Node] = None
  val expo = new ExponentialDistribution(1.0)

  def splitNode(targetDepth: Int, mutationRate: Double,nextMutNumber: Int): Tuple3[Node,Node,Int] = {
    var newLastNumber = nextMutNumber
    val rando = new Random()

    val childLeftEvents = new Array[Int](events.length)
    val childRightEvents = new Array[Int](events.length)

    events.zipWithIndex.foreach{case(evt,index) => {
      if (rando.nextDouble() <= mutationRate) {
        val length = expo.sample(1)(0).ceil.toInt
        (0 until length).foreach{ind => if (ind + index < childLeftEvents.length) childLeftEvents(ind + index) = newLastNumber}
        newLastNumber += 1
      } else {
        evt
      }
    }}

    events.zipWithIndex.foreach{case(evt,index) => {
      if (rando.nextDouble() <= mutationRate) {
        val length = expo.sample(1)(0).ceil.toInt
        (0 until length).foreach{ind => if (ind + index < childRightEvents.length) childRightEvents(ind + index) = newLastNumber}
        newLastNumber += 1
      } else {
        evt
      }
    }}

    childLeft  = Some(Node(Some(this), childLeftEvents,  depth + 1))
    childRight = Some(Node(Some(this), childRightEvents, depth + 1))
    return (childLeft.get,childRight.get,newLastNumber)
  }
}



class FakeLineageTree(depth: Int, siteCount: Int, mutationRate: Double) {

  val root = Node(None,new Array[Int](siteCount),0)
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

  // now use the final generation to estimate the distances between all nodes... we know their distances

}
