package test.scala


import java.io.PrintWriter

import main.scala._
import org.scalatest.FlatSpec
import org.scalatest.FlatSpec
import org.scalatest.Matchers

import org.scalatest.{Matchers, FlatSpec}

import scala.collection.mutable

/**
 * Created by aaronmck on 11/15/15.
 */
class DistanceMetricTest extends FlatSpec with Matchers {
  "A log2 distance matrix" should "calculate distances correctly" in {
    val events = Array[IndexedNode](FakeEvent("a"),FakeEvent("b"),FakeEvent("c"),FakeEvent("d"),FakeEvent("e"))


    val distances = new DistanceMatrix(events, FakeDistance(1))
    val minSet = distances.minimizeSet(false)
  }
}

case class FakeEvent(name: String) extends IndexedNode {
  var id = -1
  def getID(): Int = id
  def setID(idVal: Int) {id = idVal}
  def merge(otherNode: IndexedNode, branchLeft: Double, branchRight: Double, newID: Int, contrained: Boolean): IndexedNode = {
    FakeEvent("MERGED" + getID + ":" + otherNode.getID)
  }
  def compatible(otherNode: IndexedNode): Boolean = true
  def toFancyString(): String = name + ":" + id
  def newickString(distance: Double, totalDistance: Double, annotationOutput: PrintWriter): String = "NOIDEA"
  def getEventStrings(): Array[String] = Array[String](name)
  def getProportion(): Double = 0.5

  override def getCount(): Int = 100

  override def getSample(): String = "SAMPLE"
}

/**
 * 	a	b	 c	 d	e
 a  0	5	 9	 9	8
 b	5	0	 10 10	9
 c	9	10 0  8  7
 d	9	10 8  0  3
 e	8	9	 7	3	0

 */

/**
 * find the log distance of two events
 */
case class FakeDistance(garbage: Int) extends DistanceMetric[IndexedNode] {

  def distance(event1: IndexedNode, event2: IndexedNode): Double = (event1.getEventStrings()(0),event2.getEventStrings()(0)) match {
    case ("a","b") => 5.0
    case ("a","c") => 9.0
    case ("a","d") => 9.0
    case ("a","e") => 8.0

    case ("b","a") => 5.0
    case ("b","c") => 10.0
    case ("b","d") => 10.0
    case ("b","e") => 9.0

    case ("c","a") => 9.0
    case ("c","b") => 10.0
    case ("c","d") => 8.0
    case ("c","e") => 7.0

    case ("d","a") => 9.0
    case ("d","b") => 10.0
    case ("d","c") => 8.0
    case ("d","e") => 3.0

    case ("e","a") => 8.0
    case ("e","b") => 9.0
    case ("e","c") => 7.0
    case ("e","d") => 3.0

    case (_,_ )=> 0.0
  }

}