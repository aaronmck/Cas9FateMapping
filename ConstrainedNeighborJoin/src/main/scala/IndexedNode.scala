package main.scala

import java.io.PrintWriter

/**
 * generic events (nodes) for our distance calculations
 */
trait IndexedNode {
  def getID(): Int
  def setID(idVal: Int)
  def merge(otherNode: IndexedNode, branchLeft: Double, branchRight: Double, newID: Int, contrained: Boolean): IndexedNode
  def compatible(otherNode: IndexedNode): Boolean
  def toFancyString(): String
  def newickString(distance: Double, totalDistance: Double, annotationOutput: PrintWriter): String
  def getEventStrings(): Array[String]
  def getSupport(): Int
  def getSample(): String
  def addSupport(additionalSupportCount: Int)
  def countNonWT(): Int
}
