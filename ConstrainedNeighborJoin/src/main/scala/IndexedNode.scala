package main.scala

import java.io.PrintWriter

/**
 * Created by aaronmck on 11/19/15.
 */
trait IndexedNode {
  def getID(): Int
  def setID(idVal: Int)
  def merge(otherNode: IndexedNode, branchLeft: Double, branchRight: Double, newID: Int, contrained: Boolean): IndexedNode
  def compatible(otherNode: IndexedNode): Boolean
  def toFancyString(): String
  def newickString(distance: Double, totalDistance: Double, annotationOutput: PrintWriter): String
  def getEventStrings(): Array[String]
  def getProportion(): Double
  def getCount(): Int
  def getSample(): String
}
