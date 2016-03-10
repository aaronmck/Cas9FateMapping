package main.scala
import scala.collection.mutable._

/**
 * Created by aaronmck on 12/21/15.
 */
trait InputTable {
  // read in the series of events
  def getUniqueEvents(): Array[IndexedNode]
  def getUniqueMapping(): HashMap[String, String]
  def getAllEvents(): Array[IndexedNode]
  def getTargetSiteCount() : Int
  def getSampleCount() : Map[String,Int]
}
