package main.scala

import java.io.File

import scala.collection.mutable
import scala.collection.mutable.{Map, HashMap}
import scala.io.Source

/**
 * load a stats file into a series of events
 */
case class StatsFile(inputFile: File, targetSiteCount: Int = 10, filterByOccurance: Int = 10) extends InputTable {

  // read in the series of events
  var uniqueEventsMapping = new HashMap[String, IndexedNode]()

  var allEvents = Array[IndexedNode]()

  // sample counts
  val sampleToCount = new HashMap[String,Int]()

  // read in our lines
  val inputLines = Source.fromFile(inputFile).getLines()

  // find the target columns
  val header = inputLines.next().split("\t").zipWithIndex.toMap

  // find the columns that map to target outputs
  var eventColumns = header.filter{case(token,index) => token.startsWith("target")}.map{case(tk,index) => index}.toArray

  var nodeIndex = 0
  var droppedLine = 0
  
  // now make the event string from the target columns for each remaining line
  inputLines.foreach{line => {

    if ((line contains "PASS") && !(line contains "UNKNOWN") && !(line contains "WT")) {
      val sp = line.split("\t")
      val eventStrings = sp.zipWithIndex.filter { case (tk, index) => eventColumns contains index }.map { case (tk, in) => tk }
      val evt = Event(nodeIndex.toString, sp(0), 1, eventStrings)
      nodeIndex += 1
      allEvents :+= evt

      // add this event to our per-sample count
      sampleToCount(evt.getSample()) = sampleToCount.getOrElse(evt.getSample(), 0) + 1

      // a unique key of the sample and events -- we want unique events per sample
      val keyString = sp(0) + "_" + eventStrings.mkString("") // -- if we want to go back to sample specific
      //val keyString = eventStrings.mkString("")

      if (!(uniqueEventsMapping contains keyString)) {
        val isAllNone = eventStrings.map { evt => if (evt == "NONE") 0 else 1 }.sum
        if (isAllNone > 0) {
          uniqueEventsMapping(keyString) = evt
        }
      } else {
        uniqueEventsMapping(keyString).addSupport(1)
      }
    } else {
      droppedLine += 1
    }
  }}

  println("Dropped " + droppedLine + " lines")
  val uniqueEvents = uniqueEventsMapping.values.filter{event => event.getSupport > filterByOccurance}.toArray
  println("unique event size " + uniqueEvents.size)
  val uniqueEventsMappingStrings = uniqueEventsMapping.map{case(samplePlusEvent,event) => (samplePlusEvent,event.getEventStrings().mkString(""))}

  println("Processed " + allEvents.size + " events")

  // read in the series of events
  override def getUniqueEvents(): Array[IndexedNode] = uniqueEvents

  override def getUniqueMapping(): mutable.HashMap[String, String] = uniqueEventsMappingStrings

  override def getTargetSiteCount(): Int = targetSiteCount

  override def getAllEvents(): Array[IndexedNode] = allEvents

  override def getSampleCount() : Map[String,Int] = sampleToCount
}
