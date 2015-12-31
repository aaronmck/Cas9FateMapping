package main.scala

import java.io.File

import scala.collection.mutable
import scala.collection.mutable.HashMap
import scala.io.Source

/**
 * load a stats file into a series of events
 */
case class StatsFile(inputFile: File) extends InputTable {

  // read in the series of events
  var uniqueEventsMapping = new HashMap[String, IndexedNode]()
  var allEvents = Array[IndexedNode]()
  val targetSiteCount = 10

  // read in our lines
  val inputLines = Source.fromFile(inputFile).getLines()

  // find the target columns
  val header = inputLines.next().split("\t").zipWithIndex.toMap

  // find the columns that map to target outputs
  var eventColumns = header.filter{case(token,index) => token.startsWith("target")}.map{case(tk,index) => index}.toArray

  var nodeIndex = 0

  // now make the event string from the target columns for each remaining line
  inputLines.foreach{line => {
    //println(line)
    val sp = line.split("\t")
    val eventStrings = sp.zipWithIndex.filter{case(tk,index) => eventColumns contains index}.map{case(tk,in) => tk}
    val evt = Event(nodeIndex.toString, sp(0), 1, eventStrings)
    nodeIndex += 1
    allEvents :+= evt

    val keyString = sp(0) + "_" + eventStrings.mkString("")

    if (!(uniqueEventsMapping contains keyString)) {
      // println("non-duplicate event: " + keyString)
      val isAllNone = eventStrings.map { evt => if (evt == "NONE") 0 else 1 }.sum
      if (isAllNone > 0) {
        //println("all is none: " + isAllNone)
        uniqueEventsMapping(keyString) = evt
      }
    } else {
      // println("duplicate event: " + keyString)
      uniqueEventsMapping(keyString).addSupport(1)
    }

  }}
  val uniqueEvents = uniqueEventsMapping.values.toArray
  val uniqueEventsMappingStrings = uniqueEventsMapping.map{case(samplePlusEvent,event) => (samplePlusEvent,event.getEventStrings().mkString(""))}

  println("Processed " + allEvents.size + " events")

  // read in the series of events
  override def getUniqueEvents(): Array[IndexedNode] = uniqueEvents

  override def getUniqueMapping(): mutable.HashMap[String, String] = uniqueEventsMappingStrings

  override def getTargetSiteCount(): Int = targetSiteCount

  override def getAllEvents(): Array[IndexedNode] = allEvents
}
