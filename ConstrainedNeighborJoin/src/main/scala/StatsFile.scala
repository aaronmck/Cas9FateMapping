package main.scala

import java.io.File

import scala.collection.mutable.HashMap
import scala.io.Source

/**
 * load a stats file into a series of events
 */
case class StatsFile(inputFile: File) {

  // read in the series of events
  var events = Array[IndexedNode]()
  var uniqueEvents = new HashMap[String, String]()
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
    val readCount = sp(header("finalF")).toInt + sp(header("finalR")).toInt
    val readProp  = readCount.toDouble / (sp(header("initF")).toInt + sp(header("initR")).toInt).toDouble

    val eventStrings = sp.zipWithIndex.filter{case(tk,index) => eventColumns contains index}.map{case(tk,in) => tk}
    val evt = Event(nodeIndex.toString, sp(0), readCount, readProp, eventStrings)
    nodeIndex += 1
    allEvents :+= evt

    if (!(uniqueEvents contains eventStrings.mkString(""))) {
      val isAllNone = eventStrings.map { evt => if (evt == "NONE") 0 else 1 }.sum
      if (isAllNone > 0) {
        events :+= evt
      }
    }
    uniqueEvents(eventStrings.mkString("")) = "Sample_" + line.split("\t")(0)
  }}

  println("Processed " + events.size + " events")
}
