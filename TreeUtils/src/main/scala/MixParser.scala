package main.scala

import java.util

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.io.Source

/**
  * Parse the output from the PHYLIP package MIX and create
  * event strings for each of the internal nodes
  */
class MixParser(mixOutput: String, eventsToNumbers: String, treeToUse: Int) {

  // first parse out the events to number data, and make a look-up table
  val eventToNum = Source.fromFile(eventsToNumbers).getLines()
  val eventToNumHeader = eventToNum.next()

  if (eventToNumHeader != "event\tnumber\tpositions")
    throw new IllegalStateException("Unable to parse out corect header from " + eventsToNumbers + ", saw: " + eventToNumHeader)

  val numberToEvent = new mutable.HashMap[Int,String]()
  val eventToSites = new mutable.HashMap[String,Array[Int]]()

  // now load all the lines
  eventToNum.foreach{evtLine => {
    val sp = evtLine.split("\t")
    if (sp.size == 3) {
      numberToEvent(sp(1).toInt) = sp(0)
      eventToSites(sp(0)) = sp(2).split(",").map { tk => tk.toInt }
    } else {
      numberToEvent(sp(1).toInt) = sp(0)
      eventToSites(sp(0)) = Array[Int](-1)
    }
  }}

  // the header line we're looking for in the file is:
  val headerLine = "From    To     Any Steps?    State at upper node"


  val inputFile = Source.fromFile(mixOutput).getLines()

  var inGenotypeSection = false
  var currentGenotype: Option[Edge] = None
  var currentTreeNumber = 0

  var activeTree : Option[Array[Edge]] = None

  { // scope this so the temp. data structures go away

    // a mapping from the input tree number to it's annotations
    var treeToGenotypes = new mutable.HashMap[Int, ArrayBuffer[Edge]]()

    inputFile.foreach { line => {
      // skip blank lines to make this easier, also skip the weird sub-header line they provide
      if (line != "" && !line.contains("means same as in the node below it on tree")) {
        if (line.startsWith(headerLine)) {
          inGenotypeSection = true
          if (currentTreeNumber == treeToUse) {
            treeToGenotypes(currentTreeNumber) = new ArrayBuffer[Edge]()
          }
          //println("New tree" + currentTreeNumber)
        }
        else if (inGenotypeSection && !line.contains(".")) {
          inGenotypeSection = false
          if (currentGenotype.isDefined) {
            treeToGenotypes(currentTreeNumber) += currentGenotype.get
            println("tree size = " + treeToGenotypes(currentTreeNumber).size)
          }

          currentGenotype = None
          currentTreeNumber += 1
        }
        else if (inGenotypeSection && (line.contains("yes") || line.contains("no"))) {
          if (currentTreeNumber == treeToUse) {
            if (currentGenotype.isDefined) {
              treeToGenotypes(currentTreeNumber) += currentGenotype.get
            }
            val sp = line.stripPrefix(" ").split(" ")
            currentGenotype = Some(Edge(sp(0), sp(1), sp(2) == "yes", currentTreeNumber))
            currentGenotype.get.addChars(sp.slice(3, sp.size).mkString(""))
          }
        }
        else if (inGenotypeSection) {
          if (currentTreeNumber == treeToUse) {
            val sp = line.stripPrefix(" ").split(" ")
            currentGenotype.get.addChars(sp.slice(0, sp.size).mkString(""))
          }
        }
      }

    }
    }

    // close-out the remaining tree
    if (currentGenotype.isDefined) {
      treeToGenotypes(currentTreeNumber) += currentGenotype.get
      println("tree size = " + treeToGenotypes(currentTreeNumber).size)
    }
    currentGenotype = None
    currentTreeNumber += 1

    activeTree = Some(treeToGenotypes(treeToUse).toArray)
  }

  /*
   * helper functions
   */

  def lookupFroms(fromNode: String): List[Edge] = {
    if (!activeTree.isDefined)
      return List[Edge]()
    activeTree.get.filter{case(mp) => mp.from == fromNode}.toList
  }

  def lookupTos(toNode: String): Edge = {
    if (!activeTree.isDefined)
      throw new IllegalStateException("Annotations haven't been loaded!")
    val ret = activeTree.get.filter{case(mp) => mp.to == toNode}.toList
    if (ret.size != 1)
      throw new IllegalStateException("Found " + ret.size + " edges for node " + toNode)
    ret(0)
  }
}

case class Edge(from: String, to: String, changes: Boolean, treeNumber: Int) {
  var chars = new ArrayBuffer[Char]()

  def addChars(inputString: String): Unit = {
    inputString.foreach{char => if (char != ' ') chars += char}
  }

  def toFancyString = from + "," + to + "," + treeNumber + "," + chars.mkString("")
}
