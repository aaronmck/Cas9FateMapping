// --------------------------------------------------
import scala.io._
import java.io._
import scala.collection.mutable._
import scala.collection.immutable.ListMap


// --------------------------------------------------
// load up the lines from our stats/calls file
val lines = Source.fromFile(args(0)).getLines()
val header = lines.next()

val output = new PrintWriter(args(1))
val mixCommandFile = new PrintWriter(args(2))
val annotationFile = new PrintWriter(args(3))



// --------------------------------------------------
// find the event columns -- make a mapping
// from the target number to the column index, assuming they're
// labeled targetX where is X is 1,2,3,4....
// --------------------------------------------------

val targetToIndex = new LinkedHashMap[Int,Int]()
header.split("\t").zipWithIndex.foreach{case(tk,index) => {
  if (tk contains "target") {
    val col = tk.stripPrefix("target").toInt
    targetToIndex(col) = index
    println("target" + col + " position " + index)
  }
}}

// --------------------------------------------------
// load up the event columns into arrays of numeric events
// --------------------------------------------------
val eventToNumber = new HashMap[String,Int]()
val eventToCount = new HashMap[String,Int]()
val numberToEvent = new HashMap[Int,String]()

eventToNumber("NONE") = 0
eventToNumber("UNKNOWN") = 0
numberToEvent(0) = "NONE"

var nextIndex = 1

// store a mapping from the column to a list of event numbers
val columnToEventsArray = new HashMap[Int,ArrayBuilder[Int]]()
targetToIndex.foreach{case(target,column) => columnToEventsArray(target) = ArrayBuilder.make[Int]}

var linesProcessed = 0

// -------------------------------------------------------------------
//  for each line, for each target assign a number
// -------------------------------------------------------------------
lines.foreach{line => {
  linesProcessed += 1
  if ((line contains "PASS") && !((line contains "UNKNOWN") || (line contains "WT"))) {

    val sp = line.split("\t")
    targetToIndex.foreach{case(target,index) => {
      val event = sp(index)

      if (event != "UNKNOWN") { // event != "NONE" &&
        if (!(eventToNumber contains event)) {
          eventToNumber(event) = nextIndex
          numberToEvent(nextIndex) = event
          nextIndex += 1
        }
        eventToCount(event) = eventToCount.getOrElse(event,0) + 1
        columnToEventsArray(target) = columnToEventsArray.getOrElse(target,ArrayBuilder.make[Int]) += eventToNumber(event)
      }
    }}
  }
}}

// -------------------------------------------------------------------
//  now for each line, for each target assign a number
// -------------------------------------------------------------------
val lines2 = Source.fromFile(args(0)).getLines()
val header2 = lines2.next()
var index = 1

output.write(linesProcessed + "\t" + (nextIndex - 3) + "\n")
annotationFile.write("taxa\teventString\tnumberStringg\n")

lines2.foreach{line => {
  if ((line contains "PASS") && !((line contains "UNKNOWN") || (line contains "WT"))) {
    val sp = line.split("\t")
    val builder = new StringBuilder()
    val name = ("N" + index)
    index += 1
    builder ++= name + (0 until (10 - name.length)).map{i => " "}.mkString("")

    val tags = new HashMap[Int,Boolean]()
    var eventString = Array[String]()
    var numberString = Array[String]()
    targetToIndex.foreach{case(target,index) => {
      tags(eventToNumber(sp(index))) = true
      eventString :+= sp(index)
      numberString :+= eventToNumber(sp(index)).toString.reverse.padTo(5,' ').reverse
    }}

    annotationFile.write(name + "\t" + eventString.mkString("-") + "\t" + numberString.mkString("-") + "\n")

    (1 until nextIndex).foreach{index => builder ++= {if (tags contains index) "1" else "0"}}

    output.write(builder.result + "\n")
  }

}}
output.close()
