// load up a stats file and produce a distance matrix using a simple scoring scheme
import scala.io._
import java.io._
import scala.collection.mutable._

// args:
// 1: input file
// 2: the output file
// 3: the number of events to look for
val output = new PrintWriter(args(1))
val eventCount = args(2).toInt

// ------------------------------------------------------------------------------------
// our event helper class - stores a single line of our UMI stats table
// ------------------------------------------------------------------------------------
case class Event(name: String, eventStrings: Array[String], line: Int) {

  def probability(otherEvent:Event, occuranceProportion: HashMap[Int,HashMap[String,Double]]) : Double =
    if (line == otherEvent.line)
      return 0.0
    else
      eventStrings.zip(otherEvent.eventStrings).zipWithIndex.map{case((event1,event2),index) => {
        Event.probJointEvents(event1,event2,occuranceProportion(index))
      }}.product

}

object Event {
  // calculate the basic distance between two points
  def probJointEvents(event1: String, event2: String, occuranceProportion: HashMap[String,Double], epsilon: Double = 0.001 /* completely arbritrary */): Double = {
      // both not edited or both the same edit
      if (event1 == "NONE" && event2 == "NONE")
          1.0
      else if (event1 != "NONE" && event2 == event1)
          5.0
      else if (event1 == "NONE" && event2 != event1)
          0.0
      else if (event2 == "NONE" && event2 != event1)
          0.0
      else // the events should have a distance inversely proportional to occurance
          epsilon
        }
}


// ------------------------------------------------------------------------------------
// Main script
// ------------------------------------------------------------------------------------


var events = Array[Event]()
val positionToCounts = HashMap[Int,HashMap[String,Int]]()
val positionToCountsNormalized = HashMap[Int,HashMap[String,Double]]()
var totalRows = 0
(0 until eventCount).foreach{i => {
  positionToCounts(i) = new HashMap[String,Int]()
  positionToCountsNormalized(i) = new HashMap[String,Double]()
}}


// strip down a stats file to a set of events, and record the number of events per site
Source.fromFile(args(0)).getLines().drop(1).zipWithIndex.foreach{case(line,index) => {
  val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount * 2),line.split("\t").size - eventCount)
  val isAllNone = eventStrings.map{evt => if (evt == "NONE") 0 else 1}.sum

  if (isAllNone > 0) {
    eventStrings.zipWithIndex.foreach{case(edit,index) => positionToCounts(index)(edit) = positionToCounts(index).getOrElse(edit,0) + 1}
    events :+= Event(line.split("\t")(0) + "_" + line.split("\t")(1),eventStrings, index)
    totalRows += 1
  }
}}
println("loaded " + totalRows + " events")


// now normalize the positional counts by the total
positionToCounts.foreach{case(position,counts) => {
  counts.foreach{case(event,count) => positionToCountsNormalized(position)(event) = count.toDouble / totalRows.toDouble}}
}

// should we make a event list
var processedEvents2 = new HashMap[String,Boolean]()
events.foreach{event1 => {

  if (!(processedEvents2 contains event1.eventStrings.mkString(""))) {
    processedEvents2(event1.eventStrings.mkString("")) = true
  }
}}

var processedEvents = new HashMap[String,Boolean]()

// output.write("site\t" + events.map{evt => evt.name}.mkString("\t") + "\n")
output.write("\t" + processedEvents2.size + "\n")
events.foreach{event1 => {

  if (!(processedEvents contains event1.eventStrings.mkString(""))) {
    processedEvents(event1.eventStrings.mkString("")) = true
    output.write(event1.name.slice(0,10).padTo(10, " ").mkString(""))
    output.write(events.map{event2 => "%.4f" format (event1.probability(event2,positionToCountsNormalized)*100.0)}.mkString(" ") + "\n")
  }
}}
output.close()
