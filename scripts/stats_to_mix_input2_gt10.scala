// --------------------------------------------------
import scala.io._
import java.io._
import scala.collection.mutable._
import scala.collection.immutable.ListMap

// --------------------------------------------------
// load up the lines from our stats/calls file
val tearsheet = Source.fromFile(args(0)).getLines().drop(1)

val output = new PrintWriter(args(1))
val mixCommandFile = new PrintWriter(args(2))
val annotationFile = new PrintWriter(args(3))
val weightFile = new PrintWriter(args(4))
val eventsToNumbers = new PrintWriter(args(5))
val filteringNumber = 5

// --------------------------------------------------
// all reads object
// --------------------------------------------------
case class Event(events: Array[String], eventNumbers: Array[Int], count: Int, proportion: Double, sample: String, name: String) {
  val builder = new StringBuilder()

  // make a padded version of the name
  val paddedName = name + (0 until (10 - name.length)).map{i => " "}.mkString("")

  // make a binary representation of the events
  def toMixString(index: Int): String = {
    paddedName + (1 until index).map{ind => {
      if (eventNumbers contains ind)
        "1"
      else "0"
    }}.mkString("")
  }
}
val eventsBuffer = new ArrayBuffer[Event]()

// ---------------------------------------------------------
// load up the event columns into arrays of numeric events
// ---------------------------------------------------------
val eventToNumber = new HashMap[String,Int]()
val eventToCount = new HashMap[String,Int]()
val numberToEvent = new HashMap[Int,String]()
val eventToPositions = new HashMap[String,Set[Int]]()

eventToNumber("NONE") = 0
eventToNumber("UNKNOWN") = 0
numberToEvent(0) = "NONE"

var linesProcessed = 0
var nextIndex = 1

// ---------------------------------------------------------
// find the event columns -- make a mapping
// from the target number to the column index, assuming they're
// labeled targetX where is X is 1,2,3,4....
// ---------------------------------------------------------
tearsheet.foreach{line => {
  val sp = line.split("\t")
  Source.fromFile(sp(3) + "/" + sp(0) + "/" + sp(0) + ".allReadCounts").getLines().drop(1).zipWithIndex.foreach{case(line,index) => {
    val lineTks = line.split("\t")

    val eventTokens = lineTks(0).split("_").toArray

    // count the events
    eventTokens.foreach{case(event) => eventToCount(event) = eventToCount.getOrElse(event,0) + lineTks(2).toInt}

    eventTokens.zipWithIndex.foreach{case(evt,index) => {
      eventToPositions(evt) = eventToPositions.getOrElse(evt,Set[Int]()) + index
    }}

    val eventNumbers = eventTokens.map{evt => {
      if (eventToNumber contains evt) {
        eventToNumber(evt)
      } else {
        eventToNumber(evt) = nextIndex
        numberToEvent(nextIndex) = evt 
        nextIndex += 1
        eventToNumber(evt)
      }
    }}


    val evt = Event(lineTks(0).split("_").toArray, eventNumbers, lineTks(2).toInt, lineTks(3).toDouble, sp(0), "N" + linesProcessed)

    if (evt.count >= filteringNumber)
      eventsBuffer += evt
    linesProcessed += 1
  }}
}}

val keyToTag = new HashMap[String,String]()
val keyToCount = new HashMap[String,Int]()


// scale the values from counts to characters (0-9, A-Z)
val characterArray = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"

// normalize the weights to the range of values we have
val maxCount = eventToCount.values.max

// this function scales the observed counts to weights accepted by PHYLIP 
def scaleValues(value: Int, min: Int, max: Int): Char = {
  val maxLog = math.log(max)
  val valueLog = math.log(value)
  println(maxLog + " " + valueLog + " ")
  val ret = scala.math.round(((valueLog.toDouble - min.toDouble) / maxLog.toDouble) * (characterArray.length.toDouble -1.0)).toInt
  characterArray(ret)
}

// -----------------------------------------------------------------------------------
// map the each of the events to associated weights in PHYLIP space and write to disk
// -----------------------------------------------------------------------------------
val weights = (1 until nextIndex).map{
  case(index) => {
    scaleValues(eventToCount(numberToEvent(index)),0,maxCount)
  }
}.toArray

weightFile.write(weights.mkString("") + "\n")
weightFile.close()

// -----------------------------------------------------------------------------------
// for each event, output the name and events
// -----------------------------------------------------------------------------------
output.write((eventsBuffer.toArray.size + 1) + "\t" + (nextIndex - 3) + "\n")

// (events: Array[String], eventNumbers: Array[String], count: Int, proportion: Float, sample: String, name: String)
annotationFile.write("taxa\tsample\tcount\tproportion\tevent\n")
output.write("root      " + (1 until nextIndex).map{index => "0"}.mkString("") + "\n")

eventsBuffer.toArray.foreach{evt => {
  output.write(evt.toMixString(nextIndex) + "\n")
  annotationFile.write(evt.name + "\t" + evt.sample + "\t" + evt.count + "\t" + evt.proportion + "\t" + evt.events.mkString("_") + "\n")
}}

annotationFile.close()
output.close()

eventsToNumbers.write("event\tnumber\tpositions\n")
eventToNumber.foreach{case(event,number) => {
  eventToPos = eventToPositions.getOrElse(event,Set[Int]()).mkString(",")
  eventsToNumbers.write(event + "\t" + number + "\t" + eventToPos + "\n")
}}
eventsToNumbers.close()
