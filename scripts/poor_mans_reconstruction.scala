// load up a stats file and produce a distance matrix using a simple scoring scheme
import scala.io._
import java.io._
import scala.collection.mutable._

// args:
// 1: input file
// 2: the number of events to look for
val eventCount = args(1).toInt

// ------------------------------------------------------------------------------------
// our event helper class - stores a single line of our UMI stats table
// ------------------------------------------------------------------------------------
case class Event(name: String, sample: String, numberOfReads: Int, eventStrings: Array[String], distanceScore: Int = 1) {

  // make a map that contains the actual events and their positions as a way to remove duplicate - multisite events
  val eventMap = new HashMap[String,Array[Int]]()
  eventStrings.zipWithIndex.foreach{case(evt,index) => eventMap(evt) = eventMap.getOrElse(evt,Array[Int]()) :+ index}
  val containsNone = if (eventMap contains "NONE") 1 else 0

  // store our size for quick reference
  val uniqueSize = eventMap.size - containsNone

  // if and when we find a parent record it
  var probParent: Option[Array[Event]] = None

  def knownParent = if(sample.toInt < 13) 0 else math.floor((sample.toInt - 11) / 2).toInt

  // generate some output
  def toOutputString(): String = {
    val oth = if (probParent.isDefined) {
      probParent.get.map{parent => parent.eventStrings.mkString("-")}.mkString(";") + "\t" + probParent.get.map{parent => parent.sample}.mkString(";") + "\t" + knownParent
    } else "NONE\tNONE\tNONE"

    name + "\t" + uniqueSize + "\t" + sample + "\t" + numberOfReads + "\t" + eventStrings.mkString(",") + "\t" + oth
  }

  // given an array of events that have an edit distance one less than ours, find a likely parent
  def oneEdit(events: Array[Event]): Boolean = {
    var found = false
    events.foreach{evt => {
      val dst = Event.editDistance(this,evt)
      //println(eventStrings.mkString(",") + "\t" + dst + "\t" + evt.eventStrings.mkString(","))
      if (dst <= 1) {
        if (probParent.isDefined)
          probParent = Some(probParent.get ++ Array[Event](evt))
        else
          probParent = Some(Array[Event](evt))

        found = true
      }
    }}
    return found
  }
}
//val evt1 = Event("test1","test1",1,Array[String]("4D-37","79D-65","79D-65","79D-65","79D-65","NONE","NONE","NONE","NONE","NONE"))
//val evt2 = Event("test2","test2",1,Array[String]("41D-20","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))

object Event {
  def editDistance(evt1: Event, evt2: Event): Int = {
    evt1.eventMap.map{case(event1,positions1) => {
      //println(event1 + "\t" + (evt2.eventMap contains event1))
      if (evt2.eventMap contains event1) {
        val positions2 = evt2.eventMap(event1)
        positions1.map{ps => if (positions2.indexOf(ps) < 0) 1 else 0}.sum
      } else {
        1
      }
    }}.sum
  }

  // calculate the basic distance between two points
  def probJointEvents(event1: String, event2: String,
    occuranceProportion: HashMap[String,Double],
    epsilon: Double = 0.001): Double = {

      // both not edited or both the same edit
      if ((event1 == "NONE" && event2 == "NONE") || (event1 != "NONE" && event2 == event1))
          1 - epsilon
      // event 1 is NONE, event2 is not
      else if (event1 == "NONE" && event2 != event1)
          1.0 - occuranceProportion(event2)
      else if (event2 == "NONE" && event2 != event1)
          1.0 - occuranceProportion(event1)
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
val eventNumberToEvents = HashMap[Int,HashMap[String,Array[Event]]]()

(0 until eventCount + 1).foreach{i => {
  positionToCounts(i) = new HashMap[String,Int]()
  positionToCountsNormalized(i) = new HashMap[String,Double]()
  eventNumberToEvents(i) = HashMap[String,Array[Event]]()
}}

// ------------------------------------------------------------------------------------
// load events from the stats file
// ------------------------------------------------------------------------------------
Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount * 2),line.split("\t").size - eventCount)
  val isAllNone = eventStrings.map{evt => if (evt == "NONE") 0 else 1}.sum

  if (isAllNone > 0) {
    val sp = line.split("\t")
    eventStrings.zipWithIndex.foreach{case(edit,index) => positionToCounts(index)(edit) = positionToCounts(index).getOrElse(edit,0) + 1}
    val id = eventStrings.mkString("-")

    // name: String, sample: String, numberOfReads: Int, eventStrings: Array[String], distanceScore: Int = 1
    val evt = Event(sp(1), sp(0), sp(4).toInt, eventStrings)
    totalRows += 1

    val editCount = evt.uniqueSize
    if (editCount > 10)
      println("line = " + evt.toOutputString)
    eventNumberToEvents(editCount)(id) = eventNumberToEvents(editCount).getOrElse(id,Array[Event]()) :+ evt
  }
}}
println("loaded " + totalRows + " events")

// ------------------------------------------------------------------------------------
// now normalize the positional counts by the total
// ------------------------------------------------------------------------------------
positionToCounts.foreach{case(position,counts) => {
  counts.foreach{case(event,count) => positionToCountsNormalized(position)(event) = count.toDouble / totalRows.toDouble}}
}

println("unique tags:")
(1 until 10).foreach{index => println( index + " " + eventNumberToEvents(index).size)}

val output = new PrintWriter("test.txt")

//41D-20,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE,NONE
//val evt1 = Event("test1","test1",1,Array[String]("4D-37","79D-65","79D-65","79D-65","79D-65","NONE","NONE","NONE","NONE","NONE"))
//val evt2 = Event("test2","test2",1,Array[String]("41D-20","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))
//println(Event.editDistance(evt1,evt2))
// name + "\t" + sample + "\t" + numberOfReads + "\t" + eventStrings.mkString(",") + "\t" + oth
output.write("name\tevent.count\tsample\treads\tevent\tparent.event\testimated.parent\tknown.parent\n")

println("reverse mapping:")
(2 until 10).foreach{index =>
  if (eventNumberToEvents contains index) {
    val befores = eventNumberToEvents(index - 1).map{case(id,hits) => hits(0)}.toArray
    val firsts = eventNumberToEvents(index).map{case(id,hits) => hits(0)}.toArray

    print("index " + index + " count " + firsts.length + " with prob = ")

    var tt = 0
    firsts.map{case(evt) => {
      if (evt.oneEdit(befores))
        tt += 1
      output.write(evt.toOutputString + "\n")
    }}
    println(tt)
  }
}
