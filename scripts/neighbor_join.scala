import scala.io._
import java.io._
import scala.collection.mutable._
import scala.util._

// Aaron -- November 12th, 2015
//
// perform neighbor-joining with a constrained matrix - some combinations
// are not allowed, we choose only allowed combinations.  The exercise to
// determine the effect on distance estimates is left to the reader....

// our input requirements
val inputUMIfile = args(0)
val outputTree = args(1)
val outputAnnotations = args(2)
val countsOutput = args(3)

// read in the series of events
var events = Array[Event]()
var allEvents = Array[Event]()
val uniqueEvents = new HashMap[String, String]()
val targetSiteCount = 10



// ------------------------------------------------------------------------------------------------------------------------
//
// Helper class -- contains most things we know about one of our events
//
// ------------------------------------------------------------------------------------------------------------------------
case class Event(name: String, sample: String, numberOfReads: Int, eventStrings: Array[String], inferred: Boolean = false) {

  var left: Option[Event] = None
  var right: Option[Event] = None
  var branchLeft = 0.0
  var branchRight = 0.0

  // are two events compatible? they're compatible if for all sites they either:
  // 1) share the same edits
  // 2) either or both have a NONE
  // they must satisfy these criteria at every position to be compatible
  def compatible(event2: Event): Boolean = {
    // now check our compatibility
    var sm = 0
    sm += Event.scoreInternal(this, event2)
    sm += Event.scoreInternal(event2, this)

    // our compatible condition and do the reverse in case event2 has merged nodes
    sm == 0
  }

  // merge down the events to a common event on the merging event
  def merge(event2: Event, branchL: Double, branchR: Double, mergeID: Int): Event = {
    val newEvtStr = eventStrings.zip(event2.eventStrings).map{case(evt1,evt2) => {
      if (evt1 != "NONE" && evt2 != "NONE" && evt2 != evt1)
        throw new IllegalStateException("UNABLE TO MERGE " + toFancyString + " and " + event2.toFancyString)

      if ((evt1 == "NONE" && evt2 != "NONE") || (evt2 == "NONE" && evt1 != "NONE") ) {
        "NONE"
      } else {
        if (evt1 == evt2)
          evt1
        else
          "NONE"
      }
    }}.toArray

    val ret = Event("MERGED" + mergeID, // "M(" + name + "),(" + event2.name + ")",
                    "M_L" + sample + "_" + event2.sample + "J",
                    numberOfReads + event2.numberOfReads,
                    newEvtStr,
                    true)
    ret.left = Some(this)
    ret.right = Some(event2)
    ret.branchLeft = branchL
    ret.branchRight = branchR

    return ret
  }
  def toFancyString(): String = name + "\t" + sample + "\t" + numberOfReads + "\t" + eventStrings.mkString("\t") + "\t" + left.isDefined + "\t" + right.isDefined

  // newick is so weird...
  def newickString(distance: Double, annotationOutput: PrintWriter): String = {
    val leftStr =  if (left.isDefined)  left.get.newickString(branchLeft, annotationOutput)   else ""
    val rightStr = if (right.isDefined) right.get.newickString(branchRight, annotationOutput) else ""

    var newDist = if (distance.isNaN) 1.0 else distance
    val ret = name + numberOfReads + ":" + newDist
    // taxa        attrib1        attrib2
    var newSample = if (sample startsWith "M") "_" else sample
    annotationOutput.write(name + numberOfReads + "\t" + name + "\t" + newSample + "\t" + numberOfReads + "\t" + eventStrings.mkString("-") + "\n")

    //println(ret + "(" + leftStr + "," + rightStr + ")")
    val rnd = new Random()

    if (leftStr != "")
      if (rnd.nextDouble < 0.5)
        "(" + leftStr + "," + rightStr + ")" + ret
      else
        "(" + rightStr + "," + leftStr + ")" + ret
    else
      ret
  }
}

object Event {
  def scoreInternal(subNode: Event, event2: Event): Int = {
    var sm = subNode.eventStrings.zip(event2.eventStrings).map{case(evt1,evt2) =>
      if (evt1 == "NONE" || evt2 == "NONE" || evt2 == evt1)
        0
      else
        1
    }.sum
    /*if (subNode.left.isDefined)
      sm += Event.scoreInternal(subNode.left.get,event2)
    if (subNode.right.isDefined)
      sm += Event.scoreInternal(subNode.right.get,event2)*/
    return sm
  }
}


// now create a distance matrix of all of the unique samples:
// each inner array represents a single sample for ease of
// access / modification

def calculateDistances(events: Array[Event], eventToDistance: HashMap[String,Double], noneScore: Double = 0): Array[Array[Double]] = {
  var distances = Array[Array[Double]]()

  events.foreach{event => {
    val innnerDistance = new Array[Double](events.size)
    events.zipWithIndex.foreach{case(event2,index) => {
        innnerDistance(index) = event.eventStrings.zip(event2.eventStrings).map{case(evt1,evt2) => {
          if (evt1 == "NONE" && evt2 == "NONE")
            noneScore
          else if (evt1 == "NONE" && evt2 != "NONE")
            eventToDistance(evt2)
          else if (evt1 != "NONE" && evt2 == "NONE")
            eventToDistance(evt1)
          else if (evt1 != evt2)
            eventToDistance(evt1) + eventToDistance(evt2)
          else
            0
        }}.sum
    }}
    distances :+= innnerDistance
  }}
  return distances
}


// ------------------------------------------------------------------------------------------------------------------------
// parse out the event strings into events
// ------------------------------------------------------------------------------------------------------------------------
println("reading in event file...")
Source.fromFile(inputUMIfile).getLines().drop(1).foreach { line => {
  val eventStrings = line.split("\t").slice(line.split("\t").size - (targetSiteCount * 2), line.split("\t").size - targetSiteCount)
  val sp = line.split("\t")
  val eventAsString = eventStrings.mkString("_")
  val evt = Event(sp(1), sp(0), sp(4).toInt, eventStrings)
  allEvents :+= evt
  if (!(uniqueEvents contains eventAsString)) {
    val isAllNone = eventStrings.map { evt => if (evt == "NONE") 0 else 1 }.sum
    if (isAllNone > 0) {
      events :+= evt
    }
  }

  uniqueEvents(eventAsString) = "Sample_" + line.split("\t")(0)

}}

// ------------------------------------------------------------------------------------------------------------------------
// figure out a meaningful distance metric for our sequences -- if the events are
// share a really rare event, their distance is really really small.  If the event
// is common they have a medium to higher distance
// ------------------------------------------------------------------------------------------------------------------------

val totalEvents = (allEvents.size * targetSiteCount).toDouble
val eventToCount = new HashMap[String,Int]()
allEvents.foreach{evt1 => evt1.eventStrings.foreach{evtString => {
  eventToCount(evtString) = eventToCount.getOrElse(evtString,0) + 1
}}}

val eventToDistance = new HashMap[String,Double]()

// for now a simple exercise is to use the counts data -- maybe in the future we want to
// try using a the log2(counts) which would be more accurate or counts^2 which would add
// a lot more distance to common events
val countFile = new PrintWriter(countsOutput)
var minCount = Double.MaxValue
eventToCount.foreach{case(evtStr,count) =>
  if (count < minCount)
    minCount = count
}

eventToCount.foreach{case(evtStr,count) => {
  eventToDistance(evtStr) = math.log(count.toDouble)/math.log(2) + 1 // math.log(count.toDouble) / math.log(2) + 1
  countFile.write(evtStr + "\t" + (math.log(count.toDouble) / math.log(2) + 1) + "\t" + count + "\t" + totalEvents + "\n") // totalEvents
}}
//countFile.close()



var distances = calculateDistances(events,eventToDistance)
//distances.foreach{col => println(col.mkString(","))}

// keep track of mapping of event to it's distance matrix position
var indexToEvent = new HashMap[Int,Event]()
events.zipWithIndex.foreach{case(event,index) => indexToEvent(index) = event}
println("Starting nodes: " + indexToEvent.size)


println("calculating global site compatibility....")
// the compatibility test is taking a loooot of time, cache the results, so we only
// have to update it for new pairs
val compatibility = new HashMap[Event,HashMap[Event,Boolean]]()
events.foreach{event1 => {
  events.foreach{event2 => {
    if (!(compatibility contains event1))
      compatibility(event1) = new HashMap[Event,Boolean]()
    compatibility(event1)(event2) = event1.compatible(event2)
  }}
}}

// cache the minimum for each column
var columnToMinPositionTemp = ArrayBuilder.make[Int]
distances.zipWithIndex.foreach{case(column, index) => {
  var minValue = Double.MaxValue
  var minPosition = -1
  for (i <- 0 until column.size) {
    if (i != index && column(i) < minValue && compatibility(indexToEvent(index))(indexToEvent(i))) {
      minValue = column(i)
      minPosition = i
    }
  }
  columnToMinPositionTemp += minPosition
}}
var columnToMinPosition = columnToMinPositionTemp.result


var mergeID = 0


// now the main session: take the distance matrix and perform constrained neighbor joining
var still_nodes_compatible = true
val debug = false

println("Performing merges (dot per merge)")
while (still_nodes_compatible) {
  print(".") // + distances.size)
  // find a list of nodes that satisfy the neightbor joining criteria

  // find the average distances
  val mu = (0 until distances.size).map{arrayInd => {
    distances(arrayInd).zipWithIndex.map{case(datum,index) =>
      if (arrayInd == index)
        0
      else
        datum / (distances.size.toDouble - 2)
    }.sum
  }}
  if (debug)
    print("*")
  // now find the i and j for which min(Dij - Ui - Uj) <AND> the merge is compatible
  var mini = -1
  var minj = -1

  // make the minimum finding much faster -- we've cached the location of the minimum in every column (that's compatible), find the mini,minj
  // by finding the min of the columns.
  var minCol = -1
  var minValue = Double.MaxValue

  if (debug) {
    println("columnToMinPosition.size " + columnToMinPosition.size)
    println("distances.size " + distances.size)
  }
  for (i <- 0 until columnToMinPosition.size) {
      if (columnToMinPosition(i) >= 0 && columnToMinPosition(i) >= 0 && distances(i)(columnToMinPosition(i)) < minValue) {
        minValue = distances(i)(columnToMinPosition(i))
        minCol = i
      }
  }

  println("mini " + mini + " minj " + minj)
  mini = minCol
  if (mini >= 0)
    minj = columnToMinPosition(mini)


  if (debug)
    print("$")
  // did we find a pair to merge?
  if (mini >= 0 && minj >= 0) {
    val branchDisti = (0.5 * distances(mini)(minj)) + (0.5 * (mu(mini) - mu(minj)))
    val branchDistj = (0.5 * distances(mini)(minj)) + (0.5 * (mu(minj) - mu(mini)))

    if (debug) {
      println(indexToEvent(mini).sample + "\t" + indexToEvent(mini).eventStrings.mkString("\t"))
      println(indexToEvent(minj).sample + "\t" + indexToEvent(minj).eventStrings.mkString("\t"))
      println(distances(mini)(minj))
      println(mu(mini))
      println(mu(minj))
      println(0.5 * (mu(mini) - mu(minj)))
      println(0.5 * (mu(minj) - mu(mini)))
      println()
    }
    // make mini < minj
    if (mini > minj) {
      val tmp = mini
      mini = minj
      minj = tmp
    }
    if (mini == minj) {
      throw new IllegalStateException("Crap " + mini + " is equal to " + minj)
    }

    countFile.write("MERGING\t" + indexToEvent(mini).sample + "\t" + indexToEvent(minj).sample + "\t" + indexToEvent(mini).eventStrings.mkString("\t") + "\t" + indexToEvent(minj).eventStrings.mkString("\t") + "\n")

    // merge the two events
    val merged = indexToEvent(mini).merge(indexToEvent(minj), branchDisti, branchDistj, mergeID)
    mergeID += 1

    // rescale the compatibility matrix
    compatibility(merged) = new HashMap[Event,Boolean]()
    indexToEvent.foreach{case(index,event) => {
      compatibility(event)(merged) = event.compatible(merged)
      compatibility(merged)(event) = merged.compatible(event)
    }}

    if (debug)
      print("!")
    // recompute the distance for our new node and put this into the tree
    var newDistance = Array[Array[Double]]()
    var columnToMinPositionTempOld = ArrayBuilder.make[Int]
    var savedLastColumn = ArrayBuilder.make[Double] // shortcut -- we save the last row which equals the last column plus a zero on the end
    if (debug)
      println("mini " + mini + " minj " + minj + " " + minValue)
    var newIndex = 0
    for (a <- 0 until distances.size) {
      if (a != mini && a != minj) {
        var newDistanceInner = ArrayBuilder.make[Double]
        newDistanceInner ++= distances(a).slice(0,mini)
        newDistanceInner ++= distances(a).slice(mini+1,minj)
        newDistanceInner ++= distances(a).slice(minj+1,distances.length)
        val newNodeDist = ((distances(a)(mini) + distances(a)(mini) - distances(mini)(minj)) / 2.0)
        newDistanceInner += newNodeDist
        val newDist = newDistanceInner.result

        newDistance :+= newDist
        savedLastColumn += newNodeDist

        val oldValue = columnToMinPosition(a)

        var step = ""
        // fix the min distance counts -- if we used to have mini or minj as our
        // min (an unlikely tie)
        //println("before " + columnToMinPosition(a))
        //if (columnToMinPosition(a) == mini || columnToMinPosition(a) == minj) {
          step = "REFIND"
          var minV = Double.MaxValue
          var minP = -1
          for (i <- 0 until newDist.size) {
            if (i != newIndex && newDist(i) < minValue && compatibility(indexToEvent(a))(indexToEvent(i))) {
              minV = newDist(i)
              minP = i
            }
          }
          columnToMinPosition(a) = minP
        /*} else {
          columnToMinPosition(a) = columnToMinPosition(a) - (a - newIndex)
        }

        if (compatibility(merged)(indexToEvent(a)) && columnToMinPosition(a) >= 0 && distances(a)(columnToMinPosition(a)) > newNodeDist) {
          step = "final"
          columnToMinPosition(a) = newDist.length - 1
        }
        //println("after " + columnToMinPosition(a))
        */
        columnToMinPositionTempOld += columnToMinPosition(a)

        if (columnToMinPosition(a) == columnToMinPositionTempOld.result.size) {
          print("OLD value = " + oldValue + " new value " + columnToMinPosition(a) + " column " + columnToMinPositionTempOld.result.size)
          println(" vect " + newDist.size + " step " + step + " mini " + mini + " minj " + minj)
        }
        newIndex += 1
      }
    }

    savedLastColumn += 0.0// last column becomes the last row plus zero (self distance)
    val slcResult = savedLastColumn.result
    newDistance :+=slcResult

    var minVal = Double.MaxValue
    var minPos = -1
    for (i <- 0 until (slcResult.size - 1)) {
      if (slcResult(i) < minValue && compatibility(merged)(indexToEvent(i))) {
        minVal = slcResult(i)
        minPos = i
      }
    }

    columnToMinPositionTempOld += minPos
    columnToMinPosition = columnToMinPositionTempOld.result
    if (debug) {
      println()
      println("coltominpos size = " + columnToMinPosition.size)
      println("Size = " + columnToMinPosition.size)
    }

    if (debug)
      print("_")

    distances = newDistance
    if (debug)
      println("distances " + distances.size)
    var newIndexToEvent = new HashMap[Int,Event]()
    var cnt = 0
    for (i <- 0 until indexToEvent.size) {
      if (i != mini && i != minj) {
        newIndexToEvent(cnt) = indexToEvent(i)
        cnt += 1
      }
    }
    newIndexToEvent(cnt) = merged
    if (debug)
      print("_")

    indexToEvent = newIndexToEvent
    if (debug)
      print("@")
  } else {
    still_nodes_compatible = false
  }

}


indexToEvent.foreach{case(index,value) => countFile.write("FINAL\t" + value.toFancyString + "\n")}

countFile.close()

val annotationOutput = new PrintWriter(outputAnnotations)
annotationOutput.write("taxa\tname\tsample\tnumberOfReads\teventStrings\n")

val outputNewick = new PrintWriter(outputTree)
indexToEvent.foreach{case(index,value) => outputNewick.write("(" + value.newickString(1.0,annotationOutput) + ")root;\n")}
outputNewick.close()
annotationOutput.close()

println("Remaining nodes: " + indexToEvent.size)
