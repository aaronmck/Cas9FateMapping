// generate the files needed for the d3 plots
import scala.io._
import java.io._
import scala.collection.mutable._

// -----------------------------------------------------------------------------
// a container for HMIDs
// -----------------------------------------------------------------------------
sealed trait IndelType {
  def toInt(): Int
  def toStr(): String
}

case object Deletion extends IndelType { def toInt(): Int = 1; def toStr(): String = "D"}
case object Match extends IndelType { def toInt(): Int = 0; def toStr(): String = "M"}
case object Insertion extends IndelType { def toInt(): Int = 2; def toStr(): String = "I"}
case object NoneType extends IndelType { def toInt(): Int = 3; def toStr(): String = "NONE"}

// various data storage classes we have
case class Cutsite(sequence: String, start: Int, cutsite:Int)
case class Reference(name: String, sequence: String, primer1: String, primer2:String)

// HMID events
case class HMID(events: Array[Event]) {
  val stringRep = events.map{evt => evt.toStringRep()}.mkString("_")
  var count = 0
  def isWT(): Boolean = events.map{evt => if (evt.classOf == NoneType) 0 else 1}.sum == 0

  def eventToPerBase(startPosition: Int, endPosition: Int): Array[Int] = {
    val eventInts = Array.fill[Int]((endPosition - startPosition) + 1)(0)
    events.foreach{evt => {
      val adjPos = evt.position - startPosition
      if (adjPos < 0)
        println("adj position " + adjPos + " " + evt.toStringRep + " " + startPosition)
        (adjPos until (adjPos + evt.size)).foreach{pos => {
          evt.classOf match {
            case Deletion => eventInts(pos) = evt.classOf.toInt
            case Insertion if pos == adjPos => eventInts(pos) = evt.classOf.toInt
            case Insertion => eventInts(pos) = 0
            case _ => throw new IllegalStateException("Unable to match classOf")
          }
          eventInts(pos) = evt.classOf.toInt
        }}
    }}
    return eventInts
  }
}

// -----------------------------------------------------------------------------
// our main event class, which handles individual entries
// -----------------------------------------------------------------------------
case class Event(site: Int, size: Int, classOf: IndelType, bases: Option[String], position: Int) {
  def toStringRep(): String = if (classOf == NoneType) "NONE" else size + classOf.toStr() + "+" + position + (if (bases.isDefined) ("+" + bases.get) else "")
}

object Event {
  def toEvent(substr: String, site: Int): Event = {
    val tokens = substr.split("\\+")
    if (tokens.length == 1)
      return Event(site,0,NoneType,None,-1)

    val size = tokens(0).slice(0,tokens(0).length -1).toInt
    val typeOf = tokens(0).slice(tokens(0).length-1,tokens(0).length) match {
      case "D" => Deletion
      case "I" => Insertion
      case "M" => Match
      case _ => throw new IllegalStateException("Unknown type >" + tokens(0).slice(tokens(0).length-1,tokens(0).length) + "<")
    }
    Event(site,size, typeOf,if(tokens.length == 3) Some(tokens(2)) else None, tokens(1).toInt)
  }
}

// -----------------------------------------------------------------------------
// a container for HMIDs
// -----------------------------------------------------------------------------
class StatsFile(inputFile: String) {
  val statsFile = Source.fromFile(inputFile).getLines()
  val hmidCounts = new HashMap[String,HMID]()

  // setup a bunch of of ways to index the target names
  val header = statsFile.next().split("\t")
  val numberOfTargets = header.filter{tk => tk contains "target"}.foldLeft(0)((b,a) =>
    if (b > a.slice(a.length-1,a.length).toInt) b else a.slice(a.length-1,a.length).toInt)
  val targetStrings = (1 until (numberOfTargets + 1)).map{"target" + _}
  val targetToPosition = targetStrings.map{case(tg) => (tg,header.indexOf(tg))}.toMap
  val targetToNumber = targetStrings.map{case(tg) => (tg,tg.slice(tg.length-1,tg.length).toInt)}.toMap

  var totalHMIDs = 0
  // process all lines in the file
  statsFile.foreach{line => {
    if ((line contains "PASS") && !(line contains "WT") && !(line contains "UNKNOWN")) {
      val newHMID = lineToHMID(line)
      val replacementHMID = hmidCounts.getOrElse(newHMID.stringRep,newHMID)
      replacementHMID.count += 1
      hmidCounts(newHMID.stringRep) = replacementHMID
      totalHMIDs += 1
    }
  }}
  val aboveThresh = hmidCounts.filter{case(str,id) => id.count >= 10}.size
  println("Processed " + hmidCounts.size + " unique HMIDs, with " + aboveThresh + " having 10 of more occurances, from a total of " + totalHMIDs + " HMIDs in the file")
  val sortedEvents = hmidCounts.toSeq.sortBy(_._2.count).toArray.reverse
  (0 until 30).foreach{ind => println(sortedEvents(ind)._1 + "\t" + sortedEvents(ind)._2.count)}

  /** process a line into an HMID **/
  def lineToHMID(line: String): HMID = {
    val spl = line.split("\t")

    val events = new ArrayBuffer[Event]()
    targetStrings.foreach{case(tg) => spl(targetToPosition(tg)).split("\\&").foreach{subevt => events += Event.toEvent(subevt,targetToNumber(tg))}}
    HMID(events.toArray)
  }
}
// -----------------------------------------------------------------------------
// store the cutsites
// -----------------------------------------------------------------------------
class CutSiteContainer(cutSiteFile:String) {
  val csFile = Source.fromFile(cutSiteFile).getLines()
  val header = csFile.next()

  val sites = csFile.map{case(line) => {
    val sp = line.split("\t")
    println(line)
    Cutsite(sp(0), sp(1).toInt, sp(2).toInt)
  }}.toArray
}

// -----------------------------------------------------------------------------
// process the input files
// -----------------------------------------------------------------------------
val statsObj = new StatsFile(args(0))
val cutSites = new CutSiteContainer(args(5))

// -----------------------------------------------------------------------------
// now create output files
// -----------------------------------------------------------------------------
val perBaseEvents = new PrintWriter(args(1))
val occurances    = new PrintWriter(args(2))
val readCounts    = new PrintWriter(args(3))
val allEventsF    = new PrintWriter(args(4))

// -----------------------------------------------------------------------------
// first output all of the events
// -----------------------------------------------------------------------------
allEventsF.write("event\tarray\tcount\tproportion\n")
statsObj.sortedEvents.zipWithIndex.foreach{case((hmid,hmidEvents),index) =>
  allEventsF.write(hmid + "\t" + index + "\t" + hmidEvents.count + "\t" + (hmidEvents.count.toDouble / statsObj.sortedEvents.size.toDouble) + "\n")
}
allEventsF.close()

// -----------------------------------------------------------------------------
// now output the top events
// -----------------------------------------------------------------------------
readCounts.write("event\tarray\tproportion\trawCount\tWT\n")
statsObj.sortedEvents.slice(0,100).zipWithIndex.foreach{case((hmid,hmidEvents),index) =>
  readCounts.write(hmid + "\t" + index + "\t" + hmidEvents.count + "\t" + (hmidEvents.count.toDouble / statsObj.sortedEvents.size.toDouble) + "\t" + hmidEvents.isWT + "\n")
}
readCounts.close()

// -----------------------------------------------------------------------------
// output the top events as a melted string of 0s, 1s, and 2s (encoded indels)
// -----------------------------------------------------------------------------
val rangeBuffer = 20
val startPosition = cutSites.sites(0).start - rangeBuffer
val endPosition = cutSites.sites(cutSites.sites.size -1).cutsite + rangeBuffer

perBaseEvents.write("array\tposition\tevent\n")
statsObj.sortedEvents.slice(0,100).zipWithIndex.foreach{case((hmid,hmidEvents),index) =>
  hmidEvents.eventToPerBase(startPosition, endPosition).zipWithIndex.foreach{case(event,subIndex) =>
    perBaseEvents.write(index + "\t" + subIndex + "\t" + event + "\n")
  }
}
perBaseEvents.close()

// -----------------------------------------------------------------------------
// output per-base information
// -----------------------------------------------------------------------------
val insertionCounts = Array.fill[Int]((endPosition - startPosition) + 1)(0)
val deletionCounts = Array.fill[Int]((endPosition - startPosition) + 1)(0)

statsObj.sortedEvents.foreach{case(hmid,hmidEvents) =>
  hmidEvents.eventToPerBase(startPosition, endPosition).zipWithIndex.foreach{case(event,index) => event match {
    case 1 => deletionCounts(index) += 1
    case 2 => insertionCounts(index) += 1
    case _ => {}
  }}
}

occurances.write("index\tmatch\tinsertion\tdeletion\n")
val totalEvents = statsObj.sortedEvents.size
insertionCounts.zip(deletionCounts).zipWithIndex.foreach{case((ins,del),index) => {
  val matchProp = 1.0 - ((ins + del).toDouble / totalEvents.toDouble)
  occurances.write(index + "\t" + matchProp + "\t" + (ins.toDouble/totalEvents.toDouble) + "\t" + (del.toDouble/totalEvents.toDouble) + "\n")
}}
occurances.close()
