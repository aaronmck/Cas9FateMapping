// generate the files needed for the d3 plots
import scala.io._
import java.io._
import scala.collection.mutable._

// -----------------------------------------------------------------------------
// a container for HMIDs
// -----------------------------------------------------------------------------
sealed trait IndelType {
  val toInt = -1
  def toStr(): String
}

case object Deletion extends IndelType { override val toInt = 1; def toStr(): String = "D"}
case object Match extends IndelType { override val toInt = 0; def toStr(): String = "M"}
case object Insertion extends IndelType { override val toInt = 2; def toStr(): String = "I"}
case object NoneType extends IndelType { override val toInt = 3; def toStr(): String = "NONE"}

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
      if (adjPos >= 0  && adjPos + evt.size < eventInts.size) {
        (adjPos until (adjPos + evt.size)).foreach{pos => {
          evt.classOf match {
            case Deletion => eventInts(pos) = evt.classOf.toInt
            case Insertion if pos == adjPos => eventInts(pos) = evt.classOf.toInt
            case Insertion => eventInts(pos) = 0
            case _ => throw new IllegalStateException("Unable to match classOf")
          }
          eventInts(pos) = evt.classOf.toInt
        }}
      }
    }}
    return eventInts
  }
  def eventToPerBase2(startPosition: Int, endPosition: Int): Array[ETPB] = {
    var eventPBs = Array[ETPB]()
    events.foreach{evt => {
      if (evt.classOf != NoneType) {
        if (eventPBs.size == 0 && evt.position > startPosition)
          eventPBs :+= ETPB(startPosition, evt.position, 0)
        else if (eventPBs(eventPBs.size - 1).stop < evt.position)
          eventPBs :+= ETPB(eventPBs(eventPBs.size - 1).stop, evt.position, 0)

        if (evt.classOf == Deletion)
          eventPBs :+= ETPB(evt.position, evt.position + evt.size, evt.classOf.toInt)
        else
          eventPBs :+= ETPB(evt.position, evt.position + evt.size, evt.classOf.toInt)
      }
    }}
    if (eventPBs.size > 0 && eventPBs(eventPBs.size - 1).stop < endPosition) {
      eventPBs :+= ETPB(eventPBs(eventPBs.size - 1).stop, endPosition, 0)
    } else if (eventPBs.size == 0) {
      eventPBs :+= ETPB(startPosition, endPosition, 0)
    }
    return eventPBs
  }
}
case class ETPB(start: Int, stop: Int, event: Int)

// -----------------------------------------------------------------------------
// our main event class, which handles individual event entries
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
    if (b > a.stripPrefix("target").toInt) b else a.stripPrefix("target").toInt)
  println("Number of targets " + numberOfTargets)

  val targetStrings = (1 until (numberOfTargets + 1)).map{"target" + _}
  val targetToPosition = targetStrings.map{case(tg) => (tg,header.indexOf(tg))}.toMap
  val targetToNumber = targetStrings.map{case(tg) => (tg,tg.slice(tg.length-1,tg.length).toInt)}.toMap

  // process all lines in the file
  statsFile.foreach{line => {
    if ((line contains "PASS") && !(line contains "WT")) { 
      val (newHMIDString,newHMID) = lineToHMID(line)
      val replacementHMID = hmidCounts.getOrElse(newHMIDString,newHMID)
      replacementHMID.count += 1
      hmidCounts(newHMIDString) = replacementHMID
    }
  }}
  val totalHMIDs = hmidCounts.map{case(str,id) => id.count}.sum
  println("total " + totalHMIDs)
  val sortedEvents = hmidCounts.toSeq.sortBy(_._2.count).toArray.reverse
  println("total " + sortedEvents.map{case(evt,obj) => obj.count}.sum)

  /** process a line into an HMID **/
  def lineToHMID(line: String, unknownsToNone: Boolean = true): Tuple2[String,HMID] = {
    val spl = line.split("\t")

    val events = new ArrayBuffer[Event]()
    val tokens = new ArrayBuffer[String]()

    targetStrings.foreach{case(tg) => {
      val toAdd = if (spl(targetToPosition(tg)) == "UNKNOWN" && unknownsToNone) "NONE" else spl(targetToPosition(tg))
      tokens += toAdd
      spl(targetToPosition(tg)).split("\\&").foreach{
        subevt => {
          val converteSubEvt = if (subevt == "UNKNOWN" && unknownsToNone) "NONE" else subevt
          events += Event.toEvent(converteSubEvt,targetToNumber(tg)) // -- make sure we only propigate NONEs, not UNKNOWNs
        }
      }
    }}
    (tokens.mkString("_"),HMID(events.toArray))
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
val perBaseEvents = new PrintWriter(args(2))
val occurances    = new PrintWriter(args(1))
val readCounts    = new PrintWriter(args(3))
val allEventsF    = new PrintWriter(args(4))
val perBaseEventsNew = new PrintWriter(args(6))

// -----------------------------------------------------------------------------
// first output all of the events
// -----------------------------------------------------------------------------
allEventsF.write("event\tarray\tcount\tproportion\n")
statsObj.sortedEvents.zipWithIndex.foreach{case((hmid,hmidEvents),index) =>
  allEventsF.write(hmid + "\t" + index + "\t" + hmidEvents.count + "\t" + (hmidEvents.count.toDouble / statsObj.totalHMIDs.toDouble) + "\n")
}
allEventsF.close()

// -----------------------------------------------------------------------------
// now output the top events
// -----------------------------------------------------------------------------
readCounts.write("event\tarray\tproportion\trawCount\tWT\n")
statsObj.sortedEvents.slice(0,100).zipWithIndex.foreach{case((hmid,hmidEvents),index) =>
  readCounts.write(hmid + "\t" + index + "\t" + (hmidEvents.count.toDouble / statsObj.totalHMIDs.toDouble) + "\t" + hmidEvents.count + "\t" + (if(hmidEvents.isWT) 1 else 2) + "\n")
}
readCounts.close()

// -----------------------------------------------------------------------------
// output the top events as a melted string of 0s, 1s, and 2s (encoded indels)
// -----------------------------------------------------------------------------
val rangeBuffer = 20
val startPosition = cutSites.sites(0).start - rangeBuffer
val endPosition = cutSites.sites(cutSites.sites.size -1).cutsite + rangeBuffer

perBaseEvents.write("array\tposition\tevent\n")
statsObj.sortedEvents.slice(0,100).zipWithIndex.foreach{case((hmid,hmidEvents),index) => {
  hmidEvents.eventToPerBase(startPosition, endPosition).zipWithIndex.foreach{case(event,subIndex) =>
    perBaseEvents.write(index + "\t" + subIndex + "\t" + event + "\n")
  }
}}
perBaseEvents.close()

// -----------------------------------------------------------------------------
// create a simple table of the events
// -----------------------------------------------------------------------------

// ETPB(start: Int, stop: Int, event: Int)
println(startPosition + "\t" + endPosition)

perBaseEventsNew.write("array\tposition\tlength\tevent\n")
statsObj.sortedEvents.slice(0,100).zipWithIndex.foreach{case((hmid,hmidEvents),index) => {
  hmidEvents.eventToPerBase2(startPosition, endPosition).zipWithIndex.foreach{case(etpb,subIndex) =>
    perBaseEventsNew.write(index + "\t" + etpb.start + "\t" + etpb.stop + "\t" + etpb.event + "\n")
  }
}}
perBaseEventsNew.close()

// -----------------------------------------------------------------------------
// output per-base information
// -----------------------------------------------------------------------------
val insertionCounts = Array.fill[Int]((endPosition - startPosition) + 1)(0)
val deletionCounts = Array.fill[Int]((endPosition - startPosition) + 1)(0)

var totalReads = 0
statsObj.sortedEvents.foreach{case(hmid,hmidEvents) => {
  hmidEvents.eventToPerBase(startPosition, endPosition).zipWithIndex.foreach{case(event,index) => event match {
    case Deletion.toInt => deletionCounts(index) += hmidEvents.count
    case Insertion.toInt => insertionCounts(index) += hmidEvents.count
    case _ => {}
  }}
  totalReads += hmidEvents.count
}}

occurances.write("index\tmatch\tinsertion\tdeletion\n")
insertionCounts.zip(deletionCounts).zipWithIndex.foreach{case((ins,del),index) => {
  val matchProp = 1.0 - ((ins + del).toDouble / totalReads.toDouble)
  occurances.write(index + "\t" + matchProp + "\t" + (ins.toDouble/totalReads.toDouble) + "\t" + (del.toDouble/totalReads.toDouble) + "\n")
}}
occurances.close()
