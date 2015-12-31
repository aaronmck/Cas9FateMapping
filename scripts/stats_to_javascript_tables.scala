// generate the files needed for the d3 plots
import scala.io._
import java.io._
import scala.collection.mutable._

val statsFile = Source.fromFile(args(0)).getLines()
val header = statsFile.next()

val perBaseEvents = new PrintWriter(args(1))
val occurances    = new PrintWriter(args(2))
val readCounts    = new PrintWriter(args(3))
val allEventsF    = new PrintWriter(args(4))
val cutSites      = Source.fromFile(args(5)).getLines()

//val startPosition = args(5).toInt // 111
//val endPosition = args(6).toInt // 416

// load up the cutsite data and find the highest and lowest
// sites in the data, then
var cutheader = cutSites.next()
var startPosition = 1000
var endPosition = 0
cutSites.foreach{line => {
  val sp = line.split("\t")
  if (sp(1).toInt < startPosition)
    startPosition = sp(1).toInt
  if (sp(2).toInt > endPosition)
    endPosition = sp(2).toInt
}}

startPosition = math.max(0,startPosition - 20)
endPosition = endPosition + 20

println(startPosition + " , " + endPosition)


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

val posLength = endPosition - startPosition

val perBaseInsertions = new HashMap[Int,Int]()
val perBaseDeletions  = new HashMap[Int,Int]()
val perBaseMatches    = new HashMap[Int,Int]()

val eventToCount = new HashMap[String,Int]()
val eventToCigar = new HashMap[String,String]()

val plotTopXreads = 100

val matchTag = 0
val deletionTag = 1
val insertionTag = 2

/**
 * take event strings and split the data on the pluses
**/
def eventSplit(event:String): Array[Tuple3[Int,Int,Int]] = {
  var returnArray = Array[Tuple3[Int,Int,Int]]()

  var rawArray = Array[String](event)

  if (event contains "&") {
    rawArray = event.split("&")
  }

  rawArray.foreach{evt => {
    val evtsplit = evt.split("\\+")
    if (evtsplit.size > 1) {
      var position = evtsplit(1).toInt - startPosition
      var length = evtsplit(0).slice(0,evtsplit(0).length - 1).toInt

      val evtChar = evtsplit(0)(evtsplit(0).length - 1)
      if (evtChar == 'I') {
        //println(evtsplit.mkString("-"))
        returnArray :+= (position,position,insertionTag)
      } else if (evtChar == 'D')
        returnArray :+= (position,position+length,deletionTag)
      else
        println("FAILED on " + evtsplit(0)(evtsplit(0).length - 1))
    }
  }}
  returnArray
}

// for each line in the stats file, find the events, mark them on an array of events.
var totalReads = 0
statsFile.foreach{line => {
  if ((line contains "PASS") && !(line contains "WT") && !(line contains "UNKNOWN")) {
      totalReads += 1
      val sp = line.split("\t")
      val arrayPositions = new Array[Int]((endPosition - startPosition) + 1)

      var eventIntervals = Array[Tuple3[Int,Int,Int]]()

      val events = targetToIndex.foreach{case(target,index) =>
        val eps = eventSplit(sp(index))

        if (eps.size > 0) {
          eventIntervals = eventIntervals ++ eps
        }
      }

      //println(eventIntervals.mkString("--"))
      eventIntervals.foreach{interval =>
        //println(interval._1 + "\t" + interval._2 + "\t" + interval._3)
        (interval._1 until (interval._2 + 1)).foreach{indx => {
          if (indx >= 0 && indx < arrayPositions.size)
            arrayPositions(indx) = interval._3
          }
        }}

      arrayPositions.zipWithIndex.foreach{case(tag,ind) => {
        //println("TAG:" + tag + "\t" + ind)
        if (tag == matchTag)
          perBaseMatches(ind) = perBaseMatches.getOrElse(ind,0) + 1
        if (tag == deletionTag)
          perBaseDeletions(ind) = perBaseDeletions.getOrElse(ind,0) + 1
        if (tag == insertionTag)
          perBaseInsertions(ind) = perBaseInsertions.getOrElse(ind,0) + 1
      }}

      val eventString = arrayPositions.mkString("")
      //println(eventString)
      eventToCount(eventString) = eventToCount.getOrElse(eventString,0) + 1
      eventToCigar(eventString) = targetToIndex.map{case(target,index) => sp(index)}.mkString("_")
    }
}}


val allEvents = eventToCount.toSeq.sortBy(_._2).toArray
val topEvents = allEvents.slice(allEvents.size - plotTopXreads,allEvents.size).reverse
def toPct(count: Int): Double = (count.toDouble / totalReads.toDouble)

allEventsF.write("event\tarray\tcount\tproportion\n")
allEvents.zipWithIndex.foreach{case((str,count),index) => {
  allEventsF.write(eventToCigar(str) + "\t" + index + "\t" + count + "\t" + toPct(count) + "\n")
}}
allEventsF.close()


occurances.write("array\tposition\tevent\n")
readCounts.write("event\tarray\tcount\tWT\n")

topEvents.zipWithIndex.foreach{case((str,count),index) => {
  val isAllWT = eventToCigar(str).split("_").map{sat => if (sat == "NONE") 0 else 1}.sum == 0
  val outWT = if (isAllWT) 1 else 2
  readCounts.write(eventToCigar(str) + "\t" + index + "\t" + toPct(count) + "\t" + outWT + "\n")
  str.zipWithIndex.foreach{case(evt,subIndex) =>
    occurances.write(index + "\t" + subIndex + "\t" + evt + "\n")
  }
}}
occurances.close()
readCounts.close()

perBaseEvents.write("index\tmatch\tinsertion\tdeletion\n")
(0 until posLength).foreach{ind => {
  perBaseEvents.write(ind + "\t" + toPct(perBaseMatches.getOrElse(ind,0)) + "\t" +
  toPct(perBaseInsertions.getOrElse(ind,0)) + "\t" + toPct(perBaseDeletions.getOrElse(ind,0)) + "\n")
}}
perBaseEvents.close()
