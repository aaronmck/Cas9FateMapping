import scala.io._
import java.io._
import scala.collection.mutable._

// given either a stats or calls file, parse out the number of each event
// at each position, and show stats

val statsFile = Source.fromFile(args(0)).getLines()
val outputFile = new PrintWriter(args(1))

val header = statsFile.next().split("\t")

// --------------------------------------------------
// find the event columns -- make a mapping
// from the target number to the column index, assuming they're
// labeled targetX where is X is 1,2,3,4....
val targetToIndex = new HashMap[Int,Int]()
header.zipWithIndex.foreach{case(tk,index) => {
  if (tk contains "target") {
    val col = tk.stripPrefix("target").toInt
    targetToIndex(col) = index
    println("target" + col + " position " + index)
  }
}}


val matchTag = 0
val deletionTag = 1
val insertionTag = 2
val startPosition = 111
val endPosition = 416
val posLength = endPosition - startPosition

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
        returnArray :+= (position,length,insertionTag)
      } else if (evtChar == 'D')
        returnArray :+= (position,length,deletionTag)
      else
        println("FAILED on " + evtsplit(0)(evtsplit(0).length - 1))
    }
  }}
  returnArray
}

// our mapping variables to store event counts
val mapColToEventToCount = new HashMap[Int,HashMap[String,Int]]()
val mapEventToSiteCount = new HashMap[String,Int]()

// setup the hashmap for each site
targetToIndex.foreach{case(target,column) =>
  mapColToEventToCount(target) = new HashMap[String,Int]()
}

statsFile.foreach{line => {
  val sp = line.split("\t")
  val eventCounts = new HashMap[String,Int]()
  targetToIndex.foreach{case(target,index) => {
    mapColToEventToCount(target)(sp(index)) = mapColToEventToCount(target).getOrElse(sp(index),0) + 1
    eventCounts(sp(index)) = eventCounts.getOrElse(sp(index),0) + 1
  }}
  eventCounts.foreach{case(evt,cnt) => mapEventToSiteCount(evt) = cnt}
}}

outputFile.write("target\tevent\tcount\tposition\tlength\ttype\ttargetLength\n")

mapColToEventToCount.foreach{case(target,countMap) => {
  countMap.foreach{case(event,count) => {
    if (!(event contains "WT")) {
      val evt = eventSplit(event)
      evt.foreach{subEvt => {
        outputFile.write(target + "\t" + event + "\t" + count + "\t" + subEvt._1 + "\t" + subEvt._2 +
          "\t" + (if (subEvt._3 == deletionTag) "Deletion" else "Insertion") + "\t" + mapEventToSiteCount(event) + "\n")
      }}
    }}
  }
}}

outputFile.close()
