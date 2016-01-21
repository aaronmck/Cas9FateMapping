// --------------------------------------------------
// 4 gamete test script
//
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
val output2 = new PrintWriter(args(2))

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
//  for each line, for each target, add it to the columnToEventsArray
// -------------------------------------------------------------------
lines.foreach{line => {
  linesProcessed += 1
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
}}

output2.write("event\tcount\n")
var lMap = ListMap(eventToCount.toSeq.sortBy(_._2):_*)
lMap.foreach{case(event,count) => output2.write(event + "\t" + count + "\n")}
output2.close()

println("Lines processed " + linesProcessed)


output.write("event1\tevent2\t10\t01\t11\tfourGameteResult\n")

// --------------------------------------------------
// now for each event-pair across two columns, determine if two events
// satisfy the 4-gamete test
val run = new HashMap[String,Boolean]()

// loop over one column
val processedPairs = new HashMap[Int,Int]()
var totalCount = 0

columnToEventsArray.foreach{case(index,events) => {
  println(index + "\t" + events.result.size)

  // add combinations of events
  columnToEventsArray.filter(ind => ind != index).foreach{case(index2,events2) => {

    val counts = new HashMap[Int,HashMap[Int,Int]]()

    if (index > index2) {
      events.result.zip(events2.result).foreach{case(evt1Int,evt2Int) => {
        if (!(counts contains evt1Int))
          counts(evt1Int) = new HashMap[Int,Int]()
        counts(evt1Int)(evt2Int) = counts(evt1Int).getOrElse(evt2Int,0) + 1

        if (!(counts contains evt2Int))
          counts(evt2Int) = new HashMap[Int,Int]()
        counts(evt2Int)(evt1Int) = counts(evt2Int).getOrElse(evt1Int,0) + 1
      }}


      // now for each key, find out if there's a combination of both
      counts.foreach{case(e1Int,evt2Count) => {
        counts.foreach{case(e2Int,evt1Count) => {
          if (e1Int != 0 && e2Int != 0 && e1Int != e2Int && e1Int < e2Int) {
            val zeroOne = evt1Count.getOrElse(0,0)
            val oneZero = evt2Count.getOrElse(0,0)
            val oneOne  = evt1Count.getOrElse(e1Int,0)
            totalCount += 1
            val fourGamete = zeroOne > 0 && oneZero > 0 && oneOne > 0
            if (fourGamete) {
              output.write(numberToEvent(e1Int) + "\t" + numberToEvent(e2Int) +
                "\t" + oneZero + "\t" + zeroOne + "\t" + oneOne + "\t" + fourGamete + "\n")
            }
          }
        }}
      }}
    }
  }}
}}
println("total tried " + totalCount)
output.close()
