// --------------------------------------------------
// 4 gamete test script
//
// --------------------------------------------------
import scala.io._
import java.io._
import scala.collection.mutable._


// --------------------------------------------------
// load up the lines from our stats/calls file
val lines = Source.fromFile(args(0)).getLines()
val header = lines.next()

// --------------------------------------------------
// find the event columns -- make a mapping
// from the target number to the column index, assuming they're
// labeled targetX where is X is 1,2,3,4....
val targetToIndex = new LinkedHashMap[Int,Int]()
header.split("\t").zipWithIndex.foreach{case(tk,index) => {
  if (tk contains "target") {
    val col = tk.stripPrefix("target").toInt
    targetToIndex(col) = index
    println("target" + col + " position " + index)
  }
}}

// load up the event columns into arrays of numeric events
// --------------------------------------------------
val eventToNumber = new HashMap[String,Int]()
val numberToEvent = new HashMap[Int,String]()

eventToNumber("NONE") = 0
eventToNumber("UNKNOWN") = 0
var nextIndex = 1

val eventToEventCount   = new HashMap[String,HashMap[String,Int]]()
val eventCounts         = new HashMap[String,Int]()

var linesProcessed = 0

lines.foreach{line => {
  linesProcessed += 1
  val sp = line.split("\t")
  targetToIndex.foreach{case(target1,index1) => {
    targetToIndex.foreach{case(target2,index2) => {
      if (index1 != index2) {
        if (!(eventToEventCount contains sp(index1))) {
          //println(sp(index1))
          eventToEventCount(sp(index1)) = new HashMap[String,Int]()
        }
        eventToEventCount(sp(index1))(sp(index2)) = eventToEventCount(sp(index1)).getOrElse(sp(index2),0) + 1

        eventCounts(sp(index1)) = eventCounts.getOrElse(sp(index1),0) + 1
        eventCounts(sp(index2)) = eventCounts.getOrElse(sp(index2),0) + 1
      }
    }}
  }}
  if (linesProcessed % 10000 == 0)
    println("Lines processed " + linesProcessed)
}}

println("Lines processed " + linesProcessed)
val tSize = eventToEventCount.size.toDouble
println(tSize)
val output = new PrintWriter(args(1))

output.write("event1\tevent2\tevent1Prop\tevent2Prop\t01\t10\t11\n")
eventToEventCount.foreach{case(evt1,counts) => {
    counts.foreach{case(evt2,count) => {
      //println(eventCounts contains evt1)
      //println(eventCounts contains evt2)
      val evt1Prop = eventCounts(evt1).toDouble / (linesProcessed * targetToIndex.size)
      val evt2Prop = eventCounts(evt2).toDouble / (linesProcessed * targetToIndex.size)
      val oneToTwo = eventToEventCount(evt1)(evt2)
      val twoToOne = eventToEventCount(evt2)(evt1)

      output.write(evt1 + "\t" + evt2 + "\t" + evt1Prop + "\t" + evt2Prop + "\t" +
        (eventCounts(evt1) - oneToTwo) + "\t" + (eventCounts(evt2) - twoToOne) +
        "\t" + eventToEventCount(evt1)(evt2) + "\n")
    }}
}}
output.close()
