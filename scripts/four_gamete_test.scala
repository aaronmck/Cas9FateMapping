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


val output = new PrintWriter(args(1))

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

val columnToEventsArray = new HashMap[Int,ArrayBuilder[Int]]()
targetToIndex.foreach{case(target,column) => columnToEventsArray(target) = ArrayBuilder.make[Int]}

var linesProcessed = 0

lines.foreach{line => {
  linesProcessed += 1
  val sp = line.split("\t")
  targetToIndex.foreach{case(target,index) => {
    val event = sp(index)

    if (event != "NONE" && event != "UNKNOWN") {
      if (!(eventToNumber contains event)) {
        eventToNumber(event) = nextIndex
        numberToEvent(nextIndex) = event
        nextIndex += 1
      }
      columnToEventsArray(target) = columnToEventsArray.getOrElse(target,ArrayBuilder.make[Int]) += eventToNumber(event)
    }
  }}
}}

println("Lines processed " + linesProcessed)

case class FourGameteResult(evt1: String, evt2: String, zeroOne: Int = 0, oneZero: Int = 0, oneOne: Int = 0) {

  def addZeroZero(): FourGameteResult = FourGameteResult(evt1, evt2, zeroOne, oneZero, oneOne)
  def addZeroOne():  FourGameteResult = FourGameteResult(evt1, evt2, zeroOne + 1, oneZero, oneOne)
  def addOneZero():  FourGameteResult = FourGameteResult(evt1, evt2, zeroOne, oneZero + 1, oneOne)
  def addOneOne():   FourGameteResult = FourGameteResult(evt1, evt2, zeroOne, oneZero, oneOne + 1)

  def addEvent(event1: Int, event2: Int) : FourGameteResult = {
    if ((event1 == 0) && (event2 == 0))
      return addZeroZero()
    if ((event1 == 0) && !(event2 == 0))
      return addZeroOne()
    if (!(event1 == 0) && (event2 == 0))
      return addOneZero()
    if (!(event1 == 0) && !(event2 == 0))
      return addOneOne()
    println(event1 + "\t" + event2)
    return addZeroZero()
  }

  def zeroProp(total:Int): Double = (zeroOne.toDouble + oneOne.toDouble) / total.toDouble
  def oneProp(total:Int): Double = (oneZero.toDouble + oneOne.toDouble) / total.toDouble
  def gametes(total:Int): Double = (zeroOne.toDouble + oneZero.toDouble + oneOne.toDouble) / total.toDouble

}


val minCount = 20

output.write("column1\tcolumn2\tevent1\tevent2\tevent1Prop\tevent2Prop\t01\t10\t11\tgameteCount\n")

// --------------------------------------------------
// now for each event-pair across two columns, determine if two events
// satisfy the 4-gamete test
// output.write("column1\tcolumn2\tevent1\tevent2\tevent1Prop\tevent2Prop\t00\t01\t10\t11\tgameteCount\n")
val run = new HashMap[String,Boolean]()

// loop over one column
columnToEventsArray.foreach{case(index,events) => {
  val counts = new HashMap[Int,Int]()
  events.result.foreach{case(evt) => counts(evt) = counts.getOrElse(evt,0) + 1}

  // loop another column, checking that we haven
  columnToEventsArray.foreach{case(index2,events2) => {
    if (index != index2 && !(run contains (index + "," + index2)) && !(run contains (index2 + "," + index))) {
      println(index + "," + index2)

      val fourEvents = new HashMap[String,FourGameteResult]()
      events.result.zip(events2.result).foreach{case(evt1,evt2) => {
        var tag = evt1 + "," + evt2
        if (evt1 > evt2)
          tag = evt1 + "," + evt2

        fourEvents(tag) = fourEvents.getOrElse(tag,FourGameteResult(numberToEvent(evt1),numberToEvent(evt2))).addEvent(evt1,evt2)
      }}
      fourEvents.foreach{case(tg,fgr) => output.write(index + "\t" + index2 + "\t" + fgr.evt1 + "\t" + fgr.evt2 + "\t" + fgr.zeroProp(events2.result.length) + "\t" + fgr.oneProp(events2.result.length) +
        "\t" + fgr.zeroOne + "\t" + fgr.oneZero + "\t" + fgr.oneOne + "\t" + fgr.gametes(events2.result.length) + "\n")}

      run(index + "," + index2) = true
    }
  }}
}}

output.close()
