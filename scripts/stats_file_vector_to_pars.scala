// load up a stats file and produce a distance matrix using a simple scoring scheme
import scala.io._
import java.io._
import scala.collection.mutable._
import scala.math._

val output = new PrintWriter(args(1))
val eventCount = args(2).toInt
val annotations = new PrintWriter(args(3))
annotations.write("taxa\tembryo\tumi\tevents\n")

val mapping = new HashMap[String,String]()

// ------------------------------------------------------------------------------------------
// our output functions
// ------------------------------------------------------------------------------------------
def toOutputString(bitSizes: Array[Int], eventStrings: Array[String]) = eventStrings.zipWithIndex.map{case(evt,index) => toBitString(evt,bitSizes(index))}.mkString("")
def toBitString(evt: String, length: Int): String = {
  val ret = new Array[Int](length)
  ret(evt.toInt) = 1
  return(ret.mkString(""))
}

val maxByPos = new Array[Int](eventCount)

// ------------------------------------------------------------------------------------------
// for each line in the input file, pull out the events and find the count so that we
// can make a bit vector properly
// ------------------------------------------------------------------------------------------
var lineCount = 0
Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  val sp = line.split("\t")
  val events = line.split("\t").slice(line.split("\t").size - (eventCount),line.split("\t").size)
  val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount*2),line.split("\t").size - eventCount)
  annotations.write("lv" + sp(0) + "id" + lineCount + "\t" + sp(0) + "\t" + sp(1) + "\t" + eventStrings.mkString("-") + "\n")

  events.zipWithIndex.foreach{case(eventVal,index) =>
    if (eventVal.toInt > maxByPos(index))
      maxByPos(index) = eventVal.toInt
  }
  lineCount += 1
}}
annotations.close()


// ------------------------------------------------------------------------------------------
// we used to use other metrics, here just use the count
val encodingLengths = maxByPos.map{max => max + 1} // math.ceil(math.log10(max.toDouble)/math.log10(2.0)).toInt}.toArray
val indicies = encodingLengths.sum

// ------------------------------------------------------------------------------------------
// now process the input file and make a PARS file
// ------------------------------------------------------------------------------------------
var newlineCount = 0
output.write(lineCount + "\t" + encodingLengths.sum + "\n")
Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount),line.split("\t").size)
  output.write(( "lv" + line.split("\t")(0) + "id" + newlineCount).padTo(10," ").mkString + " " + toOutputString(encodingLengths,eventStrings) + "\n")
  newlineCount += 1
}}
output.close()
