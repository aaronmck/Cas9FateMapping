// load up a stats file and produce a distance matrix using a simple scoring scheme
import scala.io._
import java.io._
import scala.collection.mutable._
import scala.math._

// inputs
// 1) input melted UMI data
// 2) output file
// 3) counts of events (targets)
// 4) annotations output

val sep = "\t"
val output = new PrintWriter(args(1))
val eventCount = args(2).toInt
val annotationOutput = new PrintWriter(args(3))

def toOutputString(bitSizes: Array[Int], eventStrings: Array[String]) = eventStrings.zipWithIndex.map{case(evt,index) => toBitString(evt,bitSizes(index))}.mkString(sep)
def toBitString(evt: String, length: Int): String = {
  val ret = new Array[Int](length)
  ret(evt.toInt) = 1
  return(ret.mkString(sep))
}

val maxByPos = new Array[Int](eventCount)

def keepLine(line: String) = true // ((line startsWith "embryos_1") && (line contains "PASS"))

println("Processing input lines to find the maximum event count per position...")
var lineCount = 0
Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  if (keepLine(line)) {
    val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount),line.split("\t").size)
    eventStrings.zipWithIndex.foreach{case(eventVal,index) =>
      if (eventVal.toInt > maxByPos(index))
        maxByPos(index) = eventVal.toInt
    }
    lineCount += 1
  }
}}

val encodingLengths = maxByPos.map{max => max + 1} // math.ceil(math.log10(max.toDouble)/math.log10(2.0)).toInt}.toArray

val indicies = encodingLengths.sum
annotationOutput.write("taxa\tembryo\tumi\n")

val outputAlready2 = new HashMap[String,Boolean]()
val outputAlready = new HashMap[String,Boolean]()

println("Reprocessing...")
//output.write("sample\t" + (0 until indicies).map{ln => "index" + ln}.mkString("\t") + "\n")
Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  if (keepLine(line)) {
  val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount),line.split("\t").size)

  val outString = toOutputString(encodingLengths,eventStrings)
  outputAlready2(line.split("\t")(0) + outString) = true
}
}}

output.write("sample\t" + ((0 until encodingLengths.sum).map{st => "target" + st}).mkString("\t") + "\n")
println("Final output...")

Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  if (keepLine(line)) {
    val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount),line.split("\t").size)
    val outString = toOutputString(encodingLengths,eventStrings)
    val name = line.split("\t")(0) + "_" + line.split("\t")(1)

    if (!(outputAlready contains (line.split("\t")(0) + outString))) {
      output.write(name + "\t" + outString + "\n")
      annotationOutput.write(line.split("\t")(0) + "_" + line.split("\t")(1) + "\t" + line.split("\t")(0) + "\t" + line.split("\t")(1) + "\n")
      //output.write(line.split("\t")(0) + "\t" + toOutputString(encodingLengths,eventStrings) + "\n")
    }
    outputAlready(line.split("\t")(0) + outString) = true
  }
}}
//output.write(";\nEND;")
//output.write(";\nend;\n")
output.close()
annotationOutput.close()
