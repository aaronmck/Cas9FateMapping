// load up a stats file and produce a distance matrix using a simple scoring scheme
import scala.io._
import java.io._
import scala.collection.mutable._
import scala.math._

val output = new PrintWriter(args(1))
val eventCount = args(2).toInt
val annotationOutput = new PrintWriter(args(3))

def toOutputString(bitSizes: Array[Int], eventStrings: Array[String]) = eventStrings.zipWithIndex.map{case(evt,index) => toBitString(evt,bitSizes(index))}.mkString("")
def toBitString(evt: String, length: Int): String = {
  val ret = new Array[Int](length)
  ret(evt.toInt) = 1
  return(ret.mkString(""))
}

val maxByPos = new Array[Int](eventCount)

var lineCount = 0
Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount),line.split("\t").size)
  eventStrings.zipWithIndex.foreach{case(eventVal,index) =>
    if (eventVal.toInt > maxByPos(index))
      maxByPos(index) = eventVal.toInt
  }
  if ((line.split("\t")(0) contains "embryos_1_1") || (line.split("\t")(0) contains "embryos_1_6"))
    lineCount += 1
}}

val encodingLengths = maxByPos.map{max => max + 1} // math.ceil(math.log10(max.toDouble)/math.log10(2.0)).toInt}.toArray

// strip down a stats file to a set of events
// output.write("sample_umi\tbitarray\n")
output.write("#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=" + lineCount + " NCHAR=" + encodingLengths.sum + ";\nFORMAT DATATYPE=Standard SYMBOLS= \"0 1 2 3 4 5 6 7 8 9\" MISSING=? GAP= -;\nMATRIX\n")

//output.write("#NEXUS\nbegin data;\ndimensions ntax=" + lineCount + " nchar=" + encodingLengths.sum + ";\nformat datatype=dna interleave=yes gap=- missing=?;\nmatrix\n")
//output.write((lineCount+1) + " " + encodingLengths.sum + "\n")
// output.write((lineCount) + " " + encodingLengths.sum + "\n")
//output.write("ROOT_UNKNOWN\t" + toOutputString(encodingLengths,encodingLengths.map{st => "0".toString}.toArray) + "\n")
val indicies = encodingLengths.sum
annotationOutput.write("taxa\tembryo\tumi\n")
//output.write("sample\t" + (0 until indicies).map{ln => "index" + ln}.mkString("\t") + "\n")
Source.fromFile(args(0)).getLines().drop(1).foreach{line => {
  val eventStrings = line.split("\t").slice(line.split("\t").size - (eventCount),line.split("\t").size)
  if ((line.split("\t")(0) contains "embryos_1_1") || (line.split("\t")(0) contains "embryos_1_6"))
    output.write(line.split("\t")(0) + "_" + line.split("\t")(1) + "\t" + toOutputString(encodingLengths,eventStrings) + "\n")
    annotationOutput.write(line.split("\t")(0) + "_" + line.split("\t")(1) + "\t" + line.split("\t")(0) + "\t" + line.split("\t")(1) + "\n")
    //output.write(line.split("\t")(0) + "\t" + toOutputString(encodingLengths,eventStrings) + "\n")
}}
output.write(";\nEND;")
//output.write(";\nend;\n")
output.close()
annotationOutput.close()
