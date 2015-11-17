import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._

// check that the correct number of inputs was provided
val inputSam = Source.fromFile(args(0)).getLines()
val inputFasta= Source.fromFile(args(1)).getLines()

val inputFastqName = inputFasta.next().stripPrefix(">")
val inputFastqSeq = inputFasta.map{sp => sp}.mkString("")

val outputSam = new PrintWriter(args(2))

// output the header
outputSam.write("@HD\tVN:1.5\tSO:coordinate\n")
outputSam.write("@SQ\tSN:" + inputFastqName + "\tLN:" + inputFastqSeq.length + "\n")

// for each line process it
inputSam.foreach{line => {
  if (!(line startsWith "@")) {
    outputSam.write(line + "\n")
  }
}}
outputSam.close()
