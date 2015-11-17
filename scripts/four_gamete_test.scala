// 4 gamete test
import scala.io._
import java.io._
import scala.collection.mutable._

// val output = new PrintWriter(args(1))

val lines = Source.fromFile(args(0)).getLines()
val header = lines.next()

// setup the number of columns we expect
var columns = Array[ArrayBuilder[Int]]()
header.split("\t").drop(1).foreach{headerItem => {
  columns :+= ArrayBuilder.make[Int]
}}

lines.foreach{line => {
  line.split("\t").drop(1).zipWithIndex.foreach{case(token,index) => columns(index) += token.toInt}
}}

val outputFile = new PrintWriter(args(1))
println("processed " + columns.size + " columns")

def compareVectOneTwo(vecOne: Array[Int], vecTwo: Array[Int]): Array[Int] = {
  val result = Array[Int](0,0,0,0,0)
  vecOne.zip(vecTwo).foreach{case(indOne,indTwo) => {
    if (indOne == 0 && indTwo == 0) {
      result(0) += 1
    } else if (indOne == 0 && indTwo == 1) {
      result(1) += 1
    } else if (indOne == 1 && indTwo == 0) {
      result(2) += 1
    } else if (indOne == 1 && indTwo == 1) {
      result(3) += 1
    } else {
      result(4) += 1
    }
  }}
  result
}

var resultColumns = Array[Array[Int]]()
columns.foreach{col => {
  val rst = col.result
  if (rst.sum > 1)
    resultColumns :+= rst
}}

outputFile.write("i\tj\t0.0\t0.1\t1.0\t1.1\tfour.gamete\n")
var counts = 0
var allGTzero = 0
for (i <- 0 until resultColumns.size) {
  for (j <- 0 until resultColumns.size) {
    if (i < j) {
      val comp = compareVectOneTwo(resultColumns(i),resultColumns(j)).slice(0,4)
      counts += 1
      val allValues = comp.map{cp => if (cp == 0) 0 else 1}.sum >= 4
      outputFile.write(i + "\t" + j + "\t" + comp.mkString("\t") + "\t" + allValues + "\n")

      if (counts % 1000000 == 0)
        println("Count = " + counts)
      if (allValues)
        allGTzero += 1
    }
  }
}
println("all GT zero = " + allGTzero + " of " + counts)
outputFile.close()
