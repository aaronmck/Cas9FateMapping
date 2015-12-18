import scala.io._
import java.io._
import scala.collection.mutable._

// parse out the event length and type
def eventToLength(event:String): Tuple3[String,Int,Int] = {
    val sp = event.split("\\+")
    val typ = sp(0).slice(sp(0).length-1,sp(0).length)
    val length = sp(0).slice(0,sp(0).length -1).toInt
    val pos = sp(1).toInt
    return (typ,length,pos)
}

// a function to get the calls file, and convert to edit lengths and event lengths
def callLineToLengths(tokens: Array[String],sample: String, outputFile: PrintWriter) {
    val eventCounts = new HashMap[String,Int]()
    val eventLengths = new HashMap[String,Tuple3[String,Int,Int]]()

    tokens.foreach{x => {
        if (x != "NONE" && x != "UNKNOWN" && !(x contains "WT")) {
          x.split("&").foreach{subX => {
            eventLengths(subX) = eventToLength(subX)
            eventCounts(subX) = eventCounts.getOrElse(subX,0) + 1
          }}
        }
    }}

    eventLengths.foreach{case(event,count) =>
        outputFile.write(count._1 + "\t" + count._2 + "\t" + count._3 + "\t" + eventCounts(event) + "\t" + sample + "\n")
    }
}

def processSample(fileName:String, output: PrintWriter,sample: String) {
  val input18X = Source.fromFile(fileName).getLines()
  val header = input18X.next()
  val token1Pos = 9
  val token10Pos = 19

  input18X.zipWithIndex.foreach{case(line,ind) => {
    if ((line contains "PASS") && !(line contains "WT")) {
    val sp = line.split("\t")
    callLineToLengths(sp.slice(token1Pos,token10Pos),sample,output)
    if (ind % 100000 == 0)
      println(ind)
    }
  }}
}

val output = new PrintWriter("event_counts.txt")
output.write("type\teventLength\teventPosition\tsiteLength\tsample\n")

val t1X = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_11_23_Deep_Sequence_Dilution/data/pipeline_output/embryo_1X_1/embryo_1X_1.calls"
val t13X = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_11_23_Deep_Sequence_Dilution/data/pipeline_output/embryo_1.3X_1/embryo_1.3X_1.calls"
val t16X = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_11_23_Deep_Sequence_Dilution/data/pipeline_output/embryo_1.6X_6/embryo_1.6X_6.calls"
val t18X = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_11_23_Deep_Sequence_Dilution/data/pipeline_output/embryo_1.18X_6/embryo_1.18X_6.calls"

processSample(t1X,output,"x1")
processSample(t13X,output,"x1/3")
processSample(t16X,output,"x1/6")
processSample(t18X,output,"x1/18")

output.close()
