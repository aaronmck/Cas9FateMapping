import scala.io._
import java.io._

// given a filter criteria (a simple grep-like function on the line),
// find reads from a stats file that match the critera and align them to the reference
//
val lines = Source.fromFile(args(0)).getLines()
val filterString = args(1)
val reference = args(2)
val output = new PrintWriter(args(3))

// process the reference
val refString = Source.fromFile(reference).getLines().drop(1).mkString("")

// read the header line, and process out the token we'll be interested in
val cigarToken = "cigar"
val readToken = "final.consensus"
val header = lines.next().split("\t")
val cigarPos = header.indexOf(cigarToken)
val readPos = header.indexOf(readToken)
if (cigarPos < 0 || readPos < 0)
  throw IllegalStateException("Unable to find either the cigar token or the read token")

output.write()
lines.foreach{line => {
  if (line contains filterString) {
    val sp = line.split("\t")
    val read = sp(readPos)
    val cigar = sp(cigarPos)

    
  }
}}
