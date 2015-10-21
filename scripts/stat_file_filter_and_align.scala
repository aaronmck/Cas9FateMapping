import scala.io._
import java.io._

// given a filter criteria (a grep of the line), find reads from a stats file that
// match the critera and align them to the reference
//
val lines = Source.fromFile(args(0))
val sequenceFile = Source.fromFile(args(1))
val output = new PrintWriter(args(2))
val filterString = args(3)



output.write()
lines.foreach{line => {
  if (!(line startsWith "@")) {

  }
}}
