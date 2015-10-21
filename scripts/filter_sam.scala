/**
  *
  * given a sam file, filter out cigars that don't match the read: our S,I, and M values
  * should sum to the length of the read
  *
  */

import scala.io._
import java.io._
import java.io._
import scala.util.matching.Regex
import scala.collection._

val cigarRegex = new Regex("(\\d+)([MIDSNHP=X])")

// -----------------------------------------------------------------------------
// parameters from the command line
val inputSAM = Source.fromFile(args(0)).getLines()
val outputSAM = new PrintWriter(args(1))
val reference = args(2)

// read in the header on the fasta, get the name, and get the length of the reference
val ref = Source.fromFile(reference).getLines()
val headerLine = ref.next().stripPrefix("> ").stripPrefix(">")
val refLength = ref.map{ln => ln.length}.sum

var mismatched = 0
var matched = 0
// -----------------------------------------------------------------------------
// find all the barcodes we have, and output the cut-down read to the output file
inputSAM.foreach{line => {
  if (line startsWith "@HD") {
    outputSAM.write("@HD\tVN:1.5\n@SQ\tSN:" + headerLine + "\tLN:" + refLength + "\n")
  } else if (line startsWith "@") {
    outputSAM.write(line + "\n")
  } else {
    val sp = line.split("\t")
    var length = 0

    if (sp(5) == "*") {
      //outputFailedSAM.write(sp(0) + "\t" + length + "\t" + sp(9).length + "\n")
    } else {
      val matches = cigarRegex findAllIn sp(5)

      matches.matchData.foreach{matching => matching.subgroups(1) match {
        case "M" => length += matching.subgroups(0).toInt
        case "I" => length += matching.subgroups(0).toInt
        case "S" => length += matching.subgroups(0).toInt
        case _ => { /* do nothing */}
      }}

      if (length == sp(9).length) {
        outputSAM.write(line + "\n")
        matched += 1
      } else {
        //outputFailedSAM.write(sp(0) + "\t" + length + "\t" + sp(9).length + "\n")
        mismatched += 1
      }
    }
  }
}}

outputSAM.close()
println("matched " + matched)
println("mismatched " + mismatched)
