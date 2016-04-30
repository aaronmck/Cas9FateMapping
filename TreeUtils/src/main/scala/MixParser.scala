package main.scala

import scala.io.Source

/**
  * Parse the output from the PHYLIP package MIX and create
  * event strings for each of the internal nodes
  */
class MixParser(mixOutput: String, eventsToNumbers: String) {
  // the header line we're looking for in the file is:
  // From    To     Any Steps?    State at upper node

  // a mapping from the input tree number to it's annotations

  // first parse out the events to number data, and make a look-up table



  val inputFile = Source.fromFile(mixOutput).getLines()

  var seenHeader = false
  inputFile.foreach{line => {

  }}
}
