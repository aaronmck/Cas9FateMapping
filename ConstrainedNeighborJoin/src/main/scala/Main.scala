package main.scala


import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._

/**
 *
 * Copyright (c) 2015, aaronmck
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2.  Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */
case class Config(meltedUMIFile: File = new File(Main.NOTAREALFILENAME),
                  newickFile: File = new File(Main.NOTAREALFILENAME),
                  mergeInformation: File = new File(Main.NOTAREALFILENAME),
                  distanceMatrixFile: File = new File(Main.NOTAREALFILENAME))

object Main extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[Config]("TreeSimluator") {
    head("UMIMerge", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("meltedUMIFile") required() valueName ("<file>") action { (x, c) => c.copy(meltedUMIFile = x) } text ("(input) the UMI summary file")
    opt[File]("newickFile") required() valueName ("<file>") action { (x, c) => c.copy(newickFile = x) } text ("(output) our output Newick file")
    opt[File]("mergeInformation") required() valueName ("<file>") action { (x, c) => c.copy(mergeInformation = x) } text ("(output) annotations for the output tree (used in FigTree)")
    opt[File]("distanceMatrixFile") required() valueName ("<file>") action { (x, c) => c.copy(distanceMatrixFile = x) } text ("(output) our distance matrix file")

      // some general command-line setup stuff
      note ("processes reads with UMIs into merged reads\n")
      help ("help") text ("prints the usage information you see here")
  }

  parser.parse(args, Config()).map {
    config: Config => {
      val outputAnnotations = args(2)
      val countsOutput = args(3)
      val useConstraints = false

      // ------------------------------------------------------------------------------------------------------------------------
      // parse out the event strings into events
      // ------------------------------------------------------------------------------------------------------------------------
      println("reading in event file...")
      val statsFile = StatsFile(config.meltedUMIFile)

      // ------------------------------------------------------------------------------------------------------------------------
      // events to counts -- count the total events over all positions
      // ------------------------------------------------------------------------------------------------------------------------
      val eventCounter = new NormalizedEventCounter(statsFile,200000.0 /* 200K */)
      // ------------------------------------------------------------------------------------------------------------------------
      // setup a distance matrix
      // ------------------------------------------------------------------------------------------------------------------------
      val distances = new DistanceMatrix(statsFile.allEvents, SumLogDistance(eventCounter, 1, 2))
      distances.toDistanceFile(config.distanceMatrixFile)

      // ------------------------------------------------------------------------------------------------------------------------
      // find the minimum set of events
      // ------------------------------------------------------------------------------------------------------------------------
      println("Performing merges (dot per merge)")
      val minSet = distances.minimizeSet(useConstraints)

      // ------------------------------------------------------------------------------------------------------------------------
      // output the remaining nodes as tree
      // ------------------------------------------------------------------------------------------------------------------------
      val newickFileOutput = new PrintWriter(config.newickFile)
      val nodeStats = new PrintWriter(config.mergeInformation)

      nodeStats.write("taxa\tname\tsample\tdepth\tnumberOfReads\tevents\n")

      minSet.foreach{node => {
        newickFileOutput.write("(" + node.newickString(1.0, 0.0, nodeStats) + ")root:1.0;\n")
      }}

      newickFileOutput.close()
      nodeStats.close()
    }
  }

  // our output files
}