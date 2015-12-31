package main.scala.simulator

import main.scala.Config
import main.scala.DistanceMatrix
import main.scala.InputTable
import main.scala.Main._
import main.scala.NormalizedEventCounter
import main.scala.StatsFile
import main.scala.SumLogDistance

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
case class Config(inputReadsTable: File = new File(Simulation.NOTAREALFILENAME),
                  simulationOutput: File = new File(Simulation.NOTAREALFILENAME))

object Simulation extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[Config]("TreeSimluator") {
    head("UMIMerge", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("meltedUMIFile") required() valueName ("<file>") action { (x, c) => c.copy(inputReadsTable = x) } text ("(input) the UMI summary file")
    opt[File]("simulationOutput") required() valueName ("<file>") action { (x, c) => c.copy(simulationOutput = x) } text ("(output) out simulation results")

    // some general command-line setup stuff
    note ("take in a input table, simulate trees from that data, and output summary stastics\n")
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
      var statsFile : Option[InputTable] = None
      if (config.inputReadsTable.getAbsolutePath endsWith ".stats")
        statsFile = Some(StatsFile(config.inputReadsTable))
      else
        statsFile = Some(CallFile(config.inputReadsTable))

      // ------------------------------------------------------------------------------------------------------------------------
      // events to counts -- count the total events over all positions
      // ------------------------------------------------------------------------------------------------------------------------
      val eventCounter = new NormalizedEventCounter(statsFile.get,200000.0 /* 200K */)

      // ------------------------------------------------------------------------------------------------------------------------
      // setup a distance matrix
      // ------------------------------------------------------------------------------------------------------------------------
      val distances = new DistanceMatrix(statsFile.get.getAllEvents(), SumLogDistance(eventCounter, 1, 2))
      distances.toDistanceFile(config.distanceMatrixFile)

      // now start simulating events --

    }
  }

  // our output files
}