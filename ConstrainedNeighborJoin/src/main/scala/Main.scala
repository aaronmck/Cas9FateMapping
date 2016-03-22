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
case class Config(inputReadsTable: File = new File(Main.NOTAREALFILENAME),
                  // newickFile: File = new File(Main.NOTAREALFILENAME),
                  annotationFile: File = new File(Main.NOTAREALFILENAME),
                  noneDistance: Double = 1.0,
                  equalScore: Double = 1.0,
                  distanceMatrixFile: File = new File(Main.NOTAREALFILENAME))

object Main extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[Config]("TreeSimluator") {
    head("UMIMerge", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("meltedUMIFile") required() valueName ("<file>") action { (x, c) => c.copy(inputReadsTable = x) } text ("(input) the UMI summary file")
    //opt[File]("newickFile") required() valueName ("<file>") action { (x, c) => c.copy(newickFile = x) } text ("(output) our output Newick file")
    opt[File]("annotations") required() valueName ("<file>") action { (x, c) => c.copy(annotationFile = x) } text ("(output) annotations for the output tree (used in FigTree)")
    opt[File]("distanceMatrixFile") required() valueName ("<file>") action { (x, c) => c.copy(distanceMatrixFile = x) } text ("(output) our distance matrix file")
    opt[Double]("noneDistance") valueName ("<file>") action { (x, c) => c.copy(noneDistance = x) } text ("the distance for two wild-type cut sites when comparing events")
    opt[Double]("equalScore") valueName ("<file>") action { (x, c) => c.copy(equalScore = x) } text ("the distance for two identical edits at cut sites when comparing events")

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
      var statsFile : InputTable = StatsFile(config.inputReadsTable)

      // ------------------------------------------------------------------------------------------------------------------------
      // events to counts -- count the total events over all positions
      // ------------------------------------------------------------------------------------------------------------------------
      val eventCounter = new NormalizedEventCounter(statsFile, 200000.0 /* 200K */)

      var distMetrics = new HashMap[String,DistanceMetric[IndexedNode]]()
      distMetrics("SumLog") = SumLogDistance(eventCounter, config.noneDistance, config.equalScore, 2)
      distMetrics("AvgLog") = AvgLogDistance(eventCounter, config.noneDistance, config.equalScore, 2)
      distMetrics("Additive") = SimpleAdditiveDistance(eventCounter)

      // ------------------------------------------------------------------------------------------------------------------------
      // setup a distance matrix
      // ------------------------------------------------------------------------------------------------------------------------
      //val distMetric = SumLogDistance(eventCounter, config.noneDistance, config.equalScore, 2)
      //val distMetric = SimpleAdditiveDistance(eventCounter)
      val distMetric = Hamming(eventCounter)

      //val distMetric = Hamming(eventCounter)
      //val distMetric = SumLogDistance(eventCounter, config.noneDistance, config.equalScore, 2)
      val distances = new DistanceMatrix(statsFile.getUniqueEvents(), distMetric)
      distances.toDistanceFile(config.distanceMatrixFile)

      // ------------------------------------------------------------------------------------------------------------------------
      // find the minimum set of events
      // ------------------------------------------------------------------------------------------------------------------------
      //println("Performing merges (dot per merge)")
      //val minSet = distances.minimizeSet(useConstraints)

      distances.toAnnotationFile(config.annotationFile, statsFile)
      // ------------------------------------------------------------------------------------------------------------------------
      // output the remaining nodes as tree
      // ------------------------------------------------------------------------------------------------------------------------
      /*
      val newickFileOutput = new PrintWriter(config.newickFile)
      val nodeStats = new PrintWriter(config.mergeInformation)

      nodeStats.write("taxa\tname\tsample\tdepth\tnumberOfReads\tevents\n")

      minSet.foreach{node => {
        newickFileOutput.write("(" + node.newickString(1.0, 0.0, nodeStats) + ");\n")
      }}

      newickFileOutput.close()
      nodeStats.close()
      */
    }
  }

  // our output files
}