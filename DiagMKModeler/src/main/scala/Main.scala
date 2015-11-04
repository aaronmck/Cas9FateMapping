package main.scala


import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._

/**
 * created by aaronmck on 2/13/14
 *
 * Copyright (c) 2014, aaronmck
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
                  phylogenyOutput: File = new File(Main.NOTAREALFILENAME),
                  targetSiteCount: Int = 10
                   )


object Main extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[Config]("UMIMerge") {
    head("UMIMerge", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("meltedUMIFile") required() valueName ("<file>") action { (x, c) => c.copy(meltedUMIFile = x) } text ("the UMI summary file")
    opt[File]("phylogenyOutput") required() valueName ("<file>") action { (x, c) => c.copy(phylogenyOutput = x) } text ("where to put the output for our phlogeny")
    opt[Int]("targetSiteCount") required() action { (x, c) => c.copy(targetSiteCount = x) } text ("the length of our UMIs")


    // some general command-line setup stuff
    note("processes reads with UMIs into merged reads\n")
    help("help") text ("prints the usage information you see here")
  }

  // *********************************** Main script stuff ******************************************************
  parser.parse(args, Config()) map {
    config: Config => {

      // ------------------------------------------------------------------------------------
      // Main script
      // ------------------------------------------------------------------------------------

      val positionToCounts = HashMap[Int, HashMap[String, Int]]()
      val positionToCountsNormalized = HashMap[Int, HashMap[String, Double]]()
      var totalRows = 0
      val eventNumberToEvents = HashMap[Int, HashMap[String, Array[Event]]]()

      (0 until config.targetSiteCount + 1).foreach { i => {
        positionToCounts(i) = new HashMap[String, Int]()
        positionToCountsNormalized(i) = new HashMap[String, Double]()
        eventNumberToEvents(i) = HashMap[String, Array[Event]]()
      }
      }

      // ------------------------------------------------------------------------------------
      // now normalize the positional counts by the total
      // ------------------------------------------------------------------------------------
      positionToCounts.foreach { case (position, counts) => {
        counts.foreach { case (event, count) => positionToCountsNormalized(position)(event) = count.toDouble / totalRows.toDouble }
      }
      }

      println("unique tags:")
      (1 until 10).foreach { index => println(index + " " + eventNumberToEvents(index).size) }

      val output = new PrintWriter(config.phylogenyOutput)
      output.write(Event.headerString + "\n")

      (2 until 10).foreach { index =>
        if (eventNumberToEvents contains index) {
          val befores = eventNumberToEvents(index - 1).map { case (id, hits) => hits(0) }.toArray
          val firsts = eventNumberToEvents(index).map { case (id, hits) => hits(0) }.toArray

          print("index " + index + " count " + firsts.length + " with prob = ")

          var tt = 0
          firsts.map { case (evt) => {
            if (evt.oneEdit(befores))
              tt += 1
            output.write(evt.toOutputString + "\n")
          }
          }
          println(tt)
        }
      }
      output.close()

    }
  }

  // our output files
}