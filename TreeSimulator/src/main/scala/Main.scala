package main.scala


import main.scala.TreeSimulation.{GlobalSiteModel, TreeNode, GlobalSiteModel$}
import org.apache.commons.math3.distribution.PoissonDistribution

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
                  phylogenyOutput: File = new File(Main.NOTAREALFILENAME),
                  bestTree: File = new File(Main.NOTAREALFILENAME),
                  annotationOutput: File = new File(Main.NOTAREALFILENAME),
                  targetSiteCount: Int = 10,
                  numberOfSimulations: Int = 100,
                  oldStyleAnalysis: Boolean = true
                   )


object Main extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[Config]("TreeSimluator") {
    head("UMIMerge", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("meltedUMIFile") required() valueName ("<file>") action { (x, c) => c.copy(meltedUMIFile = x) } text ("the UMI summary file")
    opt[File]("bestTree") required() valueName ("<file>") action { (x, c) => c.copy(bestTree = x) } text ("where to put the best tree output for our phylogeny")


    // some general command-line setup stuff
    note("processes reads with UMIs into merged reads\n")
    help("help") text ("prints the usage information you see here")
  }

  parser.parse(args, Config()).map {
    config: Config => {

      // the number of markov chains we hope to simulate
      val numberOfChainsPerEventPerGeneration = 10

      // the poisson mean
      val poissonMean = 1.2
      val normalizedPoissonMean =.1 / 10.0

      // setup a global model file -- the site specific probabilities for each event over sites in the target sequence
      val (model, events) = MulitMarginalDist.createMode(config.meltedUMIFile, config.targetSiteCount)

      // now sample some of possible chains, given our events, our model, and the number of chains
      val timing = new PoissonDistribution(poissonMean)

      // a collection of event chains, find the most likely event chain for each event we have
      val chainsPerEvent = new HashMap[Int, ArrayBuffer[RootPath]]()

      (0 until numberOfChainsPerEventPerGeneration * events.size) foreach { eventID =>
        val id = eventID % numberOfChainsPerEventPerGeneration
        val buffer = chainsPerEvent.getOrElse(id, new ArrayBuffer[RootPath]())
        buffer += RootPath(events(id), timing, model)
        chainsPerEvent(id) = buffer
      }

      // now go and find the most likely chain per Event
      val bestChain = new Array[RootPath](numberOfChainsPerEventPerGeneration)

      chainsPerEvent.foreach { case (id, eventBuffer) => {
        var probString = ""
        val eventB = eventBuffer.toArray
        var bestEvent: Option[RootPath] = None
        eventB.foreach { evt => {
          probString += evt.calculateProbability(normalizedPoissonMean) + ","
          if (!bestEvent.isDefined || evt.calculateProbability(normalizedPoissonMean) > bestEvent.get.calculateProbability(normalizedPoissonMean)) {
            bestEvent = Some(evt)
          }
        }
        }
        println(bestEvent.get.calculateProbability(normalizedPoissonMean) + "\t" + bestEvent.get.toEventString + "\t" + probString)
      }
      }
    }

    // our output files
  }