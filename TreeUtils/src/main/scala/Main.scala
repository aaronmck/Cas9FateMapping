package main.scala


import beast.util._
import beast.evolution.tree._
import scala.collection.JavaConversions._
import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._
import scala.util.Random

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
case class TreeConfig(inputTree: File = new File(Main.NOTAREALFILENAME),
                      inputGenotypes: File = new File(Main.NOTAREALFILENAME),
                      inputAnnotations: File = new File(Main.NOTAREALFILENAME),
                      inputSampleToClade: File = new File(Main.NOTAREALFILENAME),
                      inputEventsToNumbers: File = new File(Main.NOTAREALFILENAME),
                      eventsToColors: File = new File(Main.NOTAREALFILENAME),
                      outputTree: File = new File(Main.NOTAREALFILENAME))


object Main extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[TreeConfig]("UMIMerge") {
    head("TreeUtils", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("inputTree") required() valueName ("<file>") action { (x, c) => c.copy(inputTree = x) } text ("first read file ")
    opt[File]("inputGenotypes") required() valueName ("<file>") action { (x, c) => c.copy(inputGenotypes = x) } text ("first read file ")
    opt[File]("inputAnnotations") required() valueName ("<file>") action { (x, c) => c.copy(inputAnnotations = x) } text ("first read file ")
    opt[File]("inputSampleToClade") required() valueName ("<file>") action { (x, c) => c.copy(inputSampleToClade = x) } text ("first read file ")
    opt[File]("inputEventsToNumbers") required() valueName ("<file>") action { (x, c) => c.copy(inputEventsToNumbers = x) } text ("first read file ")
    opt[File]("eventsToColors") valueName ("<file>") action { (x, c) => c.copy(eventsToColors = x) } text ("first read file ")
    opt[File]("outputTree") required() valueName ("<file>") action { (x, c) => c.copy(outputTree = x) } text ("second reads file")

    // some general command-line setup stuff
    note("processes reads with UMIs into merged reads\n")
    help("help") text ("prints the usage information you see here")
  }

  // *********************************** Run *******************************************************
  parser.parse(args, TreeConfig()) map { config => {
    // mixTrees: File, mixOutput: File, annotations: File, sampleToClade: File, eventsToNumbers: String
    val parser = new ParsimonyProcessor(config.inputTree,
      config.inputGenotypes,
      config.inputAnnotations,
      config.inputSampleToClade,
      config.inputEventsToNumbers,
      if (config.eventsToColors.getAbsolutePath contains NOTAREALFILENAME) None else Some(config.eventsToColors),
      config.outputTree)
  }} getOrElse {
    println("Unable to parse the command line arguments you passed in, please check that your parameters are correct")
  }

}
