/** *********************************************
  * Requires Java 1.7. On the UW machines this means
  * running at least: "module load java/7u17", or the appropriate
  * module if there's a more recent version, to load java 1.7 into
  * your environment.  The best idea is to place this into your
  * bash/shell profile
  *
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
  * @author Aaron
  * @date Decemberish, 2015
  *
  */

package org.broadinstitute.gatk.queue.qscripts

import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.function.ListWriterFunction
import org.broadinstitute.gatk.utils.commandline.Hidden
import org.broadinstitute.gatk.utils.commandline
import org.broadinstitute.gatk.engine.arguments.ValidationExclusion

import collection.JavaConversions._
import htsjdk.samtools.SAMFileReader
import htsjdk.samtools.SAMFileHeader.SortOrder
import scala.io.Source
import scala.collection.immutable._
import java.io._

/**
 * Quality control sequencing data and optionally align the data to the genome
 **/
class ScoreGuides extends QScript {
  qscript =>

  @Input(doc = "Tab-delimited bed file with target sequneces as the 4th column", fullName = "inputBED", shortName = "bed", required = true)
  var inputBED: File = _

  /** **************************************************************************
    * Optional Parameters -- control parameters for the script (alignment, etc)
    * ************************************************************************** */

  @Argument(doc = "number of targets to score per core.  The default (400) means each core takes about a couple hours to score against the whole genome", fullName = "eventsPerCore", shortName = "eventsPerCore", required = false)
  var eventsPerCore: Int = 400

  /** **************************************************************************
    * Path parameters -- where to find tools and aux. files
    * ************************************************************************** */

  @Input(doc = "The path to the tools that compares each guide against the reference and tallies hits", fullName = "tally", shortName = "tally", required = false)
  var tallyScript: File = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/try_full_scan.scala")

  @Input(doc = "The script that takes the tallied events and scores them with on and off target scores", fullName = "score", shortName = "score", required = false)
  var scoringScript: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/../codebase/scripts/score_sites.py"

  @Input(doc = "A bed file containing ALL targets in the genome", fullName = "alltargets", shortName = "alltargets", required = true)
  var allTargets: File = ""

  @Input(doc = "a temporary location to write intermediate files", fullName = "temp", shortName = "temp", required = true)
  var tmp: File = "/tmp/"

  @Input(doc = "a script to join the score files", fullName = "join", shortName = "join", required = false)
  var joinScores: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/../codebase/scripts/join_sites.py"

  @Input(doc = "the reference to scan against", fullName = "reference", shortName = "ref", required = false)
  var reference: File = "/net/shendure/vol10/nobackup/shared/genomes/human_g1k_hs37d5/hs37d5.fa"

  @Output(doc = "the output file of scored events", fullName = "output", shortName = "output", required = true)
  var outputFile: File = ""

  /** **************************************************************************
    * Global Variables
    * ************************************************************************** */

  // Gracefully hide Queue's output -- this is very important otherwise there's Queue output everywhere and it's a mess
  val queueLogDir: String = ".qlog/"


  /** **************************************************************************
    * Main script
    * ************************************************************************** */
  def script() {
      // split their target list into temporary files, each with X events, and run the discovery phase against the whole genome
      var outputFiles = Array[File]()

      // get all of the bed file entries
      val toScoresLines = Source.fromFile(inputBED).getLines().grouped(eventsPerCore)

      toScoresLines.zipWithIndex.foreach{case(collectionOfBedEvents,index) => {
        val inpF = new File(tmp.getAbsolutePath + File.separator + "temp.to.score." + index)
        val outF = new File(tmp.getAbsolutePath + File.separator + "temp.candidates." + index)

        outputFiles :+= outF

        val eventsOutput = new PrintWriter(inpF)
        collectionOfBedEvents.foreach{evtString => eventsOutput.write(evtString + "\n")}
        eventsOutput.close()

        add(FindCandidates(inpF,outF))
      }}

      var totalOutput = List[File]()

      // now score each of the returned events
      outputFiles.zipWithIndex.foreach{case(prevOutput,index) => {
        val scoredOutput = new File(tmp.getAbsolutePath + File.separator + "temp.scored." + index)

        totalOutput :+= scoredOutput
        add(ScoreCandidates(prevOutput,scoredOutput))
      }}

      // and put it all back together into a single scored table
      add(JoinScores(totalOutput,outputFile))
  }

  /** **************************************************************************
    * traits that get tacked onto runnable objects
    * ************************************************************************** */

  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 6 // set this a bit high, there's a weird java / SGE interactions that require higher mem reservations
    this.residentRequest = 6
    this.residentLimit = 6
    //this.isIntermediate = false // by default delete the intermediate files, if you want to keep output from a task set this to false in the case class
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }


  // ********************************************************************************************************
  case class FindCandidates(inputTargets: File, candidateList: File) extends ExternalCommonArgs {
    @Input(doc = "the input target file")      var input = inputTargets
    @Output(doc = "the output file to write ") var output = candidateList

    def commandLine = "scala " + tallyScript + " " + input + " " + output + " " + reference

    this.analysisName = queueLogDir + output + ".gather"
    this.jobName = queueLogDir + output + ".gather"
    this.memoryLimit = 4
  }

  // ********************************************************************************************************
  case class ScoreCandidates(inputTargets: File, scoredFile: File) extends ExternalCommonArgs {
    @Input(doc = "the input target file")      var input = inputTargets
    @Output(doc = "the output file to write ") var output = scoredFile

    def commandLine = "python " + scoringScript + " --input_sites " + input + " --output_sites " + output

    this.analysisName = queueLogDir + output + ".score"
    this.jobName = queueLogDir + output + ".score"
    this.memoryLimit = 4
  }

  // ********************************************************************************************************
  case class JoinScores(inputScores: List[File], finalScoredFile: File) extends ExternalCommonArgs {
    @Input(doc = "the input score files")      var inputFiles = inputScores
    @Output(doc = "the output file to write ") var output = finalScoredFile

    def commandLine = "python " + joinScores + " " + inputFiles.map{inF => " --input " + inF}.mkString("") + " --output " + output

    this.analysisName = queueLogDir + output + ".finalScores"
    this.jobName = queueLogDir + output + ".finalScores"
    this.memoryLimit = 4
  }
}
