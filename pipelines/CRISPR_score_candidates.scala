/** *********************************************
  * Requires Java 1.7. On the UW machines this means
  * running at least: "module load java/7u17", or the appropriate
  * module if there's a more recent version, to load java 1.7 into
  * your environment.  The best idea is to place this into your
  * bash/shell profile
  *
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
 * give a BED file containing sites in a genome of interest, find off-target hits, and score each
 * guide for it's on and off target abilities.
 *
 * The input file:
 *
 * A bed file, where the 4th column (usually a name or something) with the target sequence
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
  var tallyJar: File = new File("/net/gs/vol1/home/aaronmck/source/DeepFry/target/scala-2.11/CRISPREngineered-assembly-1.0.jar ")

  @Input(doc = "a temporary location to write intermediate files", fullName = "temp", shortName = "temp", required = true)
  var tmp: File = "/tmp/"

  @Input(doc = "the reference to scan against, split up into bins of potential off-targets", fullName = "refBed", shortName = "ref", required = false)
  var genomeFile: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_07_30_GenomeScans/data/human_bins/global.unique.txt.gz"

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

        add(FindCandidates(inpF,outF, genomeFile))
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
  case class FindCandidates(inputTargets: File, candidateList: File, genomeBed: File) extends ExternalCommonArgs {
    @Input(doc = "the input target file")      var input = inputTargets
    @Input(doc = "the genome bed file")        var gBed = genomeBed
    @Output(doc = "the output file to write ") var output = candidateList

    def commandLine = "java -Xmx3g -jar " + tallyJar + " --analysis tally --targetBed " + input + " --output " + output + " --genomeBed " + gBed

    this.analysisName = queueLogDir + output + ".gather"
    this.jobName = queueLogDir + output + ".gather"
    this.memoryLimit = 4

    /*
    --analysis tally \
    --targetBed /net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_07_30_GenomeScans/data/all3_target_set_nh.bed \
    --output /net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_07_30_GenomeScans/data/50K_scores_every_interation_guide_only_zebafish.bed \
    --genomeBed /net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_07_30_GenomeScans/data/fish_bins/global.unique.txt.gz &
    */
  }

  // ********************************************************************************************************
  case class ScoreCandidates(inputRawOTs: File, scoredBed: File) extends ExternalCommonArgs {
    @Input(doc = "the input target file")      var input = inputRawOTs
    @Output(doc = "the output file to write ") var output = scoredBed

    def commandLine = "java -Xmx3g -jar " + tallyJar + " --analysis score --inputBed " + input + " --outputBed " + output

    this.analysisName = queueLogDir + output + ".score"
    this.jobName = queueLogDir + output + ".score"
    this.memoryLimit = 4
  }

  // ********************************************************************************************************
  case class JoinScores(inputScores: List[File], finalScoredFile: File) extends ExternalCommonArgs {
    @Input(doc = "the input score files")      var inputFiles = inputScores
    @Output(doc = "the output file to write ") var output = finalScoredFile

    def commandLine = "cat " + inputFiles.mkString(" ") + " > " + output

    this.analysisName = queueLogDir + output + ".finalScores"
    this.jobName = queueLogDir + output + ".finalScores"
    this.memoryLimit = 4
  }
}
