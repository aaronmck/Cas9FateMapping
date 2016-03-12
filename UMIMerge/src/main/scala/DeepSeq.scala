package main.scala

import java.util

import _root_.aligner.{Alignment, AlignmentManager}
import _root_.utils.{CutSites}
import main.scala.stats.{StatsContainer, StatsOutput}
import main.scala.utils.{RefReadPair, ReadPair, UnmergedReadParser, ReadPairParser}

import scala.io._
import java.io._
import scala.collection.mutable._
import scala.main.SequencingRead

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
case class DeepConfig(inputFileUnmerged: File = new File(DeepSeq.NOTAREALFILENAME),
                      inputMerged: File = new File(DeepSeq.NOTAREALFILENAME),
                      outputStats: File = new File(DeepSeq.NOTAREALFILENAME),
                      cutSites: File = new File(DeepSeq.NOTAREALFILENAME),
                      primersEachEnd: File = new File(DeepSeq.NOTAREALFILENAME),
                      reference: File = new File(DeepSeq.NOTAREALFILENAME),
                      samplename: String = "TEST")


object DeepSeq extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[DeepConfig]("DeepSeq") {
    head("DeepSeq", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("inputUnmerged") required() valueName ("<file>") action { (x, c) => c.copy(inputFileUnmerged = x) } text ("unmerged reads")
    opt[File]("inputMerged") required() valueName ("<file>") action { (x, c) => c.copy(inputMerged = x) } text ("the merged read file")
    opt[File]("outputStats") required() valueName ("<file>") action { (x, c) => c.copy(outputStats = x) } text ("the output stats file")
    opt[File]("cutSites") required() valueName ("<file>") action { (x, c) => c.copy(cutSites = x) } text ("the location of the cutsites")
    opt[File]("primersEachEnd") required() valueName ("<file>") action { (x, c) => c.copy(primersEachEnd = x) } text ("the file containing the amplicon primers requred to be present, one per line, two lines total")
    opt[String]("sample") required() action { (x, c) => c.copy(samplename = x) } text ("the sample name of this run")

    // some general command-line setup stuff
    note("process aligned reads from non-UMI samples\n")
    help("help") text ("prints the usage information you see here")
  }

  // *********************************** Run *******************************************************
  // run the actual read processing -- our argument parser found all of the parameters it needed
  parser.parse(args, DeepConfig()) map {
    config: DeepConfig => {
      alignedReads(config)
    }
  } getOrElse {
    println("Unable to parse the command line arguments you passed in, please check that your parameters are correct")
  }


  /**
   * Process the merged and paired reads into the output stats file
   * @param config our config object
   *
   */
  def alignedReads(config: DeepConfig): Unit = {
    val cutsSiteObj = CutSites.fromFile(config.cutSites, 3)

    val outputStatsFile = new StatsOutput(config.outputStats,cutsSiteObj.size)

    val primers = Source.fromFile(config.primersEachEnd.getAbsolutePath).getLines().map { line => line }.toList
    if (primers.length != 2)
      throw new IllegalStateException("You should only provide a primer file with two primers")

    // store our reference and reads
    var refName = ""
    var readName = ""
    var refString = ArrayBuilder.make[String]
    var readString = ArrayBuilder.make[String]
    var inRead = false

    val mergedReadIterator = new ReadPairParser(config.inputMerged)

    println("traversing merged reads...")
    mergedReadIterator.foreach { pair => {
      printMergedRead(cutsSiteObj, outputStatsFile, pair, primers)
    }
    }

    // now process paired reads
    val firstReadIterator = new UnmergedReadParser(config.inputFileUnmerged)

    println("traversing unmerged reads...")
    firstReadIterator.foreach { twoReads => {
      printPairedRead(cutsSiteObj,
        outputStatsFile,
        twoReads,
        primers)
    }
    }
    outputStatsFile.close()
  }

  /**
   * cycle through the merged reads, creating alignments, and dumping information out to the stats file
   * @param cutsSiteObj the cut sites
   * @param outputStatsFile the stats file we're dumping output to
   * @param mergedRead the merged read - aligned read and the matching reference string
   * @param primers the primers to look for on each end of the fragment
   */
  def printMergedRead(cutsSiteObj: CutSites,
                      outputStatsFile: StatsOutput,
                      mergedRead: RefReadPair,
                      primers: List[String]): Unit = {

    val containsFwdPrimer = mergedRead.read.bases.filter(bs => bs != '-').mkString("") contains primers(0)
    val containsRevPrimer = mergedRead.read.bases.filter(bs => bs != '-').mkString("") contains primers(1)

    val baseLen = mergedRead.read.bases.map { case (ch) => if (ch == '-') 0 else 1 }.sum

    val callEvents = AlignmentManager.cutSiteEvent(mergedRead, cutsSiteObj)

    val pass = (containsFwdPrimer && containsRevPrimer && callEvents.matchingRate > .85 && callEvents.matchingBaseCount> 50 && !(callEvents.alignments.mkString("") contains "WT_"))

    outputStatsFile.outputStatEntry(StatsContainer(mergedRead.read.name, pass, containsFwdPrimer, containsRevPrimer,
      false, true, baseLen, -1, 1, 1, callEvents.matchingRate, -1.0, callEvents.matchingBaseCount, -1,
      callEvents.alignments, callEvents.basesOverTargets, None, None, Some(mergedRead.read.bases), None, None, Some(mergedRead.reference.bases)))

  }

  /**
   *
   * @param cutsSiteObj cut sites in our target design
   * @param outputStatsFile the stats output object
   * @param readPairs the read pairs
   * @param primers the primers we expect on each side of the read
   */
  def printPairedRead(cutsSiteObj: CutSites,
                      outputStatsFile: StatsOutput,
                      readPairs: ReadPair,
                      primers: List[String]): Unit = {

    val containsFwdPrimer = readPairs.pair1.read.bases.filter(bs => bs != '-').mkString("") contains primers(0)
    val containsRevPrimer = readPairs.pair2.read.bases.filter(bs => bs != '-').mkString("") contains primers(1)

    val base1Len = readPairs.pair1.read.bases.map { case (ch) => if (ch == '-') 0 else 1 }.sum
    val base2Len = readPairs.pair2.read.bases.map { case (ch) => if (ch == '-') 0 else 1 }.sum

    val callEvents = AlignmentManager.cutSiteEventsPair(readPairs.pair1, readPairs.pair2, cutsSiteObj)

    val pass = containsFwdPrimer && containsRevPrimer &&
      callEvents.matchingRate1 > .85 && callEvents.matchingRate2 > .85 &&
      callEvents.matchingBaseCount1 > 25 && callEvents.matchingBaseCount2 > 25 &&
      !(callEvents.alignments.mkString("") contains "WT_")

    outputStatsFile.outputStatEntry(StatsContainer(readPairs.pair1.read.name, pass, containsFwdPrimer, containsRevPrimer,
      false, false, base1Len, base2Len, 1, 1, callEvents.matchingRate1, callEvents.matchingRate2, callEvents.matchingBaseCount1, callEvents.matchingBaseCount2,
      callEvents.alignments, callEvents.basesOverTargets, Some(readPairs.pair1.read.bases), Some(readPairs.pair2.read.bases), None, Some(readPairs.pair1.reference.bases), Some(readPairs.pair2.reference.bases), None))

  }
}
