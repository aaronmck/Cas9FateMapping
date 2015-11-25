package main.scala

import _root_.utils.{CutSites}
import utils.Utils

import scala.io._
import java.io._
import scala.collection.mutable._
import scala.main.{ReverseReadOrientation, ForwardReadOrientation, SequencingRead}
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
case class Config(inputFileReads1: File = new File(Main.NOTAREALFILENAME),
                  inputFileReads2: File = new File(Main.NOTAREALFILENAME),
                  outputFastq1: File = new File(Main.NOTAREALFILENAME),
                  outputFastq2: File = new File(Main.NOTAREALFILENAME),
                  outputStats: File = new File(Main.NOTAREALFILENAME),
                  outputAlignments: File = new File(Main.NOTAREALFILENAME),
                  cutSites: File = new File(Main.NOTAREALFILENAME),
                  umiLength: Int = 1,
                  umiStartPos: Int = 0,
                  primersEachEnd: File = new File(Main.NOTAREALFILENAME),
                  reference: File = new File(Main.NOTAREALFILENAME),
                  samplename: String = "TEST",
                  minimumUMIReads: Int = 5,
                  minimumSurvivingUMIReads: Int = 2)



object Main extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[Config]("UMIMerge") {
    head("UMIMerge", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("inputFileReads1") required() valueName ("<file>") action { (x, c) => c.copy(inputFileReads1 = x) } text ("first read file ")
    opt[File]("inputFileReads2") required() valueName ("<file>") action { (x, c) => c.copy(inputFileReads2 = x) } text ("second reads file")
    opt[File]("outputFastq1") required() valueName ("<file>") action { (x, c) => c.copy(outputFastq1 = x) } text ("the output stats file")
    opt[File]("outputFastq2") required() valueName ("<file>") action { (x, c) => c.copy(outputFastq2 = x) } text ("the output stats file")
    opt[File]("outputStats") required() valueName ("<file>") action { (x, c) => c.copy(outputStats = x) } text ("the output stats file")
    opt[File]("cutSites") required() valueName ("<file>") action { (x, c) => c.copy(cutSites = x) } text ("the location of the cutsites")
    opt[File]("primersEachEnd") required() valueName ("<file>") action { (x, c) => c.copy(primersEachEnd = x) } text ("the file containing the amplicon primers requred to be present, one per line, two lines total")
    opt[Int]("minimumUMIReads") action { (x, c) => c.copy(minimumUMIReads = x) } text ("the minimum number of reads that each UMI should have to be considered signal and not noise")
    opt[Int]("minimumSurvivingUMIReads") action { (x, c) => c.copy(minimumSurvivingUMIReads = x) } text ("the minimum number of reads that each UMI should have post filtering")

    opt[File]("reference") required() action { (x, c) => c.copy(reference = x) } text ("the reference (as a fasta)")
    opt[Int]("umiStart") required() action { (x, c) => c.copy(umiStartPos = x) } text ("the start position, zero based, of our UMIs")
    opt[Int]("umiLength") required() action { (x, c) => c.copy(umiLength = x) } text ("the length of our UMIs")
    opt[String]("samplename") required() action { (x, c) => c.copy(samplename = x) } text ("the sample name of this run")


    // some general command-line setup stuff
    note("processes reads with UMIs into merged reads\n")
    help("help") text ("prints the usage information you see here")
  }

  // *********************************** Run *******************************************************
  // run the actual read processing -- our argument parser found all of the parameters it needed
  parser.parse(args, Config()) map {
    config: Config => {

      // our output files
      val outputFastq1File = new PrintWriter(config.outputFastq1)
      val outputFastq2File = new PrintWriter(config.outputFastq2)
      val outputStatsFile =  new PrintWriter(config.outputStats)

      // get the reference as a string
      var referenceString = ""
      Source.fromFile(config.reference).getLines().foreach { line => if (!line.startsWith(">")) referenceString += line }

      // setup clustered input of the fastq files
      // ------------------------------------------------------------------------------------------
      val forwardReads = Source.fromInputStream(Utils.gis(config.inputFileReads1.getAbsolutePath)).getLines().grouped(4)
      val reverseReads = Source.fromInputStream(Utils.gis(config.inputFileReads2.getAbsolutePath)).getLines().grouped(4)

      val primers = Source.fromFile(config.primersEachEnd.getAbsolutePath).getLines().map { line => line }.toList
      if (primers.length != 2)
        throw new IllegalStateException("You should only provide a primer file with two primers")

      // our containers for forward and reverse reads
      val umiReadsFWD =     new HashMap[String, ArrayBuffer[SequencingRead]]()
      val umiReadsRVS =     new HashMap[String, ArrayBuffer[SequencingRead]]()


      // --------------------------------------------------------------------------------
      // process the reads into bins of UMIs, keep fwd/rev reads together
      // --------------------------------------------------------------------------------
      print("Reading in sequences and parsing out UMIs (one dot per 10K reads):")
      var readsProcessed = 0
      forwardReads foreach { fGroup => {
        val rGroup = reverseReads.next()

        val umi =              fGroup(1).slice(0, config.umiLength)
        val readNoUMI =        fGroup(1).slice(config.umiLength, fGroup(1).length)
        val qualNoUMI =        fGroup(3).slice(config.umiLength, fGroup(3).length)

        val readBuilderF =     umiReadsFWD.getOrElse(umi, new ArrayBuffer[SequencingRead]())
        val readBuilderR =     umiReadsRVS.getOrElse(umi, new ArrayBuffer[SequencingRead]())

        readBuilderF +=        SequencingRead(fGroup(0),readNoUMI, qualNoUMI, ForwardReadOrientation, umi).qualityThresholdRead(3,30)
        umiReadsFWD(umi) =     readBuilderF

        readBuilderR +=        SequencingRead(rGroup(0), rGroup(1), rGroup(3), ReverseReadOrientation, umi).qualityThresholdRead(3,30)
        umiReadsRVS(umi) =     readBuilderR

        readsProcessed += 1
        if (readsProcessed % 10000 == 0)
          print(".")
      }
      }

      // --------------------------------------------------------------------------------
      // filter out UMIs that are really really abnormal in terms of size -- this will save
      // us a huge amount of computation later
      // --------------------------------------------------------------------------------
      val readCountsPerUmiFWD = ArrayBuilder.make[Double]()
      val readCountsPerUmiREV = ArrayBuilder.make[Double]()

      umiReadsFWD.foreach{case(umi,reads) => readCountsPerUmiFWD += reads.toArray.size}
      umiReadsRVS.foreach{case(umi,reads) => readCountsPerUmiREV += reads.toArray.size}

      val fwdCounts = readCountsPerUmiFWD.result()
      val revCounts = readCountsPerUmiREV.result()

      val fwdMean = fwdCounts.map{vl => vl.toDouble / fwdCounts.length.toDouble}.sum
      val revMean = revCounts.map{vl => vl.toDouble / fwdCounts.length.toDouble}.sum

      val fwdDevs = fwdCounts.map(score => (score - fwdMean) * (score - fwdMean))
      val revDevs = revCounts.map(score => (score - revMean) * (score - revMean))

      val stddevF = Math.sqrt(fwdDevs.sum / fwdCounts.size)
      val stddevR = Math.sqrt(revDevs.sum / revCounts.size)

      val stdDevLimit = 5.0

      val umisToRemove = new HashMap[String,Boolean]()
      umiReadsFWD.foreach{case(umi,reads) => {
        val sz = reads.toArray.size
        if (sz > (fwdMean + (stdDevLimit*stddevF)) || sz < (fwdMean - (2*stddevF)))
          umisToRemove(umi) = true
      }}
      umiReadsRVS.foreach{case(umi,reads) => {
        val sz = reads.toArray.size
        if (sz > (revMean + (stdDevLimit*stddevR)) || sz < (revMean - (2*stddevR)))
          umisToRemove(umi) = true
      }}

      println("Forward mean and std = (" + fwdMean + "," + stddevF + ")")
      println("Forward mean and std = (" + revMean + "," + stddevR + ")")

      umisToRemove.foreach{case(umi,status) => {
        println("removing UMI " + umi + " with forward size " + umiReadsFWD(umi).toArray.size + " and reverse size " + umiReadsRVS(umi).toArray.size)
        umiReadsFWD.remove(umi)
        umiReadsRVS.remove(umi)
      }}

      //outputStats.write("UMI\tused\tffail.reason\tfwd.reads\tprimered.fwd.reads\tfiltered.fwd.reads\trev.reads\tprimered.rev.reads\tfiltered.rev.reads\t")
      //outputStats.write("fwd.error.rate\tfwd.read\trev.error.rate\trev.read\n")
      outputStatsFile.write("UMI\tkeptPCT\tfail.reason\tinitF\tfilteredF\tfinalF\tinitR\tfilteredR\tfinalR\treadF\treadR\tfwdBasesMatching\trevBasesMatching\t")

      var passingUMI = 0
      var totalWithUMI = 0

      println("\nProcessing " + umiReadsFWD.size + " forward read UMIs and " + umiReadsRVS.size + " reverse read UMIs")
      var index = 0

      val cutsSiteObj = CutSites.fromFile(config.cutSites, 3)
      outputStatsFile.write(cutsSiteObj.windows.zipWithIndex.map{case(wds,index) => "target" + (index+1)}.mkString("\t") + "\n")
      // --------------------------------------------------------------------------------
      // for each UMI -- process its individual reads
      // --------------------------------------------------------------------------------
      umiReadsFWD.foreach { case (umi, reads) => {
        val reverseReads = umiReadsRVS(umi).toArray
        val forwardReads = reads.toArray

        if (forwardReads.length > config.minimumUMIReads && reverseReads.length > config.minimumUMIReads) {
          val res = OutputManager.mergeTogether(umi,
            forwardReads,
            reverseReads,
            cutsSiteObj,
            outputFastq1File,
            outputFastq2File,
            outputStatsFile,
            referenceString,
            config.reference,
            primers,
            config.samplename,
            config.minimumSurvivingUMIReads,
            cutsSiteObj)

          passingUMI += res
        }
        index += 1
        if (index % 100 == 0) {
          println("INFO: Processed " + index + " umis so far")
        }
      }
      }

      println("total UMIs passed:\t" + passingUMI + "\ntotal UMIs considered:\t" + umiReadsFWD.size)

      outputFastq1File.close()
      outputFastq2File.close()
      outputStatsFile.close()
    }
  } getOrElse {
    println("Unable to parse the command line arguments you passed in, please check that your parameters are correct")
  }
}
