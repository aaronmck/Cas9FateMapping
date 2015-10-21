/** *********************************************
  * Requires Java 1.7. On the UW machines this means
  * running: "module load java/7u17", or the appropriate
  * module if there's a more recent version, to load java 1.7 into
  * your environment.  The best idea is to place this into your
  * bash/shell profile
  *
  * The only critical input file to provide is the tear sheet. Here's the specification:
  *
  * a *__TAB__* separated file with the following columns (the highlighting to be clear about it being a tab)
  * column names are required and must match, e.g. library_name must be the first column.  Your header 
  * should then have the following column names:
  *
  * library_type - the type of the library, for instance DNA, RNA, NEXTERA, etc
  * sample - the sample name
  * project - the overarching project for this sequencing
  * fastq1 - the fastq 1 file - the first reads
  * fastq2 - the fastq 2 file - the second reads
  * fastqBarcode 1 - the fastq barcode file, first end
  * fastqBarcode 2 - the fastq barcode file, first end
  * output_dir - what directory to place the output into
  * sequencing_date - in the format day_month_fullyear
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
  * @date May 2nd, 2014
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
import java.io.File

/**
 * Quality control sequencing data and optionally align the data to the genome
 **/
class DNAQC extends QScript {
  qscript =>

  /** **************************************************************************
    * Required Parameters
    * ************************************************************************** */
  @Input(doc = "Tab-delimited (seperated with commas)", fullName = "input", shortName = "i", required = true)
  var input: File = _

  /** **************************************************************************
    * Optional Parameters -- control parameters for the script (alignment, etc)
    * ************************************************************************** */

  @Argument(doc = "read length (how long is each read in a pair (300 cycle kit -> 150bp reads = 150)", fullName = "readLength", shortName = "readLength", required = false)
  var readLength: Int = 250

  /** **************************************************************************
    * Path parameters -- where to find tools
    * ************************************************************************** */

  @Input(doc = "The path to the Trimmomatic adapters file", fullName = "trimadapters", shortName = "tad", required = false)
  var adaptersFile: File = new File("/net/gs/vol1/home/aaronmck/tools/trimming/Trimmomatic-0.32/adapters/TruSeq3-PE.fa")

  @Input(doc = "The path to the picard libraries", fullName = "picard", shortName = "pcd", required = false)
  var picardPath: File = "/net/gs/vol1/home/aaronmck/tools/picard/picard-tools-1.95/"

  @Input(doc = "The path to the binary of fastqc", fullName = "fastqc", shortName = "fqc", required = false)
  var fastqcPath: File = "/net/gs/vol1/home/aaronmck/tools/fastqc/FastQC/fastqc"

  @Input(doc = "The path to the jar file for trimmomatic", fullName = "trim", shortName = "tm", required = false)
  var trimmomaticPath: File = "/net/gs/vol1/home/aaronmck/tools/trimming/Trimmomatic-0.32/trimmomatic-0.32.jar"

  @Input(doc = "The path to the binary of samtools", fullName = "samtools", shortName = "smt", required = false)
  var samtoolsPath: File = "/net/gs/vol1/home/aaronmck/tools/samtools/samtools-0.1.19/samtools"

  @Input(doc = "The path to the binary of bwa", fullName = "bwa", shortName = "bwa", required = false)
  var bwaPath: File = "/net/gs/vol1/home/aaronmck/tools/bwa/bwa-0.7.9a/bwa"

  @Input(doc = "The path to the barcode splitter", fullName = "maulpath", shortName = "mlp", required = false)
  var maulPath: File = "/net/gs/vol1/home/aaronmck/source/Maul/target/scala-2.10/maul_2.10-1.0-one-jar.jar"

  @Input(doc = "The path to the seqprep tool", fullName = "seqprep", shortName = "sprep", required = false)
  var seqPrepPath: File = "/net/gs/vol1/home/aaronmck/tools/bin/SeqPrep"

  @Input(doc = "the script to analyze the edits made to a crispr target", fullName = "edits", shortName = "edits", required = false)
  var crisprPath: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_05_10_Rerun_HEK_analysis/bin/analyze_diversity_of_edits.scala"

  @Input(doc = "the script to analyze the per site edits", fullName = "posScript", shortName = "posScript", required = false)
  var positionAnalysisPath: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_05_10_Rerun_HEK_analysis/bin/sam_to_positions.scala"

  @Input(doc = "the script to analyze the per site edits", fullName = "writeweb", shortName = "wwwrite", required = false)
  var writeWebFiles: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_05_10_Rerun_HEK_analysis/bin/write_web_files.scala"

  @Input(doc = "the location of the smith waterman program", fullName = "swp", shortName = "swp", required = false)
  var smithWaterman: File = "/net/gs/vol1/home/aaronmck/tools/CRISPR/process_ssw_alignment.scala"

  @Input(doc = "the cut site analysis script", fullName = "csa", shortName = "csa", required = false)
  var cutSiteScript: File = "/net/gs/vol1/home/aaronmck/tools/cut_site/cut_site_analysis.Rmd"

  @Input(doc = "where to put the web files", fullName = "web", shortName = "www", required = false)
  var webSite: File = "/net/shendure/vol2/www/content/members/aaron/crispr/"

  @Input(doc = "where to find the diagram files", fullName = "diag", shortName = "diag", required = false)
  var diagramLoc: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_07_08_HEK4_exp3_redo/data/"

  /** **************************************************************************
    * Global Variables
    * ************************************************************************** */

  // Gracefully hide Queue's output -- this is very important otherwise there's Queue output everywhere and it's a mess
  val queueLogDir: String = ".qlog/"

  /** **************************************************************************
    * Main script
    * ************************************************************************** */
  def script() {
    // read in the tear sheet and process each sample
    parseTearSheet(input).foreach(sampleObj => {

      // the sample name we'll use for files -- the sample name plus the library name
      val sampleTag = sampleObj.sample

      // check that our basic output dir can be made
      dirOrCreateOrFail(sampleObj.outputDir, "base output directory")

      // are we running in double or single ended mode?
      if ((sampleObj.fastq1.toUpperCase != "NA" && !sampleObj.fastq1.exists()) || (sampleObj.fastq2.toUpperCase != "NA" && !sampleObj.fastq2.exists())) {
        throw new IllegalStateException("Unable to find one of the read fastqs from the sample " + sampleObj.sample + " please check the tear sheet: fq1 = " + sampleObj.fastq1 + " fq2 = " + sampleObj.fastq2)
      }

      // are we using zero, one, or two barcodes?
      val oneBarcode = sampleObj.fastqBarcode1.exists() && !sampleObj.fastqBarcode1.exists()
      val dualBarcode = sampleObj.fastqBarcode1.exists() && sampleObj.fastqBarcode2.exists()
      val pairedEnd = sampleObj.fastq2.exists()

      // **************************** Setup files and directories ******************************************
      // create the per-sample output directories, verifying that they either exist already or can be made
      val sampleOutput = dirOrCreateOrFail(new File(sampleObj.outputDir + File.separator + sampleObj.sample), "sample output directory")
      val initialFastQCDir = dirOrCreateOrFail(new File(sampleOutput + File.separator + "initial_fastQC"), "initial fastqc directory")
      val initialFastQMerged = dirOrCreateOrFail(new File(sampleOutput + File.separator + "initial_fastQC_Merged"), "initial fastqc directory")
      val postCleanFastQCDir = dirOrCreateOrFail(new File(sampleOutput + File.separator + "post_clean_fastQC"), "final fastqc directory")

      // our barcode split files
      var barcodeSplit1 = new File(sampleOutput + File.separator + sampleTag + ".barcodeSplit.fastq1.fq.gz")
      var barcodeSplit2 = new File(sampleOutput + File.separator + sampleTag + ".barcodeSplit.fastq2.fq.gz")

      val sam = new File(sampleOutput + File.separator + sampleTag + ".sam")
      val bam = new File(sampleOutput + File.separator + sampleTag + ".bam")

      val mergedReads = new File(sampleOutput + File.separator + sampleTag + ".merged.fq.gz")
      val samMerged = new File(sampleOutput + File.separator + sampleTag + ".merged.sam")
      val bamMerged = new File(sampleOutput + File.separator + sampleTag + ".merged.bam")
      
      val reheaderBam = new File(sampleOutput + File.separator + sampleTag + ".reheadered.bam") // used only when we don't acutally align
      val reheaderAddReadGroupBam = new File(sampleOutput + File.separator + sampleTag + ".readgroup.reheadered.bam") // used only when we don't acutally align
      val sortedBam = new File(sampleOutput + File.separator + sampleTag + ".sorted.bam")
      val readGroupedBam = new File(sampleOutput + File.separator + sampleTag + ".sorted.withRG.bam")
      val sortedDeDuppedBam = new File(sampleOutput + File.separator + sampleTag + ".sorted.withRG.dedupped.bam")
      val sortedDeDuppedMetrics = new File(sampleOutput + File.separator + sampleTag + ".sorted.withRG.dedupped.metrics")

      val reheaderBamMerged = new File(sampleOutput + File.separator + sampleTag + ".reheadered.merged.bam") // used only when we don't acutally align
      val reheaderAddReadGroupBamMerged = new File(sampleOutput + File.separator + sampleTag + ".readgroup.reheadered.merged.bam") // used only when we don't acutally align
      val sortedBamMerged = new File(sampleOutput + File.separator + sampleTag + ".sorted.merged.bam")
      val readGroupedBamMerged = new File(sampleOutput + File.separator + sampleTag + ".sorted.withRG.merged.bam")
      val sortedDeDuppedBamMerged = new File(sampleOutput + File.separator + sampleTag + ".sorted.withRG.dedupped.merged.bam")
      val sortedDeDuppedMetricsMerged = new File(sampleOutput + File.separator + sampleTag + ".sorted.withRG.dedupped.merged.metrics")

      // output files for Trimmomatic
      val postProcessedFastq1 = new File(sampleOutput + File.separator + sampleTag + ".cleaned.trim.fastq1.fq.gz")
      val postProcessedFastq2 = new File(sampleOutput + File.separator + sampleTag + ".cleaned.trim.fastq2.fq.gz")
      val postProcessedFastq1seqPrep = new File(sampleOutput + File.separator + sampleTag + ".cleaned.trim.sp.fastq1.fq.gz")
      val postProcessedFastq2seqPrep = new File(sampleOutput + File.separator + sampleTag + ".cleaned.trim.sp.fastq2.fq.gz")

      val postProcessedFastq1Unpaired = new File(sampleOutput + File.separator + sampleTag + ".fastq1.unpaired.fq.gz")
      val postProcessedFastq2Unpaired = new File(sampleOutput + File.separator + sampleTag + ".fastq2.unpaired.fq.gz")
      val spliceAlignments = new File(sampleOutput + File.separator + "SJ.out.tab")

      // ************************************** split the input files by barcodes, either one or two **************************************
      var barcodeFiles = if (dualBarcode) List[File](sampleObj.fastqBarcode1,sampleObj.fastqBarcode2) else List[File](sampleObj.fastqBarcode1)

      var inputFiles = if (pairedEnd) List[File](sampleObj.fastq1,sampleObj.fastq2) else List[File](sampleObj.fastq1)
      var processedFastqs = if (pairedEnd) List[File](barcodeSplit1,barcodeSplit2) else List[File](barcodeSplit1)
      val barcodeConfusion = new File(sampleOutput + File.separator + sampleTag + ".barcodeConfusion")
      val barcodeStats = new File(sampleOutput + File.separator + sampleTag + ".barcodeStats")
      val overlapFile = new File(sampleOutput + File.separator + sampleTag + ".readOverlap")

      // split by barcodes if required
      if (dualBarcode)
        add(Maul(inputFiles, barcodeFiles, List[String](sampleObj.barcode1,sampleObj.barcode2), processedFastqs, barcodeStats, barcodeConfusion,overlapFile))
      else 
        add(Maul(inputFiles, barcodeFiles, List[String](sampleObj.barcode1), processedFastqs, barcodeStats, barcodeConfusion,overlapFile))

      inputFiles = if (pairedEnd) List[File](barcodeSplit1,barcodeSplit2) else List[File](barcodeSplit1)
      
      // run an Fastqc so we know what we started with 
      // ************************************************************************************
      
      add(Fastqc(inputFiles, initialFastQCDir))

      // QC and trim the reads down, removing bad sections and remaining adapters that might have snuck through
      // ************************************************************************************
      processedFastqs = if (pairedEnd) List[File](postProcessedFastq1,postProcessedFastq2) else List[File](postProcessedFastq1)
      val unpairedOutputFiles = if (pairedEnd) List[File](postProcessedFastq1Unpaired,postProcessedFastq2Unpaired) else List[File](postProcessedFastq1Unpaired)

      // trim off the ends of the reads,and set these as inputs to our seqprep
      add(Trimmomatic(inputFiles,adaptersFile,processedFastqs,unpairedOutputFiles))
      inputFiles = processedFastqs
      processedFastqs = if (pairedEnd) List[File](postProcessedFastq1seqPrep,postProcessedFastq2seqPrep) else List[File](postProcessedFastq1seqPrep)

      add(BWAAlignPairedEndWithMEM(sampleObj.reference, processedFastqs, sam))
      add(SamToBam(sam, bam))
    
      // follow through with the rest of the pipeline for all bam files
      add(Index(bam, bam.replace(".bam", ".bai"))) // index, this greatly improves the speed of the resort call
      add(ResortSam(bam, sortedBam)) // re-sort the bam file
      add(Index(sortedBam, sortedBam.replace(".bam", ".bai")))
      
      // now analyze the edits
      val outputAlleles = new File(sampleOutput + File.separator + sampleTag + ".alleles")
      val outputControl = new File(sampleOutput + File.separator + sampleTag + ".control")
      val outputMatchedReads = new File(sampleOutput + File.separator + sampleTag + ".merged.matched.primers.bam")
      val outputMatchedReadsRG = new File(sampleOutput + File.separator + sampleTag + ".merged.matched.primers.wrg.bam")
      val cutPosOutput = new File(sampleOutput + File.separator + sampleTag + ".cutPosOutput")
      val cutDataOutput = new File(sampleOutput + File.separator + sampleTag + ".cutData")
      val outputCase = new File(sampleOutput + File.separator + sampleTag + ".case")
      val outputAnalysis = new File(sampleOutput + File.separator + sampleTag + ".analysis.html")
      val outputAnalysisFigures = new File(sampleOutput + File.separator + "figures")

      // get site specific information, regardless of the case / control status
      add(SamToPos(cutPosOutput, sampleObj.reference, bamMerged))

      // analze the edits if we're not a control
      if (sampleObj.caseControl.toUpperCase != "TRUE") {
        val otherSampleOutput = new File(sampleObj.outputDir + File.separator + sampleObj.caseControl + File.separator)
        val otherSortedDeDuppedBamMerged = new File(otherSampleOutput + File.separator + sampleObj.caseControl + ".sorted.merged.bam")
        val otherCutPos = new File(otherSampleOutput + File.separator + sampleObj.caseControl + ".cutPosOutput")

        add(AnalyzeEdits(outputAlleles, outputCase, outputControl, sampleObj.reference, sortedBamMerged, otherSortedDeDuppedBamMerged, outputMatchedReads, cutDataOutput))

        add(AddReadGroup(outputMatchedReads,
          outputMatchedReadsRG,
          ReadGroup(sampleObj.sample,sampleObj.caseControl,"ILLUMINA",sampleObj.sample+sampleObj.caseControl,sampleObj.sample,"UW","STUFF")))

        add(Index(outputMatchedReadsRG, outputMatchedReadsRG.replace(".bam", ".bai")))
        add(PostProcess(cutDataOutput, outputAnalysis, new File(sampleOutput + File.separator)))

        // setup the web configuration file
        val webLoc = webSite + File.separator + sampleObj.sample + File.separator
        val webSiteFile = webLoc + "configuration.js"

        println("Running setup " + webSiteFile)
        add(SetupWebAnalysis(sampleObj.reference, sampleObj.reference + ".cutSites", sortedBamMerged, otherSortedDeDuppedBamMerged, cutPosOutput, otherCutPos, webLoc, outputAnalysis, outputAnalysisFigures, webSiteFile, diagramLoc))
      }

    })
  }

  /** **************************************************************************
    * Helper classes and methods
    * ************************************************************************** */

  // if a directory doesn't exist, create it. Otherwise create it. if that fails, exception out.
  // return the directory as a file
  def dirOrCreateOrFail(dir: File, contentsDescription: String): File = {
    if (!dir.exists())
      if (!dir.mkdirs()) // mkdirs tries to make all the parent directories as well, if it can
        throw new IllegalArgumentException("Unable to find or create " + contentsDescription + " directory: " + dir)
      else
        println("created directory : " + dir.getAbsolutePath)
    dir
  }

  // function to check if the bwa process is likely to fail (BWA doesn't always return a non-zero exit code)
  // so try to cut down on the mistakes by verifying a few things before we begin.  Either return the BWA parameters
  // or crap out with an exception if things are setup incorrectly
  def verifyBWAParameters(bwaSetup: BWAAlignPairedEndWithMEM): BWAAlignPairedEndWithMEM = {
    if (!(new File(bwaSetup.reference + ".sa")).exists)
      throw new IllegalStateException("Unable to find the .sa file for the reference, index your reference using BWA!")
    if (!(new File(bwaSetup.reference + ".bwt")).exists)
      throw new IllegalStateException("Unable to find the .bwt file for the reference, index your reference using BWA!")
    bwaSetup
  }

  // The storage container for the data we've read from the input tear sheet
  case class SourceEntry(
    val sample: String,
    val caseControl: String,
    val reference: File,
    val outputDir: File,
    val fastq1: File,
    val fastq2: File,
    val fastqBarcode1: File,
    val fastqBarcode2: File,
    val barcode1: String,
    val barcode2: String)

  // a read group, as defined in a sam/bam file
  case class ReadGroup(val id: String, // readgroup id
    val lb: String, // the library name
    val pl: String, // platform name: almost always ILLUMINA
    val pu: String, // a platform unit, should be unique to this exact sample+run combination
    val sm: String, // the sample name
    val cn: String, // sequencing center
    val ds: String)

  // a description of the sequencing data

  /**
    * reads in the sample tear sheet, a tab separated file with the columns listed at the top of the file
   */
  def parseTearSheet(input: File): Array[SourceEntry] = {
    val inputLines = Source.fromFile(input).getLines

    // check that the header contains the correct information
    if (inputLines.next().stripLineEnd != "sample\tcontrol\treference\toutput.dir\tfastq1\tfastq2\tbarcode.fastq1\tbarcode.fastq2\tbarcode1\tbarcode2")
      throw new IllegalArgumentException("Your header doesn't seem like a correctly formatted  tear sheet!")

    return inputLines.map(line => {
      val tokens = line.split( """\s+""") // Use a regex to avoid empty tokens
      try {
        (new SourceEntry(
          tokens(0), // sample
          tokens(1), // case control status: true if control, else a name of a sample for the control
          new File(tokens(2)), // reference
          new File(tokens(3)), // output
          new File(tokens(4)), // fastq1
          new File(tokens(5)), // fastq2
          new File(tokens(6)), // barcode fastq1
          new File(tokens(7)), // barcode fastq2
          tokens(8), // barcode1
          tokens(9) // barcode2
        ))
      }
      catch {
        case e: java.lang.ArrayIndexOutOfBoundsException => {
          println("\n\n****\nunable to find all the needed columns from the input file line: " + line + "\n***\n"); throw e
        }
      }
    }).toArray
  }

   /** Turn a source entry into a read group.  The mapping is rather simple */
  def sourceToReadGroup(source: SourceEntry): ReadGroup = ReadGroup(id = source.sample, lb = source.sample, pl = "ILLUMINA", pu = source.sample, sm = source.sample, cn = "UW", ds = source.sample)

  /** **************************************************************************
    * traits that get tacked onto runnable objects
    * ************************************************************************** */

  // General arguments to non-GATK tools
  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 6 // set this a bit high, there's a weird java / SGE interactions that require higher mem reservations
    this.residentRequest = 6
    this.residentLimit = 6
    this.isIntermediate = false // by default delete the intermediate files, if you want to keep output from a task set this to false in the case class
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
    this.maxRecordsInRam = 100000
  }

  /** **************************************************************************
    * Classes (non-GATK programs)
    *
    * a note here: normally we could use the built in Picard (and GATK) functions,
    * but for some reason our SGE interacts very poorly with Java, and we need to
    * build in an even higher 'buffer' of memory (the difference between the
    * requested memory to SGE and Java's Xmx parameter)
    * ************************************************************************** */

  // ********************************************************************************************************
  case class AddReadGroup(inBam: File, outBam: File, readGroup: ReadGroup) extends ExternalCommonArgs {
    @Input(doc = "the input bam file") var input = inBam
    @Output(doc = "the output bam file") var output = outBam

    var cmd = "java -Xmx4g -jar " + picardPath + "/AddOrReplaceReadGroups.jar INPUT=" + input + " OUTPUT=" + output
    cmd += " RGID=" + readGroup.id + " RGLB=" + readGroup.lb + " RGPL=ILLUMINA RGPU=" + readGroup.pl
    cmd += " RGSM=" + readGroup.sm + " RGCN=" + readGroup.cn + " RGDS=" + readGroup.ds

    def commandLine = cmd

    this.analysisName = queueLogDir + outBam + ".rg"
    this.jobName = queueLogDir + outBam + ".rg"
    this.memoryLimit = 8
  }


  // set the max file handles to 1000, most UW machines are set to 1024 for some unknown (1970's inspired) reason
  // ********************************************************************************************************
  case class Dedup(inBam: File, outBam: File, metricsFile: File, maxFileHandles: Int = 1000) extends ExternalCommonArgs {
    @Input(doc = "the input bam file") var input = inBam
    @Output(doc = "the output bam, with duplicates marked") var output = outBam
    @Output(doc = "the metrics file, with information on the duplicate counts") var metrics = metricsFile
    @Argument(doc = "the max number of file handles we can have open at any one time") var maxFLH = maxFileHandles

    def commandLine = "java -Xmx4g -jar " + picardPath + "/MarkDuplicates.jar INPUT=" + input + " OUTPUT=" + output + " METRICS_FILE=" + metrics + " MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=" + maxFLH

    this.analysisName = queueLogDir + outBam + ".dedup"
    this.jobName = queueLogDir + outBam + ".dedup"
    this.memoryLimit = 8
  }

  // ********************************************************************************************************
  case class ResortSam(inSamBam: File, outBam: File, sortOrderP: SortOrder = SortOrder.coordinate) extends ExternalCommonArgs {
    @Input(doc = "the input bam or sam file") var input = inSamBam
    @Output(doc = "the output bam, with duplicates marked") var output = outBam
    @Argument(doc = "the sort order file") var sortP = sortOrderP

    def commandLine = "java -Xmx4g -jar " + picardPath + "/SortSam.jar INPUT=" + input + " OUTPUT=" + output + " SORT_ORDER=" + sortP

    this.analysisName = queueLogDir + outBam + ".sortSam"
    this.jobName = queueLogDir + outBam + ".sortSam"
    this.memoryLimit = 8
  }

  // our min indentical bases is set to 8, having problems with memory so I turned this up from the default 5
  // ********************************************************************************************************
  case class EstimateCompl(inputBam: File, outputComplexity: File, minIdentBases: Int = 8) extends ExternalCommonArgs {
    @Input(doc = "the input bam file") var input = inputBam
    @Output(doc = "the complexity file") var complex = outputComplexity
    @Argument(doc = "the number of bases that have to be identical at the beginning of a read to compare") var minIB = minIdentBases

    def commandLine = "java -Xmx6 -jar " + picardPath + "/EstimateLibraryComplexity.jar INPUT=" + input + " OUTPUT=" + complex + " READ_NAME_REGEX=null"
    this.isIntermediate = false
  }

  // add a header onto each of our bam files
  // ********************************************************************************************************
  case class ReplaceSamHeader(inputBam: File, outputBam: File, header: File) extends ExternalCommonArgs {
    @Input(doc = "the input bam file") var inBam = inputBam
    @Output(doc = "the complexity file") var outBam = outputBam
    @Argument(doc = "the header file") var hd = header

    def commandLine = "java -Xmx6 -jar " + picardPath + "/ReplaceSamHeader.jar INPUT=" + inBam + " OUTPUT=" + outBam + " HEADER=" + header
  }

  // ********************************************************************************************************
  case class Maul(inputFastqs: List[File], barcodeFiles: List[File], barcodes: List[String], outputFastqs: List[File], barcodeStats: File, barcodeConfusion: File, overlap: File) extends ExternalCommonArgs {
    @Input(doc = "the input fastqs") var inputFQs = inputFastqs
    @Input(doc = "the barcode files") var barcodeFls = barcodeFiles
    @Argument(doc = "the barcodes") var barCds = barcodes
    @Output(doc = "the output fastq files") var outputFQs = outputFastqs
    @Argument(doc = "the output barcode counts") var outStats = barcodeStats
    @Argument(doc = "the output barcode confusion file") var outBCC = barcodeConfusion
    @Argument(doc = "the output overlap file") var outOver = overlap

    // the base command to run
    var baseCmd = "java -Xmx4g -jar " + maulPath + " --inFQ1 " + inputFQs(0) + " --outFQ1 " + outputFQs(0)

    // figure out if we're running with one or two barcodes
    if (inputFastqs.length == 2) baseCmd += " --inFQ2 " + inputFQs(1) + " --outFQ2 " + outputFQs(1)

    // if we have an initial barecode
    if (barcodeFls.length > 0) baseCmd += " --inBarcodeFQ1 " + barcodeFls(0) + " --barcodes1 " + barCds(0)
    
    // or a second barcode
    if (barcodeFls.length > 1) baseCmd += " --inBarcodeFQ2 " + barcodeFls(1) + " --barcodes2 " + barCds(1)

    def commandLine = baseCmd + " --barcodeStatsFile " + outStats + " --barcodeStatsFileUnknown " + barcodeConfusion + 
    				  " --overlapFile " + outOver + " --readLength " + readLength

    this.isIntermediate = false
    this.memoryLimit = 8
  }

  // evaluate the quality of the reads from a fastq file
  // ********************************************************************************************************
  case class Fastqc(inputFQs: List[File], outputDir: File, isBam: Boolean = false) extends ExternalCommonArgs {
    @Input(doc = "the input files") var ins = inputFQs
    @Input(doc = "output directory") var outdir = outputDir

    def commandLine = fastqcPath + " -o " + outdir + " -f fastq " + (if (isBam) " -f bam " else "") + ins.mkString(" ")

    this.memoryLimit = 8
    this.analysisName = queueLogDir + ins(0).getName() + ".fastqc" + (if (isBam) "BAM" else "")
    this.jobName = queueLogDir + ins(0).getName() + ".fastqc" + (if (isBam) "BAM " else "")

    // DO NOT REMOVE -- Queue will delete all the output, including the output directory unless you set this.  This means there's some cleanup at the end though
    this.isIntermediate = false
  }

  // ********************************************************************************************************
  case class SamToBam(inSam: File, outBam: File) extends ExternalCommonArgs {
    @Input(doc = "the input sam file") var sam = inSam
    @Output(doc = "output bam file") var bam = outBam

    def commandLine = samtoolsPath + " view -bS " + sam + " > " + bam

    this.analysisName = queueLogDir + outBam + ".SamToBam"
    this.jobName = queueLogDir + outBam + ".SamToBam"
  }

  // ********************************************************************************************************
  case class FlagStat(inSam: File, outStats: File) extends ExternalCommonArgs {
    @Input(doc = "the input sam file") var sam = inSam
    @Output(doc = "output bam file") var stats = outStats

    def commandLine = samtoolsPath + " flagstat " + sam + " > " + outStats

    this.analysisName = queueLogDir + inSam + ".flagstat"
    this.jobName = queueLogDir + inSam + ".flagstat"
    this.isIntermediate = false
  }

  // ********************************************************************************************************
  case class Index(inBam: File, indexFile: File) extends ExternalCommonArgs {
    @Input(doc = "the input sam file") var bam = inBam
    @Output(doc = "output bam file") var ind = indexFile

    def commandLine = samtoolsPath + " index " + bam + " " + ind

    this.analysisName = queueLogDir + inBam + ".index"
    this.jobName = queueLogDir + inBam + ".index"
    this.isIntermediate = false
  }

  // convert a set of fastq files to a sam file
  // ********************************************************************************************************
  case class FastqToSam(inFastqs: List[File], outSam: File, readGroup: ReadGroup) extends ExternalCommonArgs {
    @Input(doc = "the input fastq files") var ins = inFastqs
    @Output(doc = "output sam file") var sam = outSam

    var cmdLineBuild = "java -Xmx4g -jar " + picardPath + "/FastqToSam.jar FASTQ=" + ins(0) + (if (ins.length == 2) {" FASTQ2=" + ins(1)} else "") + " OUTPUT=" + outSam
    cmdLineBuild += " READ_GROUP_NAME=" + readGroup.id + " SAMPLE_NAME=" + readGroup.sm + " LIBRARY_NAME=" + readGroup.lb
    cmdLineBuild += " PLATFORM_UNIT=" + readGroup.pu + " PLATFORM=" + readGroup.pl + " SEQUENCING_CENTER=shendure" + " SORT_ORDER=unsorted"

    def commandLine = cmdLineBuild

    this.analysisName = queueLogDir + sam.getName + ".FastqToSam"
    this.jobName = queueLogDir + sam.getName + ".FastqToSam"
  }

  // BWA mem alignment -- all in one step
  // ********************************************************************************************************
  case class BWAAlignPairedEndWithMEM(reference: File, inFastqs: List[File], outputSam: File, bwaThreads: Int = 1)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the read fastqs") var fqs = inFastqs
    @Argument(doc = "the reference fasta/fa") var ref = reference
    @Output(doc = "the output sam file") var outSam = outputSam

    def commandLine = bwaPath + " mem -t " + bwaThreads + "  " + ref + " " + fqs.mkString(" ") + " > " + outSam

    this.analysisName = queueLogDir + outSam + ".bwa_mem"
    this.jobName = queueLogDir + outSam + ".bwa_mem"
  }

  // Smith Waterman alignment
  // ********************************************************************************************************
  case class SmithWaterman(reference: File, inFastqs: List[File], outputSam: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the read fastqs") var fqs = inFastqs
    @Argument(doc = "the reference fasta/fa") var ref = reference
    @Output(doc = "the output sam file") var outSam = outputSam

    def commandLine = "scala " + smithWaterman + " " + ref + " " + fqs.mkString(" ") + " " + outSam

    this.analysisName = queueLogDir + outSam + ".sw"
    this.jobName = queueLogDir + outSam + ".sw"
  }

  // analyze the edits from the final alignments
  // ********************************************************************************************************
  case class AnalyzeEdits(outputAllele: File, outputCase: File, outputControl: File, reference: File, bamCut: File, bamControl: File, matchedPrimerReads: File, cutDataFile: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Output(doc = "the output alleles") var outa = outputAllele
    @Output(doc = "the output control data") var outcont = outputControl
    @Output(doc = "the output case data") var outcase = outputCase
    @Input(doc = "the reference") var ref = reference
    @Input(doc = "the cut bam file") var bCut = bamCut
    @Input(doc = "the control bam file") var bControl = bamControl
    @Output(doc = "the reads where both ends match the primers") var matchedReads = matchedPrimerReads
    @Output(doc = "the output cut data") var cutData = cutDataFile

    def commandLine = "scala -J-Xmx4g -cp /net/gs/vol1/home/aaronmck/source/sandbox/aaron/lib/target/scala-2.10/sequencing-lib-assembly-1.0.jar " + crisprPath + " " + outa + " " + outcase + " " + outcont + " " + ref + ".primers " + ref + ".cutSites " + bCut + " " + bControl + " " + matchedReads + " " + cutData

    this.analysisName = queueLogDir + outputAllele + ".crisprEdits"
    this.jobName = queueLogDir + outputAllele + ".crisprEdits"
    this.isIntermediate = false
  }

// analyze the edits from the final alignments
  // ********************************************************************************************************
  case class PostProcess(dataFile: File, outputFile: File, workingDir: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the data file") var dtable = dataFile
    @Input(doc = "the data file") var workdir = workingDir
    @Output(doc = "the output file") var outf = outputFile

    //def commandLine = "R -e \"library(knitr);library(knitrBootstrap);library(rmarkdown);table.data='" + dtable + "';opts_knit\\$set(base.dir = '" + workdir + "');render('/net/gs/vol1/home/aaronmck/tools/CRISPR/cut_site_analysis.Rmd','knitrBootstrap::bootstrap_document',output_file='" + outf + "')\""
    def commandLine = "R -e \"library(knitr);table.data='" + dtable + "';opts_knit\\$set(base.dir = '" + workdir + "');knit2html('/net/gs/vol1/home/aaronmck/tools/CRISPR/cut_site_analysis.Rmd',output='" + outf + "')\""

    //def commandLine = "R -e \"library(knitr);table.data='" + dtable + "';opts_knit\\$set(base.dir = '" + workdir + "');render('/net/gs/vol1/home/aaronmck/tools/CRISPR/cut_site_analysis.Rmd','knitrBootstrap::bootstrap_document',output_file='" + outf + "')\""


    this.analysisName = queueLogDir + outputFile + ".output"
    this.jobName = queueLogDir + outputFile + ".output"
    this.isIntermediate = false
  }

  // analyze the edits from the final alignments
  // ********************************************************************************************************
  case class SamToPos(cutPosOutput: File, reference: File, bamFile: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Output(doc = "cut per positions file") var cutPerFile = cutPosOutput
    @Input(doc = "the reference") var ref = reference
    @Input(doc = "the control bam file") var bam = bamFile

    def commandLine = "scala -J-Xmx4g -cp /net/gs/vol1/home/aaronmck/source/sandbox/aaron/lib/target/scala-2.10/sequencing-lib-assembly-1.0.jar " + positionAnalysisPath + " " + reference + " " + bam + " " + cutPerFile + " " + cutPerFile + ".sizes " + ref + ".primers"

    this.analysisName = queueLogDir + cutPerFile + ".crisprsites"
    this.jobName = queueLogDir + cutPerFile + ".crisprsites"
    this.isIntermediate = false
  }

  // gzip a bunch of fastqs into a single gzipped fastq
  // ********************************************************************************************************
  case class ConcatFastqs(inFastqs: List[File], outFastq: File) extends ExternalCommonArgs {
    @Input(doc = "the read fastqs") var fqs = inFastqs
    @Output(doc = "the output compressed read file") var outfq = outFastq

    def commandLine = "gzip -c " + fqs.mkString(" ") + " > " + outfq

    this.analysisName = queueLogDir + outFastq + ".concatFastqs"
    this.jobName = queueLogDir + outFastq + ".concatFastqs"
  }

  // gzip a bunch of fastqs into a single gzipped fastq
  // ********************************************************************************************************
  case class SetupWebAnalysis(reference: File,
    refCuts: File,
    caseBAM: File,
    controlBAM: File,
    caseFile: File,
    controlFile: File,
    webLocation: File,
    htmlFile: File,
    figDir: File,
    webConfigFile: File,
    diagLoc: File
  ) extends ExternalCommonArgs {

    @Input(doc = "the reference file") var ref = reference
    @Input(doc = "the reference cut locations") var refcut = refCuts
    @Input(doc = "the case bam file") var casebam = caseBAM
    @Input(doc = "the control bam file") var controlbam = controlBAM
    @Input(doc = "the case cut file") var casecut = caseFile
    @Input(doc = "the control cut file") var controlcut = controlFile
    @Input(doc = "the html analysis file") var html = htmlFile
    @Input(doc = "the diagram html and js location") var diag = diagLoc
    @Argument(doc = "the figure directory") var figures = figDir
    @Argument(doc = "output location") var outDir = webLocation
    @Output(doc = "output web configuration file") var webConfig = webConfigFile


    def commandLine = "scala -J-Xmx4g -cp /net/gs/vol1/home/aaronmck/source/sandbox/aaron/lib/target/scala-2.10/sequencing-lib-assembly-1.0.jar " + writeWebFiles + " " + ref + " " + refcut + " " + casebam + " " + controlbam + " " + " " + caseFile + " " + controlFile + " " + ref + ".primers " + outDir + " " + htmlFile + " " + figDir + " " + webConfig + " " + diag

    this.analysisName = queueLogDir + caseFile + ".cpFile"
    this.jobName = queueLogDir + caseFile + ".cpFile"
    this.isIntermediate = false
  }


  /**
   * Trimmomatic -- http://www.usadellab.org/cms/?page=trimmomatic
   *
   * trims off adapter sequence and low quality regions in paired end sequencing data
   */
  // ********************************************************************************************************
  case class Trimmomatic(inFastqs: List[File],
                         adapters: File,
                         outs: List[File],
                         outputUnpaired: List[File]
                          ) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input FASTAs") var fqs = inFastqs
    @Argument(doc = "the adapters file") var adp = adapters

    @Output(doc = "output fastas (corrected)") var fqOuts = outs
    @Output(doc = "output fastas for failed reads") var fqOutsUnpaired = outputUnpaired

    // the parameters string to tack onto the end, specifying:
    var appended = "ILLUMINACLIP:" + adapters.getAbsolutePath() + ":2:30:12" // where to find the illumina adapter file
    appended += " LEADING:3" // trim off leading bases with a quality less than 3
    appended += " TRAILING:3" // trim off trailing bases with a quality less than 3
    appended += " SLIDINGWINDOW:5:8" // remove bases with a quality less than 8 in a sliding window of size 5 (so if 5 bases have an mean qual < 8, cut the rest of the read)
    appended += " MINLEN:200" // the resulting reads must have a length of at least 36

    // setup the memory limits, high for trimmomatic (java + SGE = weird memory issues...)
    this.memoryLimit = 4
    this.residentRequest = 4
    this.residentLimit = 4

    // CMD command issued, and the hidden queue output file names
    var cmd = "java -Xmx3g -jar " + trimmomaticPath + (if (fqs.length == 2) " PE" else " SE") + " -phred33 -threads 1 " + fqs.mkString(" ")
    if (fqs.length == 2) cmd += " " + fqOuts(0) + " " + fqOutsUnpaired(0) + " " + fqOuts(1) + " " + fqOutsUnpaired(1) + " " + appended
    else cmd += " " + fqOuts(0) + " " + " " + appended

    def commandLine = cmd
    this.isIntermediate = false
    this.analysisName = queueLogDir + fqs(0) + ".trimmomatic"
    this.jobName = queueLogDir + fqs(0) + ".trimmomatic"

  }

  /**
   * seqprep -- trimmomatic's main competition, does read merging as well
   *
   * trims off adapter sequence and low quality regions in paired end sequencing data, merges on reqest
   */
  // ********************************************************************************************************
  case class SeqPrep(inFastqs: List[File], outs: List[File], outputMerged: File ) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input FASTAs") var fqs = inFastqs
    @Output(doc = "output fastas (corrected)") var fqOuts = outs
    @Output(doc = "output merged reads") var merged = outputMerged
    
    if (inFastqs.length != 2 || outs.length != 2)
      throw new IllegalArgumentException("Seqprep can only be run on paired end sequencing! for input files " + inFastqs.mkString(", "))

    this.memoryLimit = 4
    this.residentRequest = 4
    this.residentLimit = 4

    // o -> minimum overlap, n -> faction of bases that must match in overlap
    var cmd = seqPrepPath + " -o 25 -n 0.95 -f " + fqs(0) + " -r " + fqs(1) + " -1 " + fqOuts(0) + " -2 " + fqOuts(1) + " -s " + merged
    
    def commandLine = cmd
    this.isIntermediate = false
      this.analysisName = queueLogDir + fqs(0) + ".seqprep"
    this.jobName = queueLogDir + fqs(0) + ".seqprep"
  }

}


