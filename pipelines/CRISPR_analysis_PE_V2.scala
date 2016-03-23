/** *********************************************
  * Requires Java 1.7. On the UW machines this means
  * running at least: "module load java/7u17", or the appropriate
  * module if there's a more recent version, to load java 1.7 into
  * your environment.  The best idea is to place this into your
  * bash/shell profile
  *
  *
  * Copyright (c) 2014, 2015, 2016, aaronmck
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
  * @date June, 2015
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

import collection.JavaConversions._
import htsjdk.samtools.SAMFileReader
import htsjdk.samtools.SAMFileHeader.SortOrder
import scala.io.Source
import scala.collection.immutable._
import java.io.File

/**
  * Given amplicon reads from a CRISPR experiment, align and call events over those cut sites.
  *
  * PLEASE READ:
  * 
  * - This pipeline assumes your reference file, say myref.fa, has the following additional files in the same dir:
  * - myref.fa.cutSites <-- contains the CRISPR target seq, position start, the cutsite position
  * - myref.fa.primers <-- two lines, containing the positive strand primers on both ends of your amplicon. 
  * 
 **/
class DNAQC extends QScript {
  qscript =>

  /** **************************************************************************
    * Data Parameters
    * ************************************************************************** */
  @Input(doc = "Tab-delimited tear sheet containing sample information", fullName = "input", shortName = "i", required = true)
  var input: File = _

  @Input(doc = "where to put aggregated stats", fullName = "aggLocation", shortName = "agg", required = true)
  var aggregateLocation: File = _

  @Argument(doc = "the experiment name, used to generate the base output folder on the web server", fullName = "expName", shortName = "expName", required = true)
  var experimentalName: String = ""

  @Argument(doc = "do we NOT want to use Trimmomatic to clean the reads of bad sequence?", fullName = "dontTrim", shortName = "dontTrim", required = false)
  var dontTrim: Boolean = false

  /** **************************************************************************
    * Optional Parameters -- control parameters for the script (alignment, etc)
    * ************************************************************************** */

  @Argument(doc = "read length (how long is each read in a pair (300 cycle kit -> 150bp reads = 150)", fullName = "readLength", shortName = "readLength", required = false)
  var readLength: Int = 250

  @Argument(doc = "are we runing with UMI reads? If so set this value to > 0, representing the UMI length", fullName = "umi", shortName = "umi", required = false)
  var umiData: Int = 0

  @Argument(doc = "the number of UMIs required to call a successful UMI capture event, if you're using UMIs", fullName = "umiCount", shortName = "umiCount", required = false)
  var minUMIs = 10

  @Input(doc = "where to put the web files", fullName = "web", shortName = "www", required = false)
  var webSite: File = "/net/shendure/vol2/www/content/members/aaron/crispr/"

  /** **************************************************************************
    * Path parameters -- where to find tools
    * ************************************************************************** */

  @Input(doc = "The path to the Trimmomatic adapters file", fullName = "trimadapters", shortName = "tad", required = false)
  var adaptersFile: File = new File("/net/gs/vol1/home/aaronmck/tools/trimming/Trimmomatic-0.32/adapters/TruSeq3-PE.fa")

  @Input(doc = "The path to the binary of fastqc", fullName = "fastqc", shortName = "fqc", required = false)
  var fastqcPath: File = "/net/gs/vol1/home/aaronmck/tools/fastqc/FastQC/fastqc"

  @Input(doc = "The path to the jar file for trimmomatic", fullName = "trim", shortName = "tm", required = false)
  var trimmomaticPath: File = "/net/gs/vol1/home/aaronmck/tools/trimming/Trimmomatic-0.32/trimmomatic-0.32.jar"

  @Input(doc = "The path to the barcode splitter", fullName = "maulpath", shortName = "mlp", required = false)
  var maulPath: File = "/net/gs/vol1/home/aaronmck/source/Maul/target/scala-2.10/Maul-assembly-1.0.jar"

  @Input(doc = "The path to the seqprep tool", fullName = "seqprep", shortName = "sprep", required = false)
  var seqPrepPath: File = "/net/gs/vol1/home/aaronmck/tools/bin/SeqPrep"

  @Input(doc = "The path to the flash tool", fullName = "flash", shortName = "flash", required = false)
  var flashPath: File = "/net/gs/vol1/home/aaronmck/tools/bin/flash"

  @Input(doc = "the script to analyze the edits made to a crispr target", fullName = "edits", shortName = "edits", required = false)
  var crisprPath: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/analyze_diversity_of_edits.scala"

  @Input(doc = "the script to analyze the per site edits", fullName = "writeweb", shortName = "wwwrite", required = false)
  var writeWebFiles: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/scripts/write_web_files.scala"

  @Input(doc = "the location of the needleall program", fullName = "needle", shortName = "needle", required = false)
  var needlePath: File = "/net/gs/vol1/home/aaronmck/tools/bin/needleall"

  @Input(doc = "the path the javascript conversion script", fullName = "JSTable", shortName = "JSTable", required = false)
  var toJSTableScript: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/stats_to_javascript_tables2.scala"

  @Argument(doc = "zip two reads together", fullName = "ZipReads", shortName = "ZipReads", required = false)
  var zipReadsPath = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/zip_two_read_files.scala"

  @Argument(doc = "move the data over to the web location", fullName = "webpub", shortName = "webpub", required = false)
  var toWebPublishScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/push_to_web_location.scala"

  @Argument(doc = "the web location", fullName = "webpubLoc", shortName = "webpubLoc", required = false)
  var webLocation = "/net/shendure/vol2/www/content/members/aaron/fate_map/"

  @Argument(doc = "aggregate stats files (UMIed reads) together", fullName = "stats", shortName = "stats", required = false)
  var aggregateScripts = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/aggregate_stats.scala"

  @Argument(doc = "analyize stats file and generate a bunch of plots", fullName = "plotsScript", shortName = "plotsScript", required = false)
  var statsFileAnalysis = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/stats_file_analysis.py"

  @Argument(doc = "the first adapter sequence", fullName = "adaptOne", shortName = "adaptOne", required = false)
  var adapterOne = "GACCTCGAGACAAATGGCAG"

  @Argument(doc = "the second adapter sequence", fullName = "adaptTwo", shortName = "adaptTwo", required = false)
  var adapterTwo = "CGAAGCTTGAGCTCGAGATCTG"

  /** **************************************************************************
    * Global Variables
    * ************************************************************************** */

  // Gracefully hide Queue's output -- this is very important otherwise there's Queue output everywhere and it's a mess
  val queueLogDir: String = ".qlog/"

  // ** **************************************************************************
  // * Main script entry point
  // * **************************************************************************
  def script() {
    var umiStatsFiles = List[File]()
    var nonUmiStatsFiles = List[File]()

    val sampleWebBase =      dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName), "our output web publishing directory")
    val aggReportDir =       dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + "aggregateHMIDReport"), "our report output directory")
    val aggWebTreeLocation = dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + "tree"), "our output tree directory")

    // read in the tear sheet and process each sample
    parseTearSheet(input).foreach(sampleObj => {

      // the sample name we'll use for generating files
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
      

      // create the per-sample output directories, verifying that they either exist already or can be made, also make the base web dir, it's ok that this is
      // duplicated, as don't try to make it if it exists
      val webLoc =             webSite + File.separator + sampleObj.sample + File.separator
      val sampleOutput =       dirOrCreateOrFail(new File(sampleObj.outputDir + File.separator + sampleObj.sample), "sample output directory")
      val initialFastQCDir =   dirOrCreateOrFail(new File(sampleOutput + File.separator + "initial_fastQC"), "initial fastqc directory")
      val initialFastQMerged = dirOrCreateOrFail(new File(sampleOutput + File.separator + "initial_fastQC_Merged"), "initial fastqc directory")
      val postCleanFastQCDir = dirOrCreateOrFail(new File(sampleOutput + File.separator + "post_clean_fastQC"), "final fastqc directory")
      val sampleWebLocation =  dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + sampleTag), "our output web publishing directory")

      // our barcode split files
      var barcodeSplit1 = new File(sampleOutput + File.separator + sampleTag + ".barcodeSplit.fastq1.fq.gz")
      var barcodeSplit2 = new File(sampleOutput + File.separator + sampleTag + ".barcodeSplit.fastq2.fq.gz")

      val samUnmergedFasta = new File(sampleOutput + File.separator + sampleTag + ".UM.fasta")
      val samMergedFasta = new File(sampleOutput + File.separator + sampleTag + ".M.fasta")

      var mergedReads = new File(sampleOutput + File.separator + "out.extendedFrags.fastq")
      var mergedReadUnzipped = new File(sampleOutput + File.separator + sampleTag + ".merged.fq")
      var unmergedUnzipped = new File(sampleOutput + File.separator + sampleTag + ".unmerged.fq")

      // ************************************** split the input files by barcodes, either one or two **************************************
      var barcodeFiles = if (dualBarcode) List[File](sampleObj.fastqBarcode1,sampleObj.fastqBarcode2) else List[File](sampleObj.fastqBarcode1)
      var inputFiles = if (pairedEnd) List[File](sampleObj.fastq1,sampleObj.fastq2) else List[File](sampleObj.fastq1)
      var processedFastqs = if (pairedEnd) List[File](barcodeSplit1,barcodeSplit2) else List[File](barcodeSplit1)
      val barcodeConfusion = new File(sampleOutput + File.separator + sampleTag + ".barcodeConfusion")
      val barcodeStats = new File(sampleOutput + File.separator + sampleTag + ".barcodeStats")
      val overlapFile = new File(sampleOutput + File.separator + sampleTag + ".readOverlap")

      if (sampleObj.barcode1.toLowerCase != "all" && sampleObj.barcode2.toLowerCase != "all") {
        add(Maul(inputFiles, barcodeFiles, List[String](sampleObj.barcode1,sampleObj.barcode2), processedFastqs, barcodeStats, barcodeConfusion,overlapFile))
        inputFiles = List[File](barcodeSplit1,barcodeSplit2)
      }

      val cleanedFastqs = List[File](
        new File(sampleOutput + File.separator + "out.notCombined_1.fastq"),
        new File(sampleOutput + File.separator + "out.notCombined_2.fastq"))

      // run an Fastqc so we know what we started with
      // ************************************************************************************
      add(Fastqc(inputFiles, initialFastQCDir))

      val toAlignFastq1 = new File(sampleOutput + File.separator + sampleTag + ".fwd.fastq")
      val toAlignFastq2 = new File(sampleOutput + File.separator + sampleTag + ".rev.fastq")
      val toAlignStats = new File(sampleOutput + File.separator + sampleTag + ".stats")
      
      val perBaseEventFile = new File(sampleOutput + File.separator + sampleTag + ".perBase")
      val topReadFile = new File(sampleOutput + File.separator + sampleTag + ".topReadEvents")
      val topReadFileNew = new File(sampleOutput + File.separator + sampleTag + ".topReadEventsNew")
      val topReadCount = new File(sampleOutput + File.separator + sampleTag + ".topReadCounts")
      val allReadCount = new File(sampleOutput + File.separator + sampleTag + ".allReadCounts")
      val cutSites = new File(sampleObj.reference + ".cutSites")
      val unpairedReads = List[File](new File(sampleOutput + File.separator + sampleTag + ".unpaired1.fq.gz"),new File(sampleOutput + File.separator + sampleTag + ".unpaired2.fq.gz"))

      // *******************************************************************************************
      // from here on we propagate the cleanedTmp as the input to the next step
      // *******************************************************************************************
      var cleanedTmp = List[File](new File(sampleOutput + File.separator + sampleTag + ".cleanedInter1.fq.gz"),new File(sampleOutput + File.separator + sampleTag + ".cleanedInter2.fq.gz"))


      if (!dontTrim) { // yes a double negitive...sorry: if we want to trim, do this
        add(Trimmomatic(inputFiles,adaptersFile,cleanedTmp,unpairedReads))
      } else {
        // just pass through the input files as 'cleaned'
        cleanedTmp = inputFiles
      }

      // if we're a UMI run, we divert here to merge the reads by UMIs, and sub back in the merged
      // high-quality reads into the pipeline as if there hadn't been UMIs
      if (sampleObj.UMIed) {
        val toAlignFastq1 = new File(sampleOutput + File.separator + sampleTag + ".umi.fwd.fastq")
        val toAlignFastq2 = new File(sampleOutput + File.separator + sampleTag + ".umi.rev.fastq")

        // collapse the reads by their UMI and output FASTA files with the remaining quality reads
        add(UMIProcessingPaired(
          cleanedTmp(0),
          cleanedTmp(1),
          toAlignFastq1,
          toAlignFastq2,
          10,
          new File(sampleObj.reference + ".primers"),
          sampleObj.sample))

        cleanedTmp = List[File](toAlignFastq1,toAlignFastq2)
      }

      //add(SeqPrep(cleanedTmp, cleanedFastqs, mergedReads, adapterOne, adapterTwo))out.extendedFrags.
      add(Flash(cleanedTmp, cleanedFastqs, mergedReads, sampleOutput + File.separator))

      add(ZipReadFiles(cleanedFastqs(0), cleanedFastqs(1), unmergedUnzipped))
      //add(Gunzip(mergedReads, mergedReadUnzipped))

      add(NeedlemanAllFasta(sampleObj.reference, unmergedUnzipped, samUnmergedFasta, false))
      add(NeedlemanAllFasta(sampleObj.reference, mergedReads, samMergedFasta, false)) // mergedReadUnzipped

      add(ReadsToStats(samUnmergedFasta,
        samMergedFasta,
        toAlignStats,
        cutSites,
        new File(sampleObj.reference + ".primers"),
        sampleObj.sample))

      nonUmiStatsFiles :+= toAlignStats

      add(ToJavascriptTables(toAlignStats, cutSites, perBaseEventFile, topReadFile, topReadCount,allReadCount, topReadFileNew))
      add(ToWebPublish(sampleWebLocation, perBaseEventFile, topReadFileNew, topReadCount, cutSites, allReadCount))
    })

    // agg. all of the stats together into a single file
    val mergedStatsFile = new File(aggregateLocation + File.separator + "merged.stats")
    val mergedUMIInfo = new File(aggregateLocation + File.separator + "merged.info")
    add(AggregateStatsFiles(umiStatsFiles, mergedStatsFile, mergedUMIInfo))
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

  // The storage container for the data we've read from the input tear sheet
  case class SourceEntry(
    val sample: String,
    val UMIed: Boolean,
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
    if (inputLines.next().stripLineEnd != "sample\tumi\treference\toutput.dir\tfastq1\tfastq2\tbarcode.fastq1\tbarcode.fastq2\tbarcode1\tbarcode2")
      throw new IllegalArgumentException("Your header doesn't seem like a correctly formatted  tear sheet!")

    return inputLines.map(line => {
      val tokens = line.split( """\s+""") // Use a regex to avoid empty tokens
      try {
        (new SourceEntry(
          tokens(0), // sample
          tokens(1).toBoolean, // case control status: true if control, else a name of a sample for the control
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

  /** **************************************************************************
    * Classes (non-GATK programs)
    *
    * a note here: normally we could use the built in Picard (and GATK) functions,
    * but for some reason our SGE interacts very poorly with Java, and we need to
    * build in an even higher 'buffer' of memory (the difference between the
    * requested memory to SGE and Java's Xmx parameter)
    * ************************************************************************** */

  // ********************************************************************************************************
  case class ToJavascriptTables(statsFile: File, cutSites: File, perBaseEvents: File, topReads: File, topReadCounts: File, allEventCounts: File, perBaseEventsNew: File) extends ExternalCommonArgs {
    @Input(doc = "the input stats file") var stats = statsFile
    @Input(doc = "the  cut sites information file") var cuts = cutSites
    @Output(doc = "per base event counts") var perBase = perBaseEvents
    @Output(doc = "the top reads and their base by base events") var topR = topReads
    @Output(doc = "the the counts for all the top reads") var topReadC = topReadCounts
    @Output(doc = "the the counts for all the reads") var allReadC = allEventCounts
    @Output(doc = "new per base event style") var perBaseES = perBaseEventsNew

    def commandLine = "scala -J-Xmx6g " + toJSTableScript + " " + stats + " " + perBase + " " + topR + " " + topReadC + " " + allReadC + " " + cuts + " " + perBaseES

    this.analysisName = queueLogDir + stats + ".toJS"
    this.jobName = queueLogDir + stats + ".toJS"
    this.memoryLimit = 8
    this.residentRequest = 8
    this.residentLimit = 8
  }

  // ********************************************************************************************************
  case class ToWebPublish(webLocation: File, perBaseEvents: File, topReads: File, topReadCounts: File, cutSites: File, allReadCount: File) extends ExternalCommonArgs {
    @Input(doc = "per base event counts") var webL = webLocation
    @Input(doc = "per base event counts") var perBase = perBaseEvents
    @Input(doc = "the top reads and their base by base events") var topR = topReads
    @Input(doc = "the counts for all the top reads") var topReadC = topReadCounts
    @Input(doc = "the cutsites") var cuts = cutSites
    @Input(doc = "all HMIDs") var allReads = allReadCount

    def commandLine = "scala -J-Xmx1g " + toWebPublishScript + " " + webL + " " + perBase + " " + topR + " " + topReadC + " " + cuts + " " + allReads

    this.analysisName = queueLogDir + perBase + ".web"
    this.jobName = queueLogDir + perBase + ".web"
    this.memoryLimit = 2
    this.isIntermediate = false
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

  // Needleman-Wunsch aligner for single ended (or merged) reads
  // ********************************************************************************************************
  case class NeedlemanAll(reference: File, inFastq: File, outputSam: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the merged read fastq") var fq = inFastq
    @Argument(doc = "the reference fasta/fa") var ref = reference
    @Output(doc = "the output sam file") var outSam = outputSam

    def commandLine = needlePath + " -endopen 10.0 -endextend 0.01 -aformat3 sam -gapextend 0.1 -gapopen 10.0 -asequence " + ref + " -bsequence " + fq + " -outfile " + outSam

    this.analysisName = queueLogDir + outSam + ".needle"
    this.jobName = queueLogDir + outSam + ".needle"
  }

  // Needleman-Wunsch aligner for single ended (or merged) reads
  // ********************************************************************************************************
  case class NeedlemanAllFasta(reference: File, inFastq: File, outputFasta: File, reverse: Boolean) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the merged read fastq") var fq = inFastq
    @Argument(doc = "the reference fasta/fa") var ref = reference
    @Output(doc = "the output fasta file") var outFasta = outputFasta

    //-sreverse2
    def commandLine = needlePath + " -datafile /net/shendure/vol10/projects/CRISPR.lineage/nobackup/reference_data/EDNAFULL -snucleotide1 -snucleotide2 " + (if (reverse) "-sreverse2 " else " ") + "-aformat3 fasta -gapextend 0.5 -gapopen 10.0 -asequence " + ref + " -bsequence " + fq + " -outfile " + outFasta

    this.analysisName = queueLogDir + outFasta + ".needle"
    this.jobName = queueLogDir + outFasta + ".needle"
  }

  // Needleman-Wunsch aligner needs uncompressed input -- stupid aligner
  // ********************************************************************************************************
  case class Gunzip(inFile: File, outFile: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the merged read fastq") var inF = inFile
    @Output(doc = "the output sam file") var outF = outFile

    def commandLine = "gunzip -c " + inF + " > " + outFile

    this.analysisName = queueLogDir + outF + ".gunzip"
    this.jobName = queueLogDir + outF + ".gunzip"
  }

  // interleave the two read files to better keep track of the reads post alignment
  // ********************************************************************************************************
  case class ZipReadFiles(inRead1File: File, inRead2File: File, outFile: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the read 1") var inF1 = inRead1File
    @Input(doc = "the read 2") var inF2 = inRead2File
    @Output(doc = "the output sam file") var outF = outFile

    def commandLine = "scala " + zipReadsPath + " " + inF1 + " " + inF2 + " " + outF

    this.analysisName = queueLogDir + outF + ".gunzip"
    this.jobName = queueLogDir + outF + ".gunzip"
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

    def commandLine = "scala -J-Xmx4g -cp /net/gs/vol1/home/aaronmck/source/sandbox/aaron/lib/target/scala-2.11/sequencing-lib-assembly-1.0.jar " + crisprPath + " " + outa + " " + outcase + " " + outcont + " " + ref + ".primers " + ref + ".cutSites " + bCut + " " + bControl + " " + matchedReads + " " + cutData

    this.analysisName = queueLogDir + outputAllele + ".crisprEdits"
    this.jobName = queueLogDir + outputAllele + ".crisprEdits"
    this.isIntermediate = false
  }

  // aggregate stats flies down to a single file
  // ********************************************************************************************************
  case class AggregateStatsFiles(inputFiles: List[File], outputStatsFile: File, outputUmiStatsFile: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the umi input file") var inputFls = inputFiles
    @Output(doc = "the output merged stats file") var outStats = outputStatsFile
    @Output(doc = "the output umi information") var outUMIStats = outputUmiStatsFile

    def commandLine = "scala -J-Xmx8g " + aggregateScripts + " " + inputFls.mkString(",") + " " + outStats + " " + outUMIStats)

    this.analysisName = queueLogDir + outStats + ".outStats"
    this.jobName = queueLogDir + outStats + ".outStats"
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


  // 
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


    def commandLine = "scala -J-Xmx4g -cp /net/gs/vol1/home/aaronmck/source/sandbox/aaron/lib/target/scala-2.11/sequencing-lib-assembly-1.0.jar " + writeWebFiles + " " + ref + " " + refcut + " " + casebam + " " + controlbam + " " + " " + caseFile + " " + controlFile + " " + ref + ".primers " + outDir + " " + htmlFile + " " + figDir + " " + webConfig + " " + diag

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
    appended += " LEADING:10" // trim off leading bases with a quality less than 10
    appended += " TRAILING:10" // trim off trailing bases with a quality less than 10
    appended += " SLIDINGWINDOW:5:15" // remove bases with a quality less than X in a sliding window of size 
    appended += " MINLEN:75" // the resulting reads must have a length of at least X

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
   * Process the UMIs, merging and aligning reads down to a single, high-quality concensus per UMI
   */
  // ********************************************************************************************************
  case class UMIProcessingPaired(inMergedReads1: File, inMergedReads2: File, outputFASTA1: File, outputFASTA2: File,
    umiCutOff: Int, primersFile: File, sampleName: String) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input reads (fwd)") var inReads1 = inMergedReads1
    @Input(doc = "input reads (rev)") var inReads2 = inMergedReads2
    @Output(doc = "output fasta for further alignment (fwd)") var outFASTA1 = outputFASTA1
    @Output(doc = "output fasta for further alignment (rev)") var outFASTA2 = outputFASTA2
    @Argument(doc = "how many UMIs do we need to initial have to consider merging them") var umiCut = umiCutOff
    @Argument(doc = "the primers file; one line per primer that we expect to have on each end of the resulting merged read") var primers = primersFile
    @Argument(doc = "the sample name") var sample = sampleName

    var cmdString = "scala -J-Xmx13g /net/shendure/vol10/projects/CRISPR.lineage/nobackup/bin/UMIMerge.jar "
    cmdString += " --inputFileReads1 " + inReads1 + " --inputFileReads2 " + inReads2 + " --outputFastq1 " + outFASTA1 + " --outputFastq2 " + outFASTA2 
    cmdString += " --umiLength " + umiCut + " --primersEachEnd " + primers + " --samplename " + sample + " --umiStart 0 --minimumUMIReads " + minUMIs

    var cmd = cmdString

    this.memoryLimit = 16
    this.residentRequest = 16
    this.residentLimit = 16

    def commandLine = cmd
    this.isIntermediate = false
    this.analysisName = queueLogDir + inReads1 + ".umis"
    this.jobName = queueLogDir + inReads1 + ".umis"
  }

  //--inputUnmerged --inputMerged --outputStats --cutSites --primersEachEnd --sample test
  case class ReadsToStats(inputUnmerged: File,
    inputMerged: File,
    outputStats: File,
    cutSites: File,
    primersEachEnd: File,
    sampleName: String) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input reads (fwd)") var unmerged = inputUnmerged
    @Input(doc = "input reads (rev)") var merged = inputMerged
    @Input(doc = "the cutsite locations") var cutSiteFile = cutSites
    @Output(doc = "output statistics file for containing information about the UMI merging process") var outStat = outputStats
    @Argument(doc = "the primers file; one line per primer that we expect to have on each end of the resulting merged read") var primers = primersEachEnd
    @Argument(doc = "the sample name") var sample = sampleName

    var cmdString = "java -jar -Xmx2g /net/shendure/vol10/projects/CRISPR.lineage/nobackup/bin/ReadToStats.jar "
    cmdString += " --inputUnmerged " + inputUnmerged + " --inputMerged " + inputMerged + " --cutSites "
    cmdString += cutSiteFile + " --outputStats "
    cmdString += outStat + " --primersEachEnd " + primers + " --sample "
    cmdString += sample

    var cmd = cmdString

    this.memoryLimit = 6
    this.residentRequest = 6
    this.residentLimit = 6

    def commandLine = cmd
    this.isIntermediate = false
    this.analysisName = queueLogDir + outStat + ".calls"
    this.jobName = queueLogDir + outStat + ".calls"
  }

  /**
   * seqprep -- trimmomatic's main competition, does read merging as well
   *
   * trims off adapter sequence and low quality regions in paired end sequencing data, merges on reqest
   */
  // ********************************************************************************************************
  case class SeqPrep(inFastqs: List[File], outs: List[File], outputMerged: File, adapterOne: String, adapterTwo: String) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input FASTAs") var fqs = inFastqs
    @Output(doc = "output fastas (corrected)") var fqOuts = outs
    @Output(doc = "output merged reads") var merged = outputMerged

    if (inFastqs.length != 2 || outs.length != 2)
      throw new IllegalArgumentException("Seqprep can only be run on paired end sequencing! for input files " + inFastqs.mkString(", "))

    this.memoryLimit = 4
    this.residentRequest = 4
    this.residentLimit = 4

    // o -> minimum overlap, n -> faction of bases that must match in overlap
    var cmd = seqPrepPath + " -f " + fqs(0) + " -r " + fqs(1) + " -1 " + fqOuts(0) + " -2 " + fqOuts(1) + " -s " + merged + " -A " + adapterOne + " -B " + adapterTwo

    def commandLine = cmd
    this.isIntermediate = false
      this.analysisName = queueLogDir + fqs(0) + ".seqprep"
    this.jobName = queueLogDir + fqs(0) + ".seqprep"
  }

  /**
   * flash -- seems to do a much better job with read overlap assembly than SeqPrep or Trimmomatic
   * this is one of those tools where you can't specify the output file only the directory, so our 
   * output file names have to be correctly formatted for what flash will output
   */
  // ********************************************************************************************************
  case class Flash(inFastqs: List[File], outs: List[File], outputMerged: File, outputDir: File) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input FASTAs") var fqs = inFastqs
    @Output(doc = "output unmerged files") var fqOuts = outs
    @Output(doc = "output merged reads") var merged = outputMerged
    @Argument(doc = "the output directory") var outputDr = outputDir

    if (inFastqs.length != 2 || outs.length != 2)
      throw new IllegalArgumentException("Seqprep can only be run on paired end sequencing! for input files " + inFastqs.mkString(", "))

    this.memoryLimit = 4
    this.residentRequest = 4
    this.residentLimit = 4

    var cmd = flashPath + " --min-overlap 30 --max-mismatch-density 0.02 --output-directory=" + outputDr + " " + fqs(0) + " " + fqs(1) 

    def commandLine = cmd
    this.isIntermediate = false
      this.analysisName = queueLogDir + fqs(0) + ".flash"
    this.jobName = queueLogDir + fqs(0) + ".flash"
  }
}
