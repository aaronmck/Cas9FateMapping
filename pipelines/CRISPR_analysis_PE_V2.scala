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
 * Quality control sequencing data and optionally align the data to the genome
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

  @Argument(doc = "the mapping from sample name to biologically meaningful clade", fullName = "cladeAssignments", shortName = "cladeAssignments", required = false)
  var cladeAssignments: File = ""


  /** **************************************************************************
    * Optional Parameters -- control parameters for the script (alignment, etc)
    * ************************************************************************** */

  @Argument(doc = "read length (how long is each read in a pair (300 cycle kit -> 150bp reads = 150)", fullName = "readLength", shortName = "readLength", required = false)
  var readLength: Int = 250

  @Argument(doc = "are we runing with UMI reads? If so set this value to > 0, representing the UMI length", fullName = "umi", shortName = "umi", required = false)
  var umiData: Int = 0

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
  var bwaPath: File = "/net/gs/vol3/software/modules-sw/bwa/0.7.3/Linux/RHEL6/x86_64/bin/bwa"

  @Input(doc = "The path to the barcode splitter", fullName = "maulpath", shortName = "mlp", required = false)
  var maulPath: File = "/net/gs/vol1/home/aaronmck/source/Maul/target/scala-2.10/Maul-assembly-1.0.jar"

  @Input(doc = "The path to the seqprep tool", fullName = "seqprep", shortName = "sprep", required = false)
  var seqPrepPath: File = "/net/gs/vol1/home/aaronmck/tools/bin/SeqPrep"

  @Input(doc = "the script to analyze the edits made to a crispr target", fullName = "edits", shortName = "edits", required = false)
  var crisprPath: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/analyze_diversity_of_edits.scala"

  @Input(doc = "the script to analyze the per site edits", fullName = "posScript", shortName = "posScript", required = false)
  var positionAnalysisPath: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/scripts/sam_to_positions.scala"

  @Input(doc = "the script to analyze the per site edits", fullName = "writeweb", shortName = "wwwrite", required = false)
  var writeWebFiles: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/scripts/write_web_files.scala"

  @Input(doc = "the location of the smith waterman program", fullName = "swp", shortName = "swp", required = false)
  var smithWaterman: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/scripts/optimize_alignment.scala"

  @Input(doc = "the location of the script to strip UMIs", fullName = "sumi", shortName = "sumi", required = false)
  var stripUMI: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/scripts/barcode_analysis.scala"

  @Input(doc = "the location of the needleall program", fullName = "needle", shortName = "needle", required = false)
  var needlePath: File = "/net/gs/vol1/home/aaronmck/tools/bin/needleall"

  @Input(doc = "the cut site analysis script", fullName = "csa", shortName = "csa", required = false)
  var cutSiteScript: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPR/scripts/cut_site_analysis.Rmd"

  @Input(doc = "find the reference from the first reads", fullName = "cns", shortName = "cns", required = false)
  var concReadScript: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/conc_read.scala"

  @Input(doc = "filter out the reads with inconsistant cigar strings", fullName = "ft", shortName = "ft", required = false)
  var filterSamScript: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/filter_sam.scala"

  @Input(doc = "where to put the web files", fullName = "web", shortName = "www", required = false)
  var webSite: File = "/net/shendure/vol2/www/content/members/aaron/crispr/"

  @Input(doc = "where to find the diagram files", fullName = "diag", shortName = "diag", required = false)
  var diagramLoc: File = "/net/gs/vol1/home/aaronmck/source/sandbox/aaron/projects/CRISPRScripts/"

  @Input(doc = "where we can find the UMI counting script", fullName = "umiScript", shortName = "umiScript", required = false)
  var umiScript: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/umi_count.scala"

  @Input(doc = "where we can find the paired-end UMI counting script", fullName = "umiScriptPE", shortName = "umiScriptPE", required = false)
  var umiScriptPE: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/umi_count_PE.scala"

  @Input(doc = "the path to PandaSeq", fullName = "PandaSeq", shortName = "PandaSeq", required = false)
  var pandaseqPath: File = "/net/gs/vol3/software/modules-sw/PANDAseq/2.8.1/Linux/RHEL6/x86_64/bin/pandaseq"

  @Input(doc = "the path the javascript conversion script", fullName = "JSTable", shortName = "JSTable", required = false)
  var toJSTableScript: File = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/stats_to_javascript_tables.scala"

  @Argument(doc = "zip two reads together", fullName = "ZipReads", shortName = "ZipReads", required = false)
  var zipReadsPath = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/zip_two_read_files.scala"

  @Argument(doc = "move the data over to the web location", fullName = "webpub", shortName = "webpub", required = false)
  var toWebPublishScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/push_to_web_location.scala"

  @Argument(doc = "the web location", fullName = "webpubLoc", shortName = "webpubLoc", required = false)
  var webLocation = "/net/shendure/vol2/www/content/members/aaron/fate_map/"

  @Argument(doc = "perform NJ and distance matrix calculation with this tool", fullName = "NJDist", shortName = "NJDist", required = false)
  var distanceMatrix = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/ConstrainedNeighborJoin/target/scala-2.11/ConstrainedNJ-assembly-1.0.jar"

  @Argument(doc = "perform NJ and plotting in R with this script", fullName = "NJScript", shortName = "NJScript", required = false)
  var njScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/nj_from_distance_matrix.R"

  @Argument(doc = "aggregate stats files (UMIed reads) together", fullName = "stats", shortName = "stats", required = false)
  var aggregateScripts = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/aggregate_stats.scala"

  @Argument(doc = "analyize stats file and generate a bunch of plots", fullName = "plotsScript", shortName = "plotsScript", required = false)
  var statsFileAnalysis = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/stats_file_analysis.py"

  @Argument(doc = "Perform basic statistics and output plots about the reads / callsets with this script", fullName = "postProcess", shortName = "postProcess", required = false)
  var postProcess = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/post_process_alignments.R"

  @Argument(doc = "The script the generates a report at the end of the run", fullName = "reportScript", shortName = "reportScript", required = false)
  var aggReportScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/aggregate_report.Rmd"

  @Argument(doc = "the script that takes a tree file and ", fullName = "treeCopy", shortName = "treeCopy", required = false)
  var treeCopyScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/prepare_tree.scala"

  @Argument(doc = "The script that turns a distance matrix into a tree", fullName = "createTree", shortName = "cTree", required = false)
  var createTreeScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/create_tree.R"

  @Argument(doc = "The four gamete test", fullName = "fourGamete", shortName = "fourGamete", required = false)
  var fourGameteScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/four_gamete_test.scala"

  @Argument(doc = "create a file for use in PHYLIP's mix", fullName = "toMix", shortName = "toMix", required = false)
  var toMixScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/stats_to_mix_input.scala"

  /** **************************************************************************
    * Global Variables
    * ************************************************************************** */

  // Gracefully hide Queue's output -- this is very important otherwise there's Queue output everywhere and it's a mess
  val queueLogDir: String = ".qlog/"

  // ** **************************************************************************
  // * Main script
  // * **************************************************************************
  def script() {
    var umiStatsFiles = List[File]()
    var nonUmiStatsFiles = List[File]()

    val sampleWebBase =      dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName), "our output web publishing directory")
    val aggReportDir =       dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + "aggregateHMIDReport"), "our report output directory")
    val aggWebTreeLocation = dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + "tree"), "our output tree directory")

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
      

      // create the per-sample output directories, verifying that they either exist already or can be made, also make the base web dir, it's ok that this is
      // duplicated, as don't try to make it if it exists
      val webLoc =             webSite + File.separator + sampleObj.sample + File.separator
      val sampleOutput =       dirOrCreateOrFail(new File(sampleObj.outputDir + File.separator + sampleObj.sample), "sample output directory")
      val initialFastQCDir =   dirOrCreateOrFail(new File(sampleOutput + File.separator + "initial_fastQC"), "initial fastqc directory")
      val initialFastQMerged = dirOrCreateOrFail(new File(sampleOutput + File.separator + "initial_fastQC_Merged"), "initial fastqc directory")
      val postCleanFastQCDir = dirOrCreateOrFail(new File(sampleOutput + File.separator + "post_clean_fastQC"), "final fastqc directory")
      val webTreeLocation =    dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + sampleTag + File.separator + "tree"), "our output tree directory")
      val sampleWebLocation =  dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + sampleTag), "our output web publishing directory")
      val plotsDir =           dirOrCreateOrFail(new File(webLocation +  File.separator + experimentalName + File.separator + sampleTag + File.separator + "umiPlots"), "our output plot directory")

      // our barcode split files
      var barcodeSplit1 = new File(sampleOutput + File.separator + sampleTag + ".barcodeSplit.fastq1.fq.gz")
      var barcodeSplit2 = new File(sampleOutput + File.separator + sampleTag + ".barcodeSplit.fastq2.fq.gz")

      val samUnmergedFasta = new File(sampleOutput + File.separator + sampleTag + ".UM.fasta")
      val samMergedFasta = new File(sampleOutput + File.separator + sampleTag + ".M.fasta")

      var mergedReads = new File(sampleOutput + File.separator + sampleTag + ".merged.fq.gz")
      var mergedReadUnzipped = new File(sampleOutput + File.separator + sampleTag + ".merged.fq")
      var unmergedUnzipped = new File(sampleOutput + File.separator + sampleTag + ".unmerged.fq")

      val mergedUMIStripped = new File(sampleOutput + File.separator + sampleTag + ".merged.stripedUMI.fq.gz")
      val mergedUMIStrippedStats = new File(webLoc + "/UMI_stats.txt")

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

      val postProcessedFastq1 = new File(sampleOutput + File.separator + sampleTag + ".cleaned.trim.fastq1.fq.gz")
      val postProcessedFastq2 = new File(sampleOutput + File.separator + sampleTag + ".cleaned.trim.fastq2.fq.gz")
      val cleanedFastqs = List[File](postProcessedFastq1,postProcessedFastq2)

      // run an Fastqc so we know what we started with
      // ************************************************************************************
      add(Fastqc(inputFiles, initialFastQCDir))

      val toAlignFastq1 = new File(sampleOutput + File.separator + sampleTag + ".fwd.fastq")
      val toAlignFastq2 = new File(sampleOutput + File.separator + sampleTag + ".rev.fastq")
      val toAlignStats = new File(sampleOutput + File.separator + sampleTag + ".stats")
      
      val toAlignUMIStats = new File(sampleOutput + File.separator + sampleTag + ".umi.stats")

      val perBaseEventFile = new File(sampleOutput + File.separator + sampleTag + ".perBase")
      val topReadFile = new File(sampleOutput + File.separator + sampleTag + ".topReadEvents")
      val topReadCount = new File(sampleOutput + File.separator + sampleTag + ".topReadCounts")
      val allReadCount = new File(sampleOutput + File.separator + sampleTag + ".allReadCounts")
      val cutSites = new File(sampleObj.reference + ".cutSites")

      // clean up the reads as much as possible
      val unpairedReads = List[File](new File(sampleOutput + File.separator + sampleTag + ".unpaired1.fq.gz"),new File(sampleOutput + File.separator + sampleTag + ".unpaired2.fq.gz"))
      val cleanedTmp = List[File](new File(sampleOutput + File.separator + sampleTag + ".cleanedInter1.fq.gz"),new File(sampleOutput + File.separator + sampleTag + ".cleanedInter2.fq.gz"))
      add(Trimmomatic(inputFiles,adaptersFile,cleanedTmp,unpairedReads))

      if (sampleObj.UMIed) {

        add(SeqPrepNoMerge(cleanedTmp, cleanedFastqs))

        // collapse the reads by their UMI and output FASTA files with the remaining quality reads
        add(UMIProcessingPaired(postProcessedFastq1,
          postProcessedFastq2,
          toAlignFastq1,
          toAlignFastq2,
          toAlignStats,
          toAlignUMIStats,
          10,
          new File(sampleObj.reference + ".primers"),
          new File(sampleObj.reference),
          sampleObj.sample,
          new File(sampleObj.reference + ".cutSites"),
          sampleObj.UMIed))

        umiStatsFiles :+= toAlignStats 

        //add(StatsAnalysis(toAlignStats, plotsDir, sampleTag))

      } else {
        add(SeqPrep(cleanedTmp, cleanedFastqs, mergedReads))

        add(ZipReadFiles(postProcessedFastq1, postProcessedFastq2, unmergedUnzipped))
        add(Gunzip(mergedReads, mergedReadUnzipped))

        add(NeedlemanAllFasta(sampleObj.reference, unmergedUnzipped, samUnmergedFasta, false))
        add(NeedlemanAllFasta(sampleObj.reference, mergedReadUnzipped, samMergedFasta, false))

        add(ReadsToStats(samUnmergedFasta,
          samMergedFasta,
          toAlignStats,
          cutSites,
          new File(sampleObj.reference + ".primers"),
          sampleObj.sample))

        nonUmiStatsFiles :+= toAlignStats 

      }

      add(ToJavascriptTables(toAlignStats, cutSites, perBaseEventFile, topReadFile, topReadCount,allReadCount))
      add(ToWebPublish(sampleWebLocation, perBaseEventFile, topReadFile, topReadCount, cutSites, allReadCount))
      add(FourGametes(toAlignStats,
        new File(sampleOutput + File.separator + sampleTag + ".fourGametes"),
        new File(sampleOutput + File.separator + sampleTag + ".alleleCounts")))
      //add(MixFile(toAlignStats,new File(sampleOutput + File.separator + sampleTag + ".mix")))

      // create the distances files
      val treeFile = new File(sampleOutput + File.separator + sampleTag + ".tree")
      val annotationFile = new File(sampleOutput + File.separator + sampleTag + ".annotations")
      val distanceFile = new File(sampleOutput + File.separator + sampleTag + ".distance")
      add(DistanceCalculation(toAlignStats, treeFile, annotationFile, distanceFile))
    })


    // UMI set finishing
    if (umiStatsFiles.size > 0) {
      
      // agg. all of the stats together, do the distance calculation, and create a tree
      val mergedStatsFile = new File(aggregateLocation + File.separator + "merged.umi.stats")
      val mergedUMIInfo = new File(aggregateLocation + File.separator + "merged.umi.info")
      add(AggregateStatsFiles(umiStatsFiles, mergedStatsFile, mergedUMIInfo, true))

      add(AggregateReport(mergedStatsFile,aggReportDir))

      val treeFile = new File(aggregateLocation + File.separator + "merged.tree")
      val annotationFile = new File(aggregateLocation + File.separator + "merged.annotations")
      val completeAnnotationFile = new File(aggregateLocation + File.separator + "merged.complete.annotations")
      val distanceFile = new File(aggregateLocation + File.separator + "merged.distance")
      add(DistanceCalculation(mergedStatsFile, treeFile, annotationFile, distanceFile))

      // create the tree files -- 
      add(CreateTree(distanceFile,treeFile,annotationFile,completeAnnotationFile))
      add(CreateHTMLTree(treeFile,aggWebTreeLocation))

    }
    // non-UMI finishing
    else {
      // agg. all of the stats together, do the distance calculation, and create a tree
      val mergedNonUMIStatsFile = new File(aggregateLocation + File.separator + "merged.noumi.stats")
      val mergedNonUMIInfo = new File(aggregateLocation + File.separator + "merged.nonumi.info")
      add(AggregateStatsFiles(nonUmiStatsFiles, mergedNonUMIStatsFile, mergedNonUMIInfo, false))

      add(AggregateReport(mergedNonUMIStatsFile,aggReportDir))
    }
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
  case class ToJavascriptTables(statsFile: File, cutSites: File, perBaseEvents: File, topReads: File, topReadCounts: File, allEventCounts: File) extends ExternalCommonArgs {
    @Input(doc = "the input stats file") var stats = statsFile
    @Input(doc = "the  cut sites information file") var cuts = cutSites
    @Output(doc = "per base event counts") var perBase = perBaseEvents
    @Output(doc = "the top reads and their base by base events") var topR = topReads
    @Output(doc = "the the counts for all the top reads") var topReadC = topReadCounts
    @Output(doc = "the the counts for all the reads") var allReadC = allEventCounts

    def commandLine = "scala -J-Xmx6g " + toJSTableScript + " " + stats + " " + perBase + " " + topR + " " + topReadC + " " + allReadC + " " + cuts

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

  // ********************************************************************************************************
  case class SamToBam(inSam: File, outBam: File, inRef: File) extends ExternalCommonArgs {
    @Input(doc = "the input sam file") var sam = inSam
    @Argument(doc = "the reference file") var ref = inRef
    @Output(doc = "output bam file") var bam = outBam

    def commandLine = samtoolsPath + " view -bS -T " + ref + " " + sam + " > " + bam

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

    def commandLine = "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/net/gs/vol1/home/aaronmck/tools/lib/; scala -J-Xmx4g -cp /net/gs/vol1/home/aaronmck/source/sandbox/aaron/lib/target/scala-2.11/sequencing-lib-assembly-1.0.jar " + smithWaterman + " " + ref + " " + fqs.mkString(" ") + " " + ref + ".primers " + outSam

    this.analysisName = queueLogDir + outSam + ".sw"
    this.jobName = queueLogDir + outSam + ".sw"
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
    def commandLine = needlePath + " -datafile /net/shendure/vol10/projects/CRISPR.lineage/nobackup/reference_data/EDNAFULL -snucleotide1 -snucleotide2 " + (if (reverse) "-sreverse2 " else " ") + "-aformat3 fasta -gapextend 0.5 -gapopen 15.0 -asequence " + ref + " -bsequence " + fq + " -outfile " + outFasta

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

  // Run the script to remove the UMI from reads that contain them
  // ********************************************************************************************************
  case class UMIProcessPaired(inFQ1: File, inFQ2: File,  outFQ1: File, outFQ2: File, outStats: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the first read fastq") var iFQ1 = inFQ1
    @Input(doc = "the second read fastq") var iFQ2 = inFQ2
    @Output(doc = "the output first-read fastq file") var oFQ1 = outFQ1
    @Output(doc = "the output second-read fastq file") var oFQ2 = outFQ2
    @Output(doc = "the output stats file") var oStats = outStats

    def commandLine = "scala " + stripUMI + " " + iFQ1 + " " + iFQ2 + " " + umiData + " " + oFQ1 + " " + oFQ2 + " " + oStats

    this.analysisName = queueLogDir + outStats + ".stripUMIPaired"
    this.jobName = queueLogDir + outStats + ".stripUMIPaired"
  }

  // take the final stats file and generate some statistics on it
  // ********************************************************************************************************
  case class AggregateReport(inputStatsFile: File, outputDirectory: File) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the input stats file") var iStats = inputStatsFile
    @Argument(doc = "the output report directory") var oDir = outputDirectory

    val outRB = new File(outputDirectory + File.separator + "report") // our base filename 

    def commandLine = "Rscript -e \"embryo.stats = read.delim('" + iStats + "'); require(knitr); require(markdown); opts_knit\\$set(base.dir = '" + oDir + "'); knit('" + aggReportScript + "', '" + outRB + ".md'); markdownToHTML('" + outRB + ".md', '" + outRB + ".html', options=c('use_xhtml', 'base64_images'));\""

    this.analysisName = queueLogDir + iStats + ".aggReport"
    this.jobName = queueLogDir + iStats + ".aggReport"
    this.memoryLimit = 25 // set this a bit high, there's a weird java / SGE interactions that require higher mem reservations
    this.residentRequest = 25
    this.residentLimit = 25
    this.isIntermediate = false // by default delete the intermediate files, if you want to keep output from a task set this to false in the case class
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

  // analyze the edits from the final alignments
  // ********************************************************************************************************
  case class PostProcess(dataFile: File, outputFile: File, workingDir: File, umiInfo: File, concInfo: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the data file") var dtable = dataFile
    @Input(doc = "the data file") var workdir = workingDir
    @Input(doc = "the uni data file") var umi = umiInfo
    @Input(doc = "the base conc. data file") var cinfo = concInfo
    @Output(doc = "the output file") var outf = outputFile

    def commandLine = "R -e \"library(knitr);library(knitrBootstrap);library(rmarkdown);conc.data='" + cinfo.getAbsolutePath + "';target.striped='" + umi.getAbsolutePath + "';table.data='" + dtable + "';opts_knit\\$set(base.dir = '" + workdir + "');render('" + cutSiteScript + "','knitrBootstrap::bootstrap_document',output_file='" + outf + "')\""

    this.analysisName = queueLogDir + outputFile + ".output"
    this.jobName = queueLogDir + outputFile + ".output"
    this.isIntermediate = false
  }

  // do the distance calculation
  // ********************************************************************************************************
  case class DistanceCalculation(umiFileInput: File, newickOuput: File, annotationOutput: File, distanceOutput: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the umi input file") var umi = umiFileInput
    @Output(doc = "newick output file") var newick = newickOuput
    @Output(doc = "annotation output file") var annotation = annotationOutput
    @Output(doc = "the distance output file") var dist = distanceOutput

    def commandLine = "java -Xmx28g -jar " + distanceMatrix + " --meltedUMIFile " + umi + " --newickFile " + newick + " --mergeInformation " + annotation + " --distanceMatrixFile " + dist

    this.analysisName = queueLogDir + newickOuput + ".distance"
    this.jobName = queueLogDir + newickOuput + ".distance"
    this.isIntermediate = false
    this.memoryLimit = 32
    this.residentRequest = 32
    this.residentLimit = 32

  }

  // aggregate stats flies down to a single file
  // ********************************************************************************************************
  case class AggregateStatsFiles(inputFiles: List[File], outputStatsFile: File, outputUmiStatsFile: File, isUMI: Boolean)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the umi input file") var inputFls = inputFiles
    @Output(doc = "the output merged stats file") var outStats = outputStatsFile
    @Output(doc = "the output umi information") var outUMIStats = outputUmiStatsFile
    @Argument(doc = "is this a UMI or not UMI data?") var umi = isUMI

    def commandLine = "scala -J-Xmx8g " + aggregateScripts + " " + inputFls.mkString(",") + " " + outStats + " " + outUMIStats + " " + (if (umi) "UMI" else "NOUMI")

    this.analysisName = queueLogDir + outStats + ".outStats"
    this.jobName = queueLogDir + outStats + ".outStats"
    this.isIntermediate = false
  }

  // create the conc. counts at each base in the first reads
  // ********************************************************************************************************
  case class ConcReadScript(reference: File, firstRead: File, output: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Output(doc = "the output file") var out = output
    @Input(doc = "the reference") var ref = reference
    @Input(doc = "the reads") var reads = firstRead

    def commandLine = "scala -J-Xmx4g -cp /net/gs/vol1/home/aaronmck/source/sandbox/aaron/lib/target/scala-2.11/sequencing-lib-assembly-1.0.jar " + concReadScript + " " + ref + " " + reads + " " + out

    this.analysisName = queueLogDir + output + ".concread"
    this.jobName = queueLogDir + output + ".concread"
    this.isIntermediate = false
  }

  // analyze each sample for recurrent events, deletions insertion rates, etc
  // ********************************************************************************************************
  case class StatsAnalysis(statsFile: File, outputDir: File, sampleName: String)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the stats file to analyze") var statsFl = statsFile
    @Input(doc = "the output directory for our plots") var outDir = outputDir
    @Argument(doc = "our sample name") var sample = sampleName

    def commandLine = "~/ipython_env/bin/python " +
      statsFileAnalysis + " --stats_input_file " + statsFl + " --output_dir " + outputDir + " --number_of_targets 10 --sample_name " + sample

    this.analysisName = queueLogDir + outDir + ".plots"
    this.jobName = queueLogDir + outDir + ".plots"
    this.isIntermediate = false
  }

  // create a tree file using the R file
  // ********************************************************************************************************
  case class CreateTree(distanceFile: File, outputTree: File, inputAnnotations: File, outputAnnotations: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the distance file") var distFile = distanceFile
    @Input(doc = "the input annotations") var inAnnotations = inputAnnotations
    @Output(doc = "the output tree") var outTree = outputTree
    @Output(doc = "the output annotations file") var outAnnotations = outputAnnotations

    def commandLine = "Rscript " + createTreeScript + " " + distFile + " " + cladeAssignments + " " + inputAnnotations + " " + outTree + " " + outAnnotations

    this.analysisName = queueLogDir + outputTree + ".rTree"
    this.jobName = queueLogDir + outputTree + ".rTree"
    this.isIntermediate = false
  }

    // ********************************************************************************************************
  case class FourGametes(statsFile: File, fourGametesOutput: File, alleleCounts: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the stats file") var stats = statsFile
    @Output(doc = "the four gametes output") var fgOutput = fourGametesOutput
    @Output(doc = "the allele counts") var alleleC = alleleCounts

    def commandLine = "scala -J-Xmx4g " + fourGameteScript + " " + stats + " " + fgOutput + " " + alleleC

    this.analysisName = queueLogDir + fgOutput + ".fgamete"
    this.jobName = queueLogDir + fgOutput + ".fgamete"
    this.isIntermediate = false
  }


    // ********************************************************************************************************
  case class MixFile(statsFile: File, mixFile: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the stats file") var stats = statsFile
    @Output(doc = "the output mix file") var mixFl = mixFile

    def commandLine = "scala -J-Xmx4g " + toMixScript + " " + stats + " " + mixFile

    this.analysisName = queueLogDir + mixFile + ".mix"
    this.jobName = queueLogDir + mixFile + ".mix"
    this.isIntermediate = false
  }


  // copy the HTML tree file over the tree location
  // ********************************************************************************************************
  case class CreateHTMLTree(inputTree: File, treeDir: File)
    extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc = "the output tree file") var iTree = inputTree
    @Output(doc = "where to put the tree file") var tDir = treeDir

    def commandLine = "scala " + treeCopyScript + " " + iTree + " " + tDir

    this.analysisName = queueLogDir + treeDir + ".setupTree"
    this.jobName = queueLogDir + treeDir + ".setupTree"
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
    appended += " SLIDINGWINDOW:5:15" // remove bases with a quality less than 20 in a sliding window of size 3
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
   * seqprep -- trimmomatic's main competition, does read merging as well
   *
   * trims off adapter sequence and low quality regions in paired end sequencing data, merges on reqest
   */
  // ********************************************************************************************************
  case class SeqPrepNoMerge(inFastqs: List[File], outs: List[File]) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input FASTAs") var fqs = inFastqs
    @Output(doc = "output fastas (corrected)") var fqOuts = outs

    if (inFastqs.length != 2 || outs.length != 2)
      throw new IllegalArgumentException("SeqPrepNoMerge can only be run on paired end sequencing! for input files " + inFastqs.mkString(", "))

    this.memoryLimit = 4
    this.residentRequest = 4
    this.residentLimit = 4

    // o -> minimum overlap, n -> faction of bases that must match in overlap
    var cmd = seqPrepPath + " -f " + fqs(0) + " -r " + fqs(1) + " -1 " + fqOuts(0) + " -2 " + fqOuts(1)

    def commandLine = cmd
    this.isIntermediate = false
    this.analysisName = queueLogDir + fqs(0) + ".SeqPrepNoMerge"
    this.jobName = queueLogDir + fqs(0) + ".SeqPrepNoMerge"
  }

  /**
   * Process the UMIs, merging and aligning reads down to a single, high-quality concensus per UMI
   * val inputFileReads = args(0)
   * val outputFile = new PrintWriter(args(1))
   * val outputStats = new PrintWriter(args(2))
   * val umiCutoff = args(3).toInt
   * val primersEachEnd = args(4)
   * val reference = args(5)
   */
  // ********************************************************************************************************
  case class UMIProcessingPaired(inMergedReads1: File, inMergedReads2: File, outputFASTA1: File,
    outputFASTA2: File, outputStats: File, outputUMIStats: File,umiCutOff: Int, primersFile: File, reference: File, sampleName: String, cutSites: File, usingUMI: Boolean) extends CommandLineFunction with ExternalCommonArgs {

    @Input(doc = "input reads (fwd)") var inReads1 = inMergedReads1
    @Input(doc = "input reads (rev)") var inReads2 = inMergedReads2
    @Input(doc = "the cutsite locations") var cutSiteFile = cutSites
    @Output(doc = "output fasta for further alignment (fwd)") var outFASTA1 = outputFASTA1
    @Output(doc = "output fasta for further alignment (rev)") var outFASTA2 = outputFASTA2
    @Output(doc = "output statistics file for containing information about the UMI merging process") var outStat = outputStats
    @Output(doc = "output statistics file of UMI merging process") var outUMIStat = outputUMIStats
    @Argument(doc = "how many UMIs do we need to initial have to consider merging them") var umiCut = umiCutOff
    @Argument(doc = "the primers file; one line per primer that we expect to have on each end of the resulting merged read") var primers = primersFile
    @Argument(doc = "the reference file") var ref = reference
    @Argument(doc = "the sample name") var sample = sampleName
    @Argument(doc = "using UMI") var umied = usingUMI

    var cmdString = "scala -J-Xmx13g /net/shendure/vol10/projects/CRISPR.lineage/nobackup/bin/UMIMerge.jar "
    cmdString += " --inputFileReads1 " + inReads1 + " --inputFileReads2 " + inReads2 + " --outputFastq1 " + outFASTA1 + " --outputFastq2 " + outFASTA2 + " --outputStats "
    cmdString += outStat + " --umiLength " + umiCut + " --primersEachEnd " + primers + " --reference " + ref + " --samplename "
    cmdString += sample + " --umiStart 0 --cutSites " + cutSiteFile + " --outputUMIStats " + outUMIStat

    var cmd = cmdString

    this.memoryLimit = 16
    this.residentRequest = 16
    this.residentLimit = 16

    def commandLine = cmd
    this.isIntermediate = false
    this.analysisName = queueLogDir + outStat + ".umis"
    this.jobName = queueLogDir + outStat + ".umis"
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
    var cmd = seqPrepPath + " -f " + fqs(0) + " -r " + fqs(1) + " -1 " + fqOuts(0) + " -2 " + fqOuts(1) + " -s " + merged

    def commandLine = cmd
    this.isIntermediate = false
      this.analysisName = queueLogDir + fqs(0) + ".seqprep"
    this.jobName = queueLogDir + fqs(0) + ".seqprep"
  }
}
