
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
  * make trees for samples
  * 
  **/
class DNAQC extends QScript {
  qscript =>

  /** **************************************************************************
    * Data Parameters
    * ************************************************************************** */
  @Input(doc = "Tab-delimited tear sheet containing sample information", fullName = "input", shortName = "i", required = true)
  var input: File = _

  /** **************************************************************************
    * Path parameters -- where to find tools
    * ************************************************************************** */
  @Input(doc = "The path to a functional installation of the scala tool", fullName = "scala", shortName = "scala", required = false)
  var scalaPath: File = new File("/net/gs/vol1/home/aaronmck/tools/bin/scala")

  @Input(doc = "convert a stats file into a mix ready file", fullName = "tomix", shortName = "tomix", required = false)
  var toMixScript: File = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/stats_to_mix_input_embryos.scala")

  @Input(doc = "run's the mix program", fullName = "runmix", shortName = "runmix", required = false)
  var runMixScript: File = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/control_mix_program.scala")

  @Input(doc = "generate the tree from the mix output and annotations", fullName = "totree", shortName = "totree", required = false)
  var generateTree: File = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/TreeUtils/target/scala-2.11/TreeUtils-assembly-1.1.jar")

  @Argument(doc = "where we drop out final tree files", fullName = "treedir", shortName = "treedir", required = false)
  var treedir: File = new File("/net/shendure/vol2/www/content/members/aaron/fate_map/all_trees/tree_data/")

  /** **************************************************************************
    * Global Variables
    * ************************************************************************** */

  // Gracefully hide Queue's output -- this is very important otherwise there's Queue output everywhere and it's a mess
  val queueLogDir: String = ".qlog/"

  // ** **************************************************************************
  // * Main script entry point
  // * **************************************************************************
  def script() {
    var statsFiles = List[File]()

    // read in the tear sheet and process each sample
    parseTearSheet(input).foreach(sampleObj => {

      val mixConsumingFile = sampleObj.outputDir + "/" + sampleObj.name + ".mix"
      val annotationsFile = sampleObj.outputDir + "/" + sampleObj.name + ".annotations"
      val samplesToClade = sampleObj.outputDir + "/" + sampleObj.name + ".sampleToClade"
      val weightsFile = sampleObj.outputDir + "/" + sampleObj.name + ".weights"
      val eventsToNumbers = sampleObj.outputDir + "/" + sampleObj.name + ".eventsToNumber"
      val filterNumber = sampleObj.outputDir + "/" + sampleObj.name + ".filterNumber"
      val mixFile = sampleObj.outputDir + "/outfile"
      val mixTree = sampleObj.outputDir + "/outtree"
      val outputTree =  treedir + "/" + sampleObj.name + ".json"

      add(CreateMixFiles(sampleObj.allReadsFile, mixConsumingFile, annotationsFile, weightsFile, eventsToNumbers, sampleObj.name, samplesToClade, filterNumber))

      add(RunMix(mixConsumingFile, weightsFile, mixTree, mixFile, sampleObj.outputDir))

      add(CreateTree(mixTree, mixFile, annotationsFile, samplesToClade, eventsToNumbers, outputTree))

    })
  }

  case class SourceEntry(name: File, allReadsFile: File, outputDir: File)

  /**
    * reads in the sample tear sheet, a tab separated file with the columns listed at the top of the file
    */
  def parseTearSheet(input: File): Array[SourceEntry] = {
    val inputLines = Source.fromFile(input).getLines

    // check that the header contains the correct information
    if (inputLines.next().stripLineEnd != "sample\tallReadsFile\tworkingDir")
      throw new IllegalArgumentException("Your header doesn't seem like a correctly formatted tear sheet!")

    return inputLines.map(line => {
      val tokens = line.split( """\s+""") // Use a regex to avoid empty tokens
      (new SourceEntry(
        tokens(0), // sample
        new File(tokens(1)), // all reads file
        new File(tokens(2)) // output dir file
        ))
    }).toArray
  }

  /** **************************************************************************
    * traits that get tacked onto runnable objects
    * ************************************************************************** */

  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 6 // set this a bit high, there's a weird java / SGE interactions that require higher mem reservations
    this.residentRequest = 6
    this.residentLimit = 6
    this.isIntermediate = false // by default delete the intermediate files, if you want to keep output from a task set this to false in the case class
  }

  /** **************************************************************************
    *  helper functions
    * ************************************************************************** */

  // ********************************************************************************************************
  case class CreateMixFiles(allReadsFile: File, mixOutput: File, annotations: File, weights: File, eventsToNumbers: File, sampleName: String, sampleToClade: File, filteringNumber: File) extends ExternalCommonArgs {
    @Input(doc = "the all reads file input") var allR = allReadsFile
    @Output(doc = "mix output file") var mixF = mixOutput
    @Output(doc = "the annotations output file") var annot = annotations
    @Output(doc = "the weights for annotations") var wghts = weights
    @Output(doc = "the events to numbers file") var eventsToNum = eventsToNumbers
    @Output(doc = "filtering output file") var filterNum = filteringNumber
    @Output(doc = "sample to clade file") var sample2c = sampleToClade
    @Argument(doc = "our sample name to output") var sample = sampleName

    def commandLine = scalaPath + " -J-Xmx6g " + toMixScript + " " + allR + " " + mixF + " " + annot + " " + wghts + " " + eventsToNum + " " + sample + " " + sample2c + " " + filterNum

    this.analysisName = queueLogDir + allR + ".mixFiles"
    this.jobName = queueLogDir + allR + ".mixFiles"
    this.memoryLimit = 8
    this.residentRequest = 8
    this.residentLimit = 8
  }

  // ********************************************************************************************************
  case class RunMix(mixFile: File, weightsFile: File, fakeTreeOutput: File, fakeOutputFile: File, outputDir: File) extends ExternalCommonArgs {
    @Input(doc = "the mix file") var mxFile = mixFile
    @Input(doc = "the weights file") var wgts = weightsFile
    @Output(doc = "fake tree") var ftree = fakeTreeOutput
    @Output(doc = "fake output file") var foutput = fakeOutputFile

    def commandLine = scalaPath + " -J-Xmx2g " + runMixScript + " " + mxFile + " " + wgts + " " + outputDir

    this.analysisName = queueLogDir + mxFile + ".mixed"
    this.jobName = queueLogDir + mxFile + ".mixed"
    this.memoryLimit = 1
    this.residentRequest = 1
    this.residentLimit = 1
  }

  // ********************************************************************************************************
  case class CreateTree(mixTree: File, mixGenotypes: File, annotations: File, sampleToClade: File, eventToNumber: File, outputTree: File) extends ExternalCommonArgs {
    @Input(doc = "the newick tree from mix") var tree = mixTree
    @Input(doc = "the genotypes from mix") var geno = mixGenotypes
    @Input(doc = "the annotations for each allele") var annot = annotations
    @Input(doc = "the mapping of sample to clade") var stoc = sampleToClade
    @Input(doc = "the allele event to the mix assigned number") var eventToNum = eventToNumber
    @Output(doc = "the output tree we make") var outTree = outputTree

    def commandLine = "java -Xmx5g -jar " + generateTree + " --inputTree " + tree + " --inputGenotypes " + geno + " --inputAnnotations " + annot + " --inputSampleToClade " + sampleToClade +
      " --inputEventsToNumbers " + eventToNumber + " --outputTree " + outTree

    this.analysisName = queueLogDir + outTree + ".tree"
    this.jobName = queueLogDir + outTree + ".tree"
    this.memoryLimit = 6
    this.residentRequest = 6
    this.residentLimit = 6
  }
}
