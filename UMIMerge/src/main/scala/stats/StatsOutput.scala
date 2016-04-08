package main.scala.stats

import java.io._

/**
 * this class handles outputting a statistics file for runs of the umi/non-umi amplicon sequencing
 */
class StatsOutput(outputFile: File, targetCount: Int, outputReadStrings: Boolean = true) {

  // setup the output file
  val output = new PrintWriter(outputFile.getAbsolutePath)

  // output the header
  output.write(StatsOutput.standardHeaderItems(outputReadStrings) +
    StatsOutput.sepString +
    StatsOutput.createTargetString(targetCount) +
    StatsOutput.endline)

  /**
   * output a stats entry object to the stats file
    *
    * @param container a filled out container with read statistics
   */
  def outputStatEntry(container: StatsContainer): Unit = {
    output.write(container.outputString(outputReadStrings) + StatsOutput.endline)
  }

  def close() = output.close()
}


/**
 * our static methods for generating non-stateful aspects of the stats file
 */
object StatsOutput {
  def sepString = "\t"

  def endline = "\n"

  // return the standard header items that all stats files have
  val standardHeaderItemsArray = Array[String]("readName", "keep", "conflict", "hasForwardPrimer", "hasReversePrimer"
    , "umi", "merged", "mergedReadLen", "read1len",
    "read2len", "finalReadCount1", "finalReadCount2", "matchRate1",
    "matchRate2", "alignedBases1", "alignedBases2")

  // if we're outputing read alignment strings we include these entries
  val standardHeaderItemsArrayAddition = Array[String]("fwdRead","revRead","mergedRead","fwdReadRef","revReadRef","mergedReadRef")

  // get the output header strings
  def standardHeaderItems(outputReadStrings: Boolean) = {
    if (outputReadStrings) {
      (standardHeaderItemsArray ++ standardHeaderItemsArrayAddition).mkString(sepString)
    } else {
      (standardHeaderItemsArray).mkString(sepString)
    }
  }

  // create the target specific strings at the top of the file
  def createTargetString(numberOfTargets: Int): String = {
    (1 until (numberOfTargets + 1)).map { i => "target" + i }.toArray[String] ++
      (1 until (numberOfTargets + 1)).map { i => "sequence" + i }.toArray[String]
  }.mkString(sepString)
}

// case class that must be filled out to output data into a stats file.  It enforces type checking
// which is nice, and formats the strings for output
case class StatsContainer(name: String, keep: Boolean, isConflicted: Boolean, hasFWDPrimer: Boolean, hasREVPrimer: Boolean, isUMI: Boolean,
                          isMerged: Boolean, read1Len: Int, read2Len: Int, finalRead1Count: Int,
                          finalRead2Count: Int, matchRate1: Double, matchRate2: Double, alignedBases1: Int,
                          alignedBases2: Int, targetEvents: Array[String], targetSequences: Array[String],
                          alignedFWDRead: Option[String], alignedRVSRead: Option[String], mergedRead: Option[String],
                          alignedFWDReadRef: Option[String], alignedRVSReadRef: Option[String], mergedReadRef: Option[String]) {

  // check that the target events and target sequence counts are the same
  if (targetEvents.size != targetSequences.size)
    throw new IllegalStateException("The target event array and target sequence array are not the same length: " + targetEvents.size + " and " + targetSequences.size)

  // create an output string version of the target data
  val targetEventString = targetEvents.mkString(StatsOutput.sepString)
  val targetSequenceString = targetSequences.mkString(StatsOutput.sepString)

  // get rid of bad characters in any read name
  def convertName(readName: String): String = readName.replace(' ', '_')

  // convert the keep boolean to string
  def convertKeep(keepB: Boolean): String = if (keepB) "PASS" else "FAIL"

  // convert the merged boolean to string
  def convertMerged(mergedB: Boolean): String = if (mergedB) "MERGED" else "PAIR"

  // convert the umi boolean tag to a string
  def convertUMI(umiB: Boolean): String = if (umiB) "UMI" else "AMPLICON"

  // convert the umi boolean tag to a string
  def convertConflicted(conflicted: Boolean): String = if (conflicted) "CONFLICTED" else "CONSISTENT"

  def createReadStrings(): String = {
    val rdFwd = alignedFWDRead.getOrElse("NA")
    val rdRev = alignedRVSRead.getOrElse("NA")
    val rdMrg = mergedRead.getOrElse("NA")
    val rdFwdRef = alignedFWDReadRef.getOrElse("NA")
    val rdRevRef = alignedRVSReadRef.getOrElse("NA")
    val rdMrgRef = mergedReadRef.getOrElse("NA")
    rdFwd + StatsOutput.sepString + rdRev + StatsOutput.sepString + rdMrg + StatsOutput.sepString +
      rdFwdRef + StatsOutput.sepString + rdRevRef + StatsOutput.sepString + rdMrgRef
  }

  // our final output string
  def outputString(outputReadStrings: Boolean) =
    convertName(name) + StatsOutput.sepString +
      convertKeep(keep) + StatsOutput.sepString +
      convertConflicted(isConflicted) + StatsOutput.sepString +
      hasFWDPrimer + StatsOutput.sepString +
      hasREVPrimer + StatsOutput.sepString +
      convertUMI(isUMI) + StatsOutput.sepString +
      convertMerged(isMerged) + StatsOutput.sepString +
      read1Len + StatsOutput.sepString +
      read2Len + StatsOutput.sepString +
      finalRead1Count + StatsOutput.sepString +
      finalRead2Count + StatsOutput.sepString +
      finalRead2Count + StatsOutput.sepString +
      matchRate1 + StatsOutput.sepString +
      matchRate2 + StatsOutput.sepString +
      alignedBases1 + StatsOutput.sepString +
      alignedBases2 +
      (if (outputReadStrings) {StatsOutput.sepString + createReadStrings} else {""}) + StatsOutput.sepString +
      targetEventString + StatsOutput.sepString + targetSequenceString
}