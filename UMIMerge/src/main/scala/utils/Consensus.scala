package main.scala

import aligner.Clustalo
import main.scala.utils.Utils

import scala.collection.mutable._
import scala.main._

/**
 * Functional transforms that revolve around merging reads down to consensus, filtering out poor hits, etc
 */
object Consensus {

  /**
   * develop a consensus over the read -- look at each aligned position and find the most common character and it's proportion
   * return the consensus read and proportion of bases at each position that match the called base
   *
   * @param reads the reads to make a consensus over
   * @return a tuple, representing the concensus read and the per-base error rate over positions in that read
   */
  def consensus(reads: Array[SequencingRead]): SequencingRead = {
    // each string should be the same length coming back from CLUSTAL
    var highestOccuringBase = Array[Char]()
    var highestOccuringBaseProp = Array[Double]()

    // get our read length, if we have any reads
    val readLength = if (reads.length > 0) reads(0).bases.length else 0

    (0 until readLength).foreach { index => {
      var bases = new HashMap[Char, Int]()

      // add reads base at the current position
      reads.foreach { read => bases(read.bases(index)) = bases.getOrElse(read.bases(index), 0) + 1}

      var maxChar = 'U'
      var total = 0
      bases.foreach { case (base, count) => {
        if (maxChar == 'U' || count > bases(maxChar)) maxChar = base
        if (count == bases(maxChar) && maxChar == '-') maxChar = base // favor a base over a dash
        total += count
      }
      }

      highestOccuringBase :+= maxChar
      highestOccuringBaseProp :+= (if (bases contains maxChar) bases(maxChar).toDouble / total.toDouble else 0.0)
    }
    }

    SequencingRead("Consensus",
      highestOccuringBase.mkString(""),
      "H" * highestOccuringBase.length, // highestOccuringBaseProp.map{pr => Utils.probabilityOfErrorToPhredChar(1.0 - pr)}.mkString(""),
      if (reads.size == 0) ForwardReadOrientation else reads(0).readOrientation,
      "")
  }

  /**
   * do a first pass consensus, then remove any reads that are X bases (proportionally) different
   * than the consensus read, then make a new consensus returns a tuple-4: the consensus, it's average
   * per-base error rate, number of reads used in the final assembly, and the proportion retained of the total
   *
   * @param reads the reads, with appropriate forward and reverse correctly annotationed
   * @param differenceThreshold the amount of divergence we allow
   * @return
   */
  def removeMismatchedReads(reads: Array[SequencingRead], differenceThreshold: Double = 0.9, minReadLength: Int = 50): Array[SequencingRead] = {

    val consensusValue = consensus(reads)

    var newReads = Array[SequencingRead]()
    reads.foreach { read => {
      val differences = consensusValue.bases.zip(read.bases).map { case (concBase, readBase) => if (readBase != concBase) 1 else 0 }.sum.toDouble / read.bases.length.toDouble
      if (read.bases.length > minReadLength && differences <= differenceThreshold)
        newReads :+= read
    }}
    newReads
  }

  /**
   * split a pile of reads out into a Tuple3 of fwd, rev, and reference reads
   * @param reads a mixed array of reads
   * @return a tuple3 of forward reads, reverse reads, and reference reads
   */
  def splitReadsToFwdRevRefArrays(reads: Array[SequencingRead], throwExceptionIfNoRef: Boolean = false): (Array[SequencingRead], Array[SequencingRead], Array[SequencingRead]) = {
    // split the reads into forward and reverse reads, and develop a consensus from each
    val fwdReads = reads.filter { rd => rd.readOrientation == ForwardReadOrientation }.toArray
    val revReads = reads.filter { rd => rd.readOrientation == ReverseReadOrientation }.toArray
    val refRead = reads.filter { rd => rd.readOrientation == ReferenceRead }.toArray

    if (throwExceptionIfNoRef && refRead.size != 1)
      throw new IllegalStateException("Unable to find one and only reference sequence in the reads, instead we saw " + refRead.size)
    (fwdReads, revReads, refRead)
  }

  /**
   * given forward and backwards reads (throw in a reference read to be sure), call cigar events from it
   * @param reads the reads to develop a consensus from, including a reference read
   * @return a read with the dashes split-out
   */
  def dualConsensus(reads: Array[SequencingRead], minReadLength: Int, minMeanQualScore: Double): Array[SequencingRead] = {

    val (fwdReads: Array[SequencingRead], revReads: Array[SequencingRead], refRead: Array[SequencingRead]) = splitReadsToFwdRevRefArrays(reads)

    // get the forward and reverse consensus for each read
    Array[SequencingRead](consensus(fwdReads), consensus(revReads))
  }

  /**
   * given forward and backwards reads (throw in a reference read to be sure), call cigar events from it
   * @param reads the reads to develop a consensus from, including a reference read
   * @return a read with the dashes split-out

  def callEvents(reads: Array[SequencingRead], minReadLength: Int, minMeanQualScore: Double): (SequencingRead, SequencingRead, Array[CigarHit]) = {

    val (fwdReads: Array[SequencingRead], revReads: Array[SequencingRead], refRead: Array[SequencingRead]) = splitReadsToFwdRevRefArrays(reads)

    // get the forward and reverse consensus for each read
    val fwdConcensus = consensus(refRead ++ fwdReads)
    val revConcensus = consensus(refRead ++ revReads)

    // merge the two consensus reads with the reference, and call events
    val readsToMerge: Array[SequencingRead] = Array[SequencingRead](fwdConcensus, revConcensus)
    val merging = Clustalo(Consensus.prepareReferenceConsensus(readsToMerge, refRead(0), minReadLength, minMeanQualScore))

    // now call CIGAR strings using the reads
    val fwdCigar = fwdConcensus.callCigarsAgainstOtherRead(refRead(0))
    val revCigar = revConcensus.callCigarsAgainstOtherRead(refRead(0))

    // merge the two call sets down to a common set of events spanning the whole reference
    //val

    return (fwdConcensus, revConcensus, fwdCigar ++ revCigar)
  } */

  /**
   * given the reference sequence as the first read in a multiread alignment, strip off columns before the reference's
   * first base (read into the adapter), after the reference last base (into the other adapter)
   * @param reads the reads, of which the reference read is the first
   * @return an array of reads with adapter and other non-sense stripped off the ends
   */
  def referenceStrip(reads: Array[SequencingRead]): Array[SequencingRead] = {
    if (reads.length < 3)
      throw new IllegalStateException("Unable to process reads when theres less than 3 reads")

    val referenceRead = reads(0)

    // find the first non-dash base in the reference
    val firstAndLastColumns = referenceRead.firstAndLastActualBases()

    reads.map { rd => rd.slice(firstAndLastColumns._1, firstAndLastColumns._2) }.toArray
  }


  /**
   * do the consensus matching with the reference string
   *
   * @param fwdReads forward reads as an array of SequencingRead
   * @param revReads reverse reads
   * @param reference the reference as a string, which we make into a SequencingRead
   * @return an array of sequencing reads with the reference in it
   *
   */
  def prepareReferenceConsensus(fwdReads: Array[SequencingRead],
                                revReads: Array[SequencingRead],
                                reference: String,
                                minimumReadLength: Int,
                                minimumQualScore: Double): Array[SequencingRead] = {

    // reverse-complement each of the reads in the reverse orientation
    val reverseRead = revReads.map { read => SequencingRead.reverseComplement(read) }

    // make a read of the reference
    val refQuals = ("H" * reference.length).mkString("")
    val referenceRead = SequencingRead("reference", reference, refQuals, ReferenceRead, "UNKNOWN")

    // check that our reads meet the minimum thresholds to be merged
    val (fwdFiltered, revFiltered) = checkReadsMetMinimumInputQuality(fwdReads, reverseRead, minimumReadLength, minimumQualScore)

    // create a consensus of reference, forward, and reverse reads
    Array[SequencingRead](referenceRead) ++ fwdFiltered ++ revFiltered
  }

  /**
   * do the consensus matching with the reference string
   *
   * @param reads all reads as an array of SequencingRead
   * @param referenceRead the reference read
   * @return an array of sequencing reads with the reference in it
   */
  def prepareReferenceConsensus(reads: Array[SequencingRead],
                                referenceRead: SequencingRead,
                                minimumReadLength: Int,
                                minimumQualScore: Double): Array[SequencingRead] = {
    // check that our reads meet the minimum thresholds to be merged
    val filteredReads = checkReadsMetMinimumInputQuality(reads, minimumReadLength, minimumQualScore)

    // create a consensus of reference, forward, and reverse reads
    Array[SequencingRead](referenceRead) ++ filteredReads
  }

  /**
   * do the consensus matching with the reference string
   *
   * @param reads all reads as an array of SequencingRead
   * @return an array of sequencing reads with the reference in it
   */
    def prepareConsensus(reads: Array[SequencingRead],
                         minimumReadLength: Int,
                         minimumQualScore: Double): Array[SequencingRead] = {
    // check that our reads meet the minimum thresholds to be merged
    val filteredReads =   checkReadsMetMinimumInputQuality(reads, minimumReadLength, minimumQualScore)

    // create a consensus of reference, forward, and reverse reads
    filteredReads
  }

  /**
   * do the consensus matching with the reference string
   *
   * @param reads all reads as an array of SequencingRead
   * @param reference the reference as a string, which we make into a SequencingRead
   * @return an array of sequencing reads with the reference in it
   */
  def prepareReferenceConsensus(reads: Array[SequencingRead],
                                reference: String,
                                minimumReadLength: Int,
                                minimumQualScore: Double): Array[SequencingRead] = {
    // make a read of the reference
    val refQuals = ("H" * reference.length).mkString("")
    val referenceRead = SequencingRead("reference", reference, refQuals, ReferenceRead, "UNKNOWN")

    // check that our reads meet the minimum thresholds to be merged
    val filteredReads = checkReadsMetMinimumInputQuality(reads, minimumReadLength, minimumQualScore)

    // create a consensus of reference, forward, and reverse reads
    Array[SequencingRead](referenceRead) ++ filteredReads
  }

  /**
   * a chance to filter out any reads that may challenge us latter: too short, really poor quality, etc.  We do this on the
   * whole set at once to keep the forward and reverse reads in sync., otherwise we risk merging unmatched reads later on
   *
   * @param fwdReads the forward reads
   * @param revReads the reverse reads
   * @return a tuple of the matched forward and reverse reads
   */
  def checkReadsMetMinimumInputQuality(fwdReads: Array[SequencingRead],
                                       revReads: Array[SequencingRead],
                                       minimumReadLength: Int,
                                       minimumMeanQualScore: Double): Tuple2[Array[SequencingRead], Array[SequencingRead]] = {
    var retFwd = new ArrayBuffer[SequencingRead]()
    var retRev = new ArrayBuffer[SequencingRead]()

    fwdReads.zip(revReads).foreach { case (fwd, rev) => {
      var valid = true
      if (fwd.length < minimumReadLength || rev.length < minimumReadLength)
        valid = false
      if (fwd.intQuals.sum.toDouble / fwd.length.toDouble < minimumMeanQualScore || rev.intQuals.sum.toDouble / rev.length.toDouble < minimumMeanQualScore)
        valid = false

      if (valid) {
        retFwd += fwd
        retRev += rev
      }
    }
    }
    return ((retFwd.toArray, retRev.toArray))
  }

  /**
   * a chance to filter out any reads that may challenge us latter: too short, really poor quality, etc.  We do this on the
   * whole set at once to keep the forward and reverse reads in sync., otherwise we risk merging unmatched reads later on
   *
   * @param reads the reverse reads
   * @return a tuple of the matched forward and reverse reads
   */
  def checkReadsMetMinimumInputQuality(reads: Array[SequencingRead],
                                       minimumReadLength: Int,
                                       minimumMeanQualScore: Double): Array[SequencingRead] = {
    var retReads = new ArrayBuffer[SequencingRead]()

    reads.foreach { case (read) => {
      val quals = read.intQuals.zipWithIndex.filter { case (qual, index) => read.quals(index) != '-' }.map { case (q, i) => q }.toArray
      val qualSum = quals.sum.toDouble
      var valid = true
      if (read.length < minimumReadLength) {
        //println("SHORT: Dropping read " + read.umi + " quals " + read.intQuals.mkString(",") + " with min qual " + (qualSum / quals.length.toDouble) + " old len " + read.intQuals.length)
        valid = false
      }

      if (qualSum / quals.length.toDouble < minimumMeanQualScore) {
        //println("QUAL: Dropping read " + read.umi + " quals " + read.intQuals.mkString(",") + " with min qual " + (qualSum / quals.length.toDouble) + " old len " + read.intQuals.length)
        valid = false
      }

      if (valid) {
        retReads += read
      }
    }
    }
    val retReadsArray = retReads.toArray
    if (retReadsArray.size == 0)
      println("WARNING: we filtered " + reads.size + " down to " + retReadsArray.size)


    return (retReads.toArray)
  }

}
