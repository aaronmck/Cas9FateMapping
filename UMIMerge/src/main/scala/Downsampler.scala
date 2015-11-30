package main.scala

import scala.collection.mutable.HashMap
import scala.main.{SequencingReadQualOrder, SequencingRead}
import scala.util.Sorting

/**
 * Created by aaronmck on 11/29/15.
 */
object Downsampler {
  /**
   * downsample
   * @param reads
   * @param returnSize
   * @return
   */
  def downsample(reads: Array[SequencingRead], returnSize: Int): Array[SequencingRead] = {
    val newReads = reads

    Sorting.quickSort(newReads)(SequencingReadQualOrder)

    return newReads.slice(reads.size - returnSize,reads.size)
  }

}
