package aligner

import scala.main.{ConsensusRead, SequencingRead}

/**
 * Created by aaronmck on 11/17/15.
 */
trait Aligner {
  def alignTo(reads: Array[SequencingRead], ref: Option[String], debug: Boolean = false): Array[SequencingRead]


}

object Aligner {
  /**
   * create a sequencing read from a name a base string, just to save time in the alignment return
   * @param name the name of the read
   * @param bases the bases of the read
   * @return a sequencing read
   */
  def SequencingReadFromNameBases(name: String, bases: String): SequencingRead = {
    return SequencingRead(name, bases, "H" * bases.length, ConsensusRead, "UNKNOWN")
  }
}