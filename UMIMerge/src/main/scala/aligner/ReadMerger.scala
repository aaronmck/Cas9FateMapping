package main.scala.aligner

import aligner.Aligner
import aligner.aligner.Waterman

import scala.main.SequencingRead

/**
 * rely on needle to do the alignment of the two reads, and post-process for successful alignment
 */
case class ReadMergerRet(read: SequencingRead, matches: Int, overlap: Int, orientationOK: Boolean)

object ReadMerger {
  def mergeRead(read1: SequencingRead, read2: SequencingRead): Option[ReadMergerRet] = {
    val alignments = Waterman.alignTo(Array[SequencingRead](read1, read2), None, false, 100.0, 0.1)
    if (alignments.size != 2)
      return None

    var overlap = 0
    var matches = 0
    var output = Array[Char]()

    // sanity check that the alignment put read1 first, then some merged, then some read2
    var situationFWDOK = false
    var situationREVOK = false
    var situationFWDSet = false
    var situationREVSet = false

    alignments(0).bases.zip(alignments(1).bases).foreach { case (base1, base2) => (base1, base2) match {
      case (b1, b2) if b1 == '-' && b2 == '-' => {/* keep going */}
      case ('-',b2) => {if (!situationFWDSet) {situationFWDSet = true; situationFWDOK = false; }}
      case (b1, '-') => {if (!situationFWDSet) {situationFWDSet = true; situationFWDOK = true}}
      case _ => {}
    }
    }
    alignments(0).bases.zip(alignments(1).bases).reverse.foreach { case (base1, base2) => {
      (base1, base2) match {
        case (b1, b2) if b1 == '-' && b2 == '-' => {/* keep going */}
        case ('-',b2) => {if (!situationREVSet) {situationREVSet = true; situationREVOK = true; }}
        case (b1, '-') => {if (!situationREVSet) {situationREVSet = true; situationREVOK = false; }}
        case _ => {}
      }
    }}


    alignments(0).bases.zip(alignments(1).bases).foreach { case (base1, base2) => (base1, base2) match {
      case (b1, b2) if b1 == '-' && b2 != '-' => output :+= b2
      case (b1, b2) if b1 != '-' && b2 == '-' => output :+= b1
      case (b1, b2) if b1 != '-' && b2 != '-' && b1 == b2 => {
        output :+= b1
        overlap += 1
        matches += 1
      }
      case (b1, b2) if b1 != '-' && b2 != '-' && b1 != b2 => {
        output :+= b1
        overlap += 1
      }
    }
    }

    println(read1.umi + "\t" + overlap + "\t" + matches + "\t" + alignments(0).bases + "\t" + alignments(1).bases + "\t" + situationFWDOK + " " + situationREVOK + " " + situationFWDSet + " " + situationREVSet)

    val returnRead = Aligner.SequencingReadFromNameBases("Consensus", output.filter { b => b != '-' }.mkString(""))
    return Some(ReadMergerRet(returnRead, matches, overlap, situationFWDOK & situationREVOK & situationFWDSet & situationREVSet))
  }
}
