package main.scala.aligner

import aligner.Aligner
import aligner.aligner.Waterman

import scala.main.SequencingRead

/**
 * Created by aaronmck on 12/2/15.
 */
object ReadMerger {
  def mergeRead(read1: SequencingRead, read2: SequencingRead): Option[Tuple3[SequencingRead,Int,Int]] = {
    val alignments = Waterman.alignTo(Array[SequencingRead](read1,read2),None,false,6.0,1.0)
    if (alignments.size != 2)
      return None

    var overlap = 0
    var matches = 0
    var output = Array[Char]()

    alignments(0).bases.zip(alignments(1).bases).foreach{case(base1,base2) => (base1,base2) match {
      case(b1,b2) if b1 == '-' && b2 != '-' => output :+= b2
      case(b1,b2) if b1 != '-' && b2 == '-' => output :+= b1
      case(b1,b2) if b1 != '-' && b2 != '-' && b1 == b2 => {
        output :+= b1
        overlap += 1
        matches += 1
      }
      case(b1,b2) if b1 != '-' && b2 != '-' && b1 != b2 => {
        output :+= b1
        overlap += 1
      }
    }}


    val returnRead = Aligner.SequencingReadFromNameBases("Consensus", output.filter{b => b!='-'}.mkString(""))
    return Some(returnRead,matches,overlap)
  }
}
