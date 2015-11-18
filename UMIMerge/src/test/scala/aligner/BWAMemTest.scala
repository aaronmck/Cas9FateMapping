package aligner

import java.io.File

import org.scalatest.{Matchers, FlatSpec}

import scala.main.{ForwardReadOrientation, SequencingRead}

/**
 * Created by aaronmck on 11/17/15.
 */
class BWAMemTest extends FlatSpec with Matchers {
  val reference = new File("./test_files/ref.fa")

  "BWA" should "align correctly" in {
    val read1 = SequencingRead("read1",
      "ATCTCGAGCTCAAGCTTCGGATACGATACGCGCACGCTATGGAGTCGACACGACTCGCGC",
      "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH",
      ForwardReadOrientation,
      "GGGGGGGGGG")

    val bwa = MAFFT(Array[SequencingRead](read1), reference)
    bwa.returnReads.foreach(read => {
      println("read 1 " + read.name)
    })
  }
}