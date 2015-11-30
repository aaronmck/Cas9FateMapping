package test.scala

import aligner.MAFFT2
import main.scala.Downsampler
import utils.CutSites

import org.scalatest.{Matchers, FlatSpec}

import scala.main.{ForwardReadOrientation, SequencingRead}

/**
 * Created by aaronmck on 11/29/15.
 */
class DownsamplerTest extends FlatSpec with Matchers {

  "A DownsamplerTest" should "find the right subset correctly" in {
    val readSeq =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq2 =  "TTTTTTTT--------------TTTTTTTT"
    val readQual =  "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    val readQual2 = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
    val testRead = SequencingRead("FAKE1", readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead("FAKE2", readSeq2, readQual2, ForwardReadOrientation, "TTTT")
    val testRead3 = SequencingRead("FAKE3", readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val ret = Downsampler.downsample(Array[SequencingRead](testRead,testRead2,testRead3,testRead3,testRead3), 2)
    ret.size should be (2)
    ret(0).averageQual() should be (40.0)
    ret(1).averageQual() should be (40.0)
  }

  "A DownsamplerTest2" should "should find one indel right" in {
    val readSeq =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq2 =  "TTTTTTTT--------------TTTTTTTT"
    val readQual =  "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    val readQual2 = "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"
    val testRead = SequencingRead("FAKE1", readSeq, readQual, ForwardReadOrientation, "TTTT")

    val ret = Downsampler.downsample(Array[SequencingRead](testRead), 1)
    ret.size should be (1)
    ret(0).averageQual() should be (33.0)
  }
}
