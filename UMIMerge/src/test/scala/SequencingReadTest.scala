import main.scala.utils.Utils
import org.scalatest.junit.JUnitSuite

import scala.main.{ForwardReadOrientation, SequencingRead}
import org.scalatest._

/**
 * Created by aaronmck on 10/22/15.
 */
class SequencingReadTest extends FlatSpec with Matchers {
  val readName = "TestRead1"

  "A SequencingRead" should "find its non-dash length correctly" in {
    val readSeq = "AAAAAAAA-------"
    val readQual = "HHHHHHHHHHHHHHH"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    testRead.trueEnd() should be (7)
  }

  "A SequencingRead" should "not report a smaller subset when there isnt one" in {
    val readSeq = "AAAAAAAAAAAAAAA"
    val readQual = "HHHHHHHHHHHHHHH"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    testRead.trueEnd() should be (readSeq.length -1)
  }

  "A SequencingRead" should "report a full dash read appropriately as zero" in {
    val readSeq =  "---------------"
    val readQual = "HHHHHHHHHHHHHHH"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    testRead.trueEnd() should be (-1)
  }

  "A SequencingRead" should "slice down a read to the right size" in {
    val readSeq =  "---------------"
    val readQual = "HHHHHHHHHHHHHHH"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT").slice(0,5)
    testRead.length should be (5)
  }

  "A SequencingRead" should "should quality threshold itself down when qual drops off in a one base window" in {
    val readSeq =  "AAAAAAAAAAAAAAA"
    val readQual = "HHHHHHHHH$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT").qualityThresholdRead(1,10)
    testRead.bases should be ("AAAAAAAAA")
  }

  "A SequencingRead" should "should quality threshold itself down when qual drops off in a two base window" in {
    val readSeq =  "123456789123456"
    val readQual = "HHHHHHHHH$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT").qualityThresholdRead(2,10)
    testRead.bases should be ("123456789")
  }

  "A SequencingRead" should "should trim some real garbage" in {
    // ...
    val readSeq =  "CAATTTCCCTACTCAGATCTCGACCTCGAGACAGGTTTGGAGCGAGATTGATAAAGTGCCGACGGTCACGCTTGCGATCTCGTATGCCGTCTTCTGCTTGTAAAAAAAAACTATCTCTCCCTCTACTTATTATTCTTTATATTTTATCTCCCCCTTCACTTCACCAACATTTATCTACTAACTTAACTCAAAAGCTTTTTAATATATATTCTTAAACAAACATTATCTATTCATCTCACAATTCAACACAAGTTTA"
    val readQual = "CCCCCFFFFGGGGGGGFFFGCG,CCGFG7@FDGGGEGGGGGGGGGFGDFE<FFFFFGFFGGGGGGGGEGGGGFEFCCF<FGGGGCEF9CCFGG?FCFG9E,,B,A,@+=>+,,,,,8,,,,,,,:,,,,<,<,A,,<8,,,8,8,,,,,,,3+++,8,,,,3,,,>>++,,33,,,,,:,,,,8>,,,,,37,,,,,,,,,2,72;;2,,,,,,,,522,,/,2=,,,22,,+,+,,,25*0+4,,++2*1*++++"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT").qualityThresholdRead(5,20)
    // println("test read = " + testRead.bases)
    testRead.bases should be ("CAATTTCCCTACTCAGATCTCGACCTCGAGACAGGTTTGGAGCGAGATTGATAAAGTGCCGACGGTCACGCTTGCGATCTCGTATGCCGTCTTCTGCTTGTAAAAAA")
  }

  "A SequencingRead" should "find the reverse primer correctly in a read" in {
    val readSeq =  "TGTCTCGAGGTCGAGAATTCGGATCCTCCCCCTCGGGCCGCAGCGCCGTTCCATCACCTTCAGCCCCTGTGCCAATCCCCCACTTCCAGCCTCCGTGCCTTACCGCCTCCTCCGGCCACAGAGCCATACCTCCACTCCCAGCTGCAGTGCTTATCCACCACCTGCAGCCACAGTGCCTACCGCCCCCCCCAGCCCCCGTGCCATTCCACCACCCCTAACCCCCCTCCTAAACCTCCCCCCCCACTCCCCGAGCCCC"
    val readQual = "CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGFGGGGGGGGGGGFGGGGGGFGGGGGGGGGGEGGFGGGGGGFGGGGGGGFEFGGFFGDFEDFFDEFFFDC5ECF7CC*8;EEC=DG5@FG5@>5E5CCEEFEEBB7925=EEECC2CFEGC(688?2<-66.97<6:(3(((4::42>0(-,345>?,8<:(49261(.2((-((((-(("
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT").qualityThresholdRead(3,10)
    testRead.startsWithPrimer(Utils.reverseComplement("GAATTCTCGACCTCGAGACA")) should be (true)
  }

  "A SequencingRead" should "should avoid single bad bases when in a 3-window mode" in {
    val readSeq =  "123456789123456"
    val readQual = "HH$HH$HHH$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT").qualityThresholdRead(3,10)
    testRead.bases should be ("123456789")
  }

  "A SequencingRead" should "should catch double bad bases when in a 3-window mode" in {
    val readSeq =  "123456789123456"
    val readQual = "HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT").qualityThresholdRead(3,18)
    testRead.bases should be ("12345")
  }

  "A SequencingRead" should "should calculate a basic zero distance correctly" in {
    val readSeq =  "AATTAATTAATTAAT"
    val readQual = "HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")

    testRead.distance(testRead2) should be (0)
  }

  "A SequencingRead" should "should calculate a half-ish distance correctly" in {
    val readSeq =   "AAAAAAAATTTTTTTT"
    val readSeq2 =  "TTTTTTTTTTTTTTTT"
    val readQual =  "HH$HH$H$H$$$$$$H"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual, ForwardReadOrientation, "TTTT")

    testRead.distance(testRead2) should be (.5)
  }

  "A SequencingRead" should "should calculate a double read length distance correctly" in {
    val readSeq =   "TTTTTTTTTTTTTTT"
    val readSeq2 =  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readQual =  "HH$HH$H$H$$$$$$"
    val readQual2 =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    testRead.distance(testRead2) should be (.5)
  }

  "A SequencingRead" should "should find a cigar of one correctly" in {
    val readSeq =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq2 =  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readQual =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val readQual2 =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    //Array[Tuple2[CigarEvent,Int]]
    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (1)
    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (readSeq.length)
  }

  "A SequencingRead" should "should find one indel right" in {
    val readSeq =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq2 =  "TTTTTTTT--------------TTTTTTTT"
    val readQual =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val readQual2 =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (3)
    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (8)

    cigars(1).event should be (scala.main.Insertion)
    cigars(1).length should be (14)

    cigars(2).event should be (scala.main.Match)
    cigars(2).length should be (8)
  }

  "A SequencingRead" should "should find a bracketed indel with mismatches right" in {
    val readSeq =   "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq2 =  "TTTTTTTTA------------ATTTTTTTT"
    val readQual =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val readQual2 =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (3)
    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (8)

    cigars(1).event should be (scala.main.Insertion)
    cigars(1).length should be (14)

    cigars(2).event should be (scala.main.Match)
    cigars(2).length should be (8)
  }

  "A SequencingRead" should "should find a bracketed deletion with mismatches right" in {
    val readSeq2  =  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq   =  "TTTTTTTTA------------ATTTTTTTT"
    val readQual =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val readQual2 =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (3)
    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (8)

    cigars(1).event should be (scala.main.Deletion)
    cigars(1).length should be (14)

    cigars(2).event should be (scala.main.Match)
    cigars(2).length should be (8)
  }

  "A SequencingRead" should "should sort out a weird mismatch event(s) the right way" in {
    val readSeq2  =  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq   =  "TTTTTTTTATTTTTTTTTTTTATTTTTTTT"
    val readQual =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val readQual2 =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (1)
    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (readSeq.length)
  }

  "A SequencingRead" should "should sort out a weird mismatch event(s)v2 the right way" in {
    val readSeq2  =  "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    val readSeq   =  "ATTTTTTTATTTAATTTTTTTATTTTAATT"
    val readQual =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val readQual2 =  "HH$HH$H$H$$$$$$HH$HH$H$H$$$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (1)
    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (readSeq.length)
  }

  "A SequencingRead" should "should sort out a starting deletion the right way" in {
    val readSeq   =  "----TTTTA"
    val readSeq2  =  "TTTTTTTTT"
    val readQual =   "HH$HH$H$H"
    val readQual2 =  "HH$HH$H$H"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (1)
    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (5)
  }

  "A SequencingRead" should "should sort out a ending deletion the right way" in {
    val readSeq   =  "TTTTATTTA----"
    val readSeq2  =  "TTTTTTTTTTTTT"
    val readQual =   "HH$HH$H$H$$$$"
    val readQual2 =  "HH$HH$H$H$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val testRead2 = SequencingRead(readName, readSeq2, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(testRead2)
    cigars.size should be (1)

    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (9)
  }

  "A SequencingRead" should "should find our position in the read correctly" in {
    val readSeq   =  "---------TTTTATTTA----"
    val refSeq  =    "TTTTTTTTTTTTTTTTTTTTTT"
    val readQual =   "HH$HH$H$H$$$$HH$HH$$$$"
    val readQual2 =  "HH$HH$H$H$$$$H$H$H$$$$"
    val testRead = SequencingRead(readName, readSeq, readQual, ForwardReadOrientation, "TTTT")
    val ref = SequencingRead(readName, refSeq, readQual2, ForwardReadOrientation, "TTTT")

    val cigars = testRead.callCigarsAgainstOtherRead(ref)
    cigars.size should be (1)

    cigars(0).event should be (scala.main.Match)
    cigars(0).length should be (9)
    cigars(0).position should be (9)
  }
}
