import java.io.File

import aligner.{Aligner, MAFFT, AlignmentManager, AlignmentManager$, Alignment}
import org.scalatest.{Matchers, FlatSpec}
import utils.CutSites

import scala.main.{ForwardReadOrientation, SequencingRead}

/**
 * Created by aaronmck on 11/18/15.
 */
class AlignmentManagerTest extends FlatSpec with Matchers {
  val readName = "TestRead1"


  "A MAFFT" should "find basic deletion correctly" in {
    val ref =     "AAATAAAAT"
    val readFwd = "AAAA-AAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = AlignmentManager.callEdits(ref,readFwd,1)
    testCalls.length should be (3)
    testCalls(0).cigarCharacter should be ("M")
    testCalls(0).readBase should be ("AAAA")
    testCalls(0).refBase should be ("AAAT")
    testCalls(0).refBase should be ("AAAT")

    testCalls(1).cigarCharacter should be ("D")
    testCalls(1).readBase should be ("-")
    testCalls(1).refBase should be ("A")

    testCalls(2).cigarCharacter should be ("M")
    testCalls(2).readBase should be ("AAAA")
    testCalls(2).refBase should be ("AAAT")
  }

  "A MAFFT" should "find multibase deletion correctly" in {
    val ref =     "AAATAAAAA"
    val readFwd = "AAAA---AA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = AlignmentManager.callEdits(ref,readFwd,1)
    testCalls.length should be (3)
    testCalls(0).cigarCharacter should be ("M")
    testCalls(0).readBase should be ("AAAA")
    testCalls(0).refBase should be ("AAAT")

    testCalls(1).cigarCharacter should be ("D")
    testCalls(1).readBase should be ("---")
    testCalls(1).refBase should be ("AAA")

    testCalls(2).cigarCharacter should be ("M")
    testCalls(2).readBase should be ("AA")
    testCalls(2).refBase should be ("AA")
  }

  "A MAFFT" should "find multi-deletions correctly" in {
    val ref =     "AAATAAAAT"
    val readFwd = "--AA---AA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = AlignmentManager.callEdits(ref,readFwd,1)
    testCalls.length should be (3)

    testCalls(0).cigarCharacter should be ("M")
    testCalls(0).readBase should be ("AA")
    testCalls(0).refBase should be ("AT")
    testCalls(0).refPos should be (2)

    testCalls(1).cigarCharacter should be ("D")
    testCalls(1).readBase should be ("---")
    testCalls(1).refBase should be ("AAA")
    testCalls(1).refPos should be (4)

    testCalls(2).cigarCharacter should be ("M")
    testCalls(2).readBase should be ("AA")
    testCalls(2).refBase should be ("AT")
    testCalls(2).refPos should be (7)
  }

  "A MAFFT" should "find basic insertion correctly" in {
    val ref =     "AAAT-AAAA"
    val readFwd = "AAAATAAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = AlignmentManager.callEdits(ref,readFwd,1)
    testCalls.length should be (3)
    testCalls(0).cigarCharacter should be ("M")
    testCalls(0).readBase should be ("AAAA")
    testCalls(0).refBase should be ("AAAT")
    testCalls(0).refPos should be (0)

    testCalls(1).cigarCharacter should be ("I")
    testCalls(1).readBase should be ("T")
    testCalls(1).refBase should be ("-")
    testCalls(1).refPos should be (4)

    testCalls(2).cigarCharacter should be ("M")
    testCalls(2).readBase should be ("AAAA")
    testCalls(2).refBase should be ("AAAA")
    testCalls(2).refPos should be (4)
  }

  "A MAFFT" should "find offsets correctly" in {
    val ref =     "AAATAAAAA"
    val readFwd = "----TAAAA"

    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = AlignmentManager.callEdits(ref,readFwd,1)
    testCalls.length should be (1)
    testCalls(0).cigarCharacter should be ("M")
    testCalls(0).readBase should be ("TAAAA")
    testCalls(0).refBase should be  ("AAAAA")
    testCalls(0).refPos should be (4)

  }

  "A MAFFT" should "merge an event and a non-event correctly" in {
    val ref1 =     "AAATAAAAA"
    val readFwd1 = "AAAATAAAA"

    val ref2 =     "AAATAAAAA"
    val readFwd2 = "AAAAT-AAA"

    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls1 = AlignmentManager.callEdits(ref1,readFwd1,1)
    val testCalls2 = AlignmentManager.callEdits(ref2,readFwd2,1)

    testCalls1.length should be (1)
    testCalls2.length should be (3)

    val combined = AlignmentManager.editsToCutSiteCalls(List[List[Alignment]](testCalls1,testCalls2),cutSites)
    combined._2.size should be (1)
    combined._2(0) should be ("1D+5")
  }

  "A MAFFT" should "merge an dual-event and a non-event correctly" in {
    val ref1 =     "AAATAAAAAAAATAAAAA"
    val readFwd1 = "AAAATAAAAAAAATAAAA"

    //              012345678901234567
    val ref2 =     "AAATAAAAAAAATAAAAA"
    val readFwd2 = "AAAAT-AAAAAAAT--AA"

    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7),(12,14,16)))
    val testCalls1 = AlignmentManager.callEdits(ref1,readFwd1,1)
    val testCalls2 = AlignmentManager.callEdits(ref2,readFwd2,1)

    testCalls1.length should be (1)
    testCalls2.length should be (5)

    val combined = AlignmentManager.editsToCutSiteCalls(List[List[Alignment]](testCalls1,testCalls2),cutSites)
    combined._2.size should be (2)
    println(combined._2.mkString(")("))
    combined._2(0) should be ("1D+5")
    combined._2(1) should be ("2D+14")
  }

  "A MAFFT" should "merge an a collision correctly" in {
    val ref1 =     "AAATAAAAAAAATAAAAA"
    val readFwd1 = "AAAAT-AAAAAAAT--AA"

    //              012345678901234567
    val ref2 =     "AAATAAAAAAAATAAAAA"
    val readFwd2 = "AAAAT-AAAAAAAT--AA"

    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7),(12,14,16)))
    val testCalls1 = AlignmentManager.callEdits(ref1,readFwd1,1)
    val testCalls2 = AlignmentManager.callEdits(ref2,readFwd2,1)

    testCalls1.length should be (5)
    testCalls2.length should be (5)

    val combined = AlignmentManager.editsToCutSiteCalls(List[List[Alignment]](testCalls1,testCalls2),cutSites)
    combined._2.size should be (2)
    println(combined._2.mkString(")("))
    combined._2(0) should be ("1D+5")
    combined._2(1) should be ("2D+14")
  }

  "A MAFFT" should "work with real data4" in {
    val ref =     "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNCTTCCTCCAGCTCTTCAGCTCGTCTCTCCAGCAGTTCCCCCGAGTCTGCACCTCCCCAGAAGTCCTCCAGTCCAAACGCTGCTGTCCAGTCTGGCCCGGCGACGGCTCTGTGTGCGGCGTCCAGTCAGGTCGAGGGTTCTGTCAGGACGTCCTGGTGTCCGACCTTCCCAACGGGCCGCAGTATCCTCACTCAGGAGTGGACGATCGAGAGCGATGGCCTTTAGTGTTTTACAACCAAACCTGCCAGTGCGCCGGAAACTACATGGGGTTTGATTGCGGCGAATGCAAGTTCGGCTTCTTCGGTGCCAACTGCGCAGAGAGACGCGAGTCTGTGCGCAGAAATATATTCCAGCTGTCCACTACCGAGAGGCAGAGGTTCATCTCGTACCTAAATCTGGCCAAAACCACCATAAGCCCCGATTATATGATCGTAACAGGAACGTACGCGCAGATGAACAACGGCTCCACGCCAATGTTCGCCAACATCAGTGTGTACGATTTATTCGTGTGGATGCATTATTACGTGTCCCGGGACGCTCTGCTCGGTGGGCCTGGGAATGTGTGGGCTGATAGATCGGAAGAGCACACGTCTGAACT"
    val readFwd =  "CTTCCTCCAGCTCTTCAGCTCGTCTCTCCAGCAGTTCCCCCGAGTCTGCACCACCCAAACGCTGCTGTCCAGTCTGGCCCGGCGACGGCTCCGTGTGCGGCGTCCAGTCAGGTCGAGGGTTCTGTCAGGACGTCCTGGTGTCCGACCTTCCCAACGGGCCGCAGTATCCTCACTCAGTGGACGATCGAGAGCGATGGCCTTTAGTGTTAACACC"

    val readRev = "ATCAGCCCACACATTCCCAGGCCCACCGAGCAGAGCGTCCCGGGACACGTAATAATGCATCCACATAAATCGTACACACTGATGTTGGCGAACATTGGCGTGGAGATTGTTC"

    val fRead = Aligner.SequencingReadFromNameBases("fwd",readFwd)
    val rRead = Aligner.SequencingReadFromNameBases("rev",readRev)

    val cutsSiteObj = CutSites.fromFile(new File("test_files/TYRFull.fasta.cutSites"), 3)

    rRead.reverseCompAlign = true
    println(AlignmentManager.cutSiteEvents("testUMI", ref, fRead, rRead, cutsSiteObj, 10, true)._3.mkString("<->"))
  }

}
