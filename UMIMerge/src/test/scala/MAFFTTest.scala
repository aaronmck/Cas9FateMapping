import java.io.File

import aligner.MAFFT2
import org.scalatest.{Matchers, FlatSpec}
import utils.CutSites

import scala.main.{ForwardReadOrientation, SequencingRead}

/**
 * Created by aaronmck on 11/18/15.
 */
class MAFFTTest extends FlatSpec with Matchers {
  val readName = "TestRead1"


  "A MAFFT" should "find basic deletion correctly" in {
    val ref =     "AAATAAAAT"
    val readFwd = "AAAA-AAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = MAFFT2.callEdits(ref,readFwd,1)
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
    val testCalls = MAFFT2.callEdits(ref,readFwd,1)
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
    val testCalls = MAFFT2.callEdits(ref,readFwd,1)
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
    val testCalls = MAFFT2.callEdits(ref,readFwd,1)
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
    val testCalls = MAFFT2.callEdits(ref,readFwd,1)
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
    val testCalls1 = MAFFT2.callEdits(ref1,readFwd1,1)
    val testCalls2 = MAFFT2.callEdits(ref2,readFwd2,1)

    testCalls1.length should be (1)
    testCalls2.length should be (3)

    val combined = MAFFT2.combineTwoReadEdits(testCalls1,testCalls2,cutSites)
    combined.size should be (1)
    combined(0) should be ("1D+5")
  }

  "A MAFFT" should "merge an dual-event and a non-event correctly" in {
    val ref1 =     "AAATAAAAAAAATAAAAA"
    val readFwd1 = "AAAATAAAAAAAATAAAA"

    //              012345678901234567
    val ref2 =     "AAATAAAAAAAATAAAAA"
    val readFwd2 = "AAAAT-AAAAAAAT--AA"

    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7),(12,14,16)))
    val testCalls1 = MAFFT2.callEdits(ref1,readFwd1,1)
    val testCalls2 = MAFFT2.callEdits(ref2,readFwd2,1)

    testCalls1.length should be (1)
    testCalls2.length should be (5)

    val combined = MAFFT2.combineTwoReadEdits(testCalls1,testCalls2,cutSites)
    combined.size should be (2)
    println(combined.mkString(")("))
    combined(0) should be ("1D+5")
    combined(1) should be ("2D+14")
  }

  "A MAFFT" should "merge an a collision correctly" in {
    val ref1 =     "AAATAAAAAAAATAAAAA"
    val readFwd1 = "AAAAT-AAAAAAAT--AA"

    //              012345678901234567
    val ref2 =     "AAATAAAAAAAATAAAAA"
    val readFwd2 = "AAAAT-AAAAAAAT--AA"

    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7),(12,14,16)))
    val testCalls1 = MAFFT2.callEdits(ref1,readFwd1,1)
    val testCalls2 = MAFFT2.callEdits(ref2,readFwd2,1)

    testCalls1.length should be (5)
    testCalls2.length should be (5)

    val combined = MAFFT2.combineTwoReadEdits(testCalls1,testCalls2,cutSites)
    combined.size should be (2)
    println(combined.mkString(")("))
    combined(0) should be ("1D+5")
    combined(1) should be ("2D+14")
  }

  // CTTCCTCCAGCTCTTCAGCTCGTCTCTCCAGCAGTTCCCCCGAGTCTGCACCTCCCGCCAACATCACCAACATCAGTGTGTACGATTTATTCGTGTGGATGCATTATTACGTGTCCCGGGACGCTCTCCTCGGTGGGCCTGGGAATGTGTGGGCTGATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGGTTGGAAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGCACTCACAAGCTGCACAAACAGTAGCGCGCAAATGAGAGAAGG
  // ATCAGCCCACACATTCCCAGGCCCACCGAGGAGAGCGTCCCGGGACACGTAATAATGCATCCACACGAATAAATCGTACACACTGATGTTGGTGATGTTGGCGGGAGGTGCAGACTCGGGGGAACTGCTGGAGAGACGAGCTGAAGAGCTGGAGGAAGACCCCCAGCACTGTCTCTTATACACATCTGACGCTGCCGACGAAGACGGTCA

  /*"A MAFFT" should "work with real data" in {
    val ref =     "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNCTTCCTCCAGCTCTTCAGCTCGTCTCTCCAGCAGTTCCCCCGAGTCTGCACCTCCCCAGAAGTCCTCCAGTCCAAACGCTGCTGTCCAGTCTGGCCCGGCGACGGCTCTGTGTGCGGCGTCCAGTCAGGTCGAGGGTTCTGTCAGGACGTCCTGGTGTCCGACCTTCCCAACGGGCCGCAGTATCCTCACTCAGGAGTGGACGATCGAGAGCGATGGCCTTTAGTGTTTTACAACCAAACCTGCCAGTGCGCCGGAAACTACATGGGGTTTGATTGCGGCGAATGCAAGTTCGGCTTCTTCGGTGCCAACTGCGCAGAGAGACGCGAGTCTGTGCGCAGAAATATATTCCAGCTGTCCACTACCGAGAGGCAGAGGTTCATCTCGTACCTAAATCTGGCCAAAACCACCATAAGCCCCGATTATATGATCGTAACAGGAACGTACGCGCAGATGAACAACGGCTCCACGCCAATGTTCGCCAACATCAGTGTGTACGATTTATTCGTGTGGATGCATTATTACGTGTCCCGGGACGCTCTGCTCGGTGGGCCTGGGAATGTGTGGGCTGATAGATCGGAAGAGCACACGTCTGAACT"
    val readFwd =  "CTTCCTCCAGCTCTTCAGCTCGTCTCTCCAGCAGTTCCCCCGAGTCTGCACCTCCCGCCAACATCACCAACATCAGTGTGTACGATTTATTCGTGTGGATGCATTATTACGTGTCCCGGGACGCTCTCCTCGGTGGGCCTGGGAATGTGTGGGCTGATAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGGTTGGAAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAGCACTCACAAGCTGCACAAACAGTAGCGCGCAAATGAGAGAAGG"

    val readRev = "ATCAGCCCACACATTCCCAGGCCCACCGAGGAGAGCGTCCCGGGACACGTAATAATGCATCCACACGAATAAATCGTACACACTGATGTTGGTGATGTTGGCGGGAGGTGCAGACTCGGGGGAACTGCTGGAGAGACGAGCTGAAGAGCTGGAGGAAGACCCCCAGCACTGTCTCTTATACACATCTGACGCTGCCGACGAAGACGGTCA"

    val fRead = MAFFT2.SequencingReadFromNameBases("fwd",readFwd)
    val rRead = MAFFT2.SequencingReadFromNameBases("rev",readRev)

    val cutsSiteObj = CutSites.fromFile(new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_18_Full_TYR_sequencing/./data/TYRFull.fasta.cutSites"), 3)

    MAFFT2.alignTo(Array[SequencingRead](fRead),Some(ref), false,true)
    MAFFT2.alignTo(Array[SequencingRead](rRead),Some(ref), true ,true)
    MAFFT2.cutSiteEvents(ref, fRead, rRead, cutsSiteObj, true)
  }

  "A MAFFT" should "work with real data2" in {
    val ref =     "ATCTCGAGCTCAAGCTTCGGGCACTGCGGCTGGAGGTGGGGGAAAGGCACTGAGACTGGGGGTGGAGGTTTAGCAGTGCGGCTAGAGGTGGTGGAATGGCACTGGGGCTGGGGGAGGCGGTTAGGCACTGTGGCTGCAGGTGGTGGATAAGCACTGCAGCTGGGAGTGGAGGTATGGCTCTGTGGCCGGAGGAGGCGGTAAGGCACGGAGGCTGGAAGTGGGGGATTGGCACAGGGGCTGAAGGTGATGGAACGGCGCTGCGGCCCGAGGGGGAGGATCCGAATTCTCGACCTCGAGACAAGATCGGAAGAGCACACGTCTGAACT"
    val readFwd =  "ACTCAGATCTCGAGCTCAAGCTTCGGGCACTGCGGCTGGGGGAAAGGCACTGAGACTGGGGGTGGAGGTTTAGCAGTGCGGCTAGAGGTGGTGGAATGGCACTGGGGCTGGGGGAGGCGGTTAGGCACTGTGGCTGCAGGTGGTGGATAAGCACTGCAGCTGGGAGTGGAGGTATGGCTCTGTGGCCGGAGGAGGCGGTAAGGCACGGAGGCTGGAAGTGGGGGATTGGCACAGGGGCTGA"

    val readRev = "TGTCTCGAGGTCGAGAATTCGGATCCTCCCCCTCGGGCCGCAGCGCCGTTCCATCACCTTCAGCCCCTGTGCCAATCCCCCACTTCCAGCCTCCGTGCCTTACCGCCTCCTCCGGCCACAGAGCCATACCTCC"

    val fRead = MAFFT2.SequencingReadFromNameBases("fwd",readFwd)
    val rRead = MAFFT2.SequencingReadFromNameBases("rev",readRev)

    val cutsSiteObj = CutSites.fromFile(new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_18_Full_TYR_sequencing/./data/TYRFull.fasta.cutSites"), 3)

    MAFFT2.alignTo(Array[SequencingRead](fRead),Some(ref), false,true)
    MAFFT2.alignTo(Array[SequencingRead](rRead),Some(ref), true ,true)
    MAFFT2.cutSiteEvents(ref, fRead, rRead, cutsSiteObj, 10, true)
  }

  "A MAFFT" should "work with real data3" in {
    val ref =     "ATCTCGAGCTCAAGCTTCGGGCACTGCGGCTGGAGGTGGGGGAAAGGCACTGAGACTGGGGGTGGAGGTTTAGCAGTGCGGCTAGAGGTGGTGGAATGGCACTGGGGCTGGGGGAGGCGGTTAGGCACTGTGGCTGCAGGTGGTGGATAAGCACTGCAGCTGGGAGTGGAGGTATGGCTCTGTGGCCGGAGGAGGCGGTAAGGCACGGAGGCTGGAAGTGGGGGATTGGCACAGGGGCTGAAGGTGATGGAACGGCGCTGCGGCCCGAGGGGGAGGATCCGAATTCTCGACCTCGAGACAAGATCGGAAGAGCACACGTCTGAACT"
    val readFwd =  "ACTCAGATCTCGAGCTCAAGCTTCGGGCACTGCGGCTGGGGGAAAGGCACTGAGACTGGGGGTGGAGGTTTAGCAGTGCGGCTAGAGGTGGTGGAATGGCACTGGGGCTGGGGGAGGCGGTTAGGCACTGTGGCTGCAGGTGGTGGATAAGCACTGCAGCTGGGAGTGGAGGTATGGCTCTGTGGCCGGAGGAGGCGGTAAGGCACGGAGGCTGGAAGTGGGGGATTGGCACAGGGGCTGAAGGTGATGG"

    val readRev = "TGTCTCGAGGTCGAGAATTCGGATCCTCCCCCTCGGGCCGCAGCGCCGTTCCATCACCTTCAGCCCCTGTGCCAATCCCCCACTTCCAGCCTCCGTGCCTTACCGCCTCCTCCGGCCACAGAGCCATACCTCCACTCCCAGCTGCAGTGCTTATCCCCACC"

    val fRead = MAFFT2.SequencingReadFromNameBases("fwd",readFwd)
    val rRead = MAFFT2.SequencingReadFromNameBases("rev",readRev)

    val cutsSiteObj = CutSites.fromFile(new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_18_Full_TYR_sequencing/./data/TYRFull.fasta.cutSites"), 3)

    MAFFT2.alignTo(Array[SequencingRead](fRead),Some(ref), false,true)
    MAFFT2.alignTo(Array[SequencingRead](rRead),Some(ref), true ,true)
    MAFFT2.cutSiteEvents(ref, fRead, rRead, cutsSiteObj, 10, true)
  }*/

  "A MAFFT" should "work with real data4" in {
    val ref =     "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNNCTTCCTCCAGCTCTTCAGCTCGTCTCTCCAGCAGTTCCCCCGAGTCTGCACCTCCCCAGAAGTCCTCCAGTCCAAACGCTGCTGTCCAGTCTGGCCCGGCGACGGCTCTGTGTGCGGCGTCCAGTCAGGTCGAGGGTTCTGTCAGGACGTCCTGGTGTCCGACCTTCCCAACGGGCCGCAGTATCCTCACTCAGGAGTGGACGATCGAGAGCGATGGCCTTTAGTGTTTTACAACCAAACCTGCCAGTGCGCCGGAAACTACATGGGGTTTGATTGCGGCGAATGCAAGTTCGGCTTCTTCGGTGCCAACTGCGCAGAGAGACGCGAGTCTGTGCGCAGAAATATATTCCAGCTGTCCACTACCGAGAGGCAGAGGTTCATCTCGTACCTAAATCTGGCCAAAACCACCATAAGCCCCGATTATATGATCGTAACAGGAACGTACGCGCAGATGAACAACGGCTCCACGCCAATGTTCGCCAACATCAGTGTGTACGATTTATTCGTGTGGATGCATTATTACGTGTCCCGGGACGCTCTGCTCGGTGGGCCTGGGAATGTGTGGGCTGATAGATCGGAAGAGCACACGTCTGAACT"
    val readFwd =  "CTTCCTCCAGCTCTTCAGCTCGTCTCTCCAGCAGTTCCCCCGAGTCTGCACCACCCAAACGCTGCTGTCCAGTCTGGCCCGGCGACGGCTCCGTGTGCGGCGTCCAGTCAGGTCGAGGGTTCTGTCAGGACGTCCTGGTGTCCGACCTTCCCAACGGGCCGCAGTATCCTCACTCAGTGGACGATCGAGAGCGATGGCCTTTAGTGTTAACACC"

    val readRev = "ATCAGCCCACACATTCCCAGGCCCACCGAGCAGAGCGTCCCGGGACACGTAATAATGCATCCACATAAATCGTACACACTGATGTTGGCGAACATTGGCGTGGAGATTGTTC"

    val fRead = MAFFT2.SequencingReadFromNameBases("fwd",readFwd)
    val rRead = MAFFT2.SequencingReadFromNameBases("rev",readRev)

    val cutsSiteObj = CutSites.fromFile(new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_10_22_Cell_Culture_Lineage/data/HEK4_t1_20.fa.cutSites"), 3)

    MAFFT2.alignTo(Array[SequencingRead](fRead),Some(ref), false,true)
    MAFFT2.alignTo(Array[SequencingRead](rRead),Some(ref), true ,true)
    MAFFT2.cutSiteEvents(ref, fRead, rRead, cutsSiteObj, 10, true)
  }

  /**
   * OLD TESTS
   * simple, single-event tests


  "A MAFFT" should "find basic deletion correctly" in {
    val readRef = "AAAAAAAAA"
    val readFwd = "AAAA-AAAA"
    val readRev = "AAAA-AAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (1)
    testCalls._2(0) should be ("1D-4")
  }

  "A MAFFT" should "find basic 2-base deletion correctly" in {
    val readRef = "AAAAAAAAA"
    val readFwd = "AAAA--AAA"
    val readRev = "AAAA--AAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (1)
    testCalls._2(0) should be ("2D-4")
  }

  "A MAFFT" should "find basic 1-base insertion correctly" in {
    val readRef = "AAAA-AAAA"
    val readFwd = "AAAAAAAAA"
    val readRev = "AAAAAAAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (1)
    testCalls._2(0) should be ("1I-4")
  }


  "A MAFFT" should "find basic 2-base insertion correctly" in {
    val readRef = "AAAA--AAA"
    val readFwd = "AAAAAAAAA"
    val readRev = "AAAAAAAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (1)
    testCalls._2(0) should be ("2I-4")
  }


  /**
   * two site tests
   */

  "A MAFFT" should "find a 2-site basic deletion correctly" in {
    // position:   012345678901234567
    val readRef = "AAAAAAAAAAAAAAAAAA"
    val readFwd = "AAAA-AAAAAAAAAAAAA"
    val readRev = "AAAA-AAAAAAAAAAAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7), (12,14,16)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (2)
    testCalls._2(0) should be ("1D-4")
    testCalls._2(1) should be ("")
  }

  "A MAFFT" should "find a paired 2-site basic deletion correctly" in {
    // position:   012345678901234567
    val readRef = "AAAAAAAAAAAAAAAAAA"
    val readFwd = "AAAA-AAAAAAAAA-AAA"
    val readRev = "AAAA-AAAAAAAAA-AAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7), (12,14,16)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (2)
    testCalls._2(0) should be ("1D-4")
    testCalls._2(1) should be ("1D-14")
  }

  "A MAFFT" should "find a 2-site mixed insertion/deletion correctly" in {
    // position:   012345678901234567
    val readRef = "AAAAAAAAAAAAAA-AAA"
    val readFwd = "AAAA-AAAAAAAAAAAAA"
    val readRev = "AAAA-AAAAAAAAAAAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7), (12,14,16)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (2)
    testCalls._2(0) should be ("1D-4")
    testCalls._2(1) should be ("1I-14")
  }

  "A MAFFT" should "find a 2-site paired insertions correctly" in {
    // position:   012345456789123456
    val readRef = "AAAA-AAAAAAAAA-AAA"
    val readFwd = "AAAAAAAAAAAAAAAAAA"
    val readRev = "AAAAAAAAAAAAAAAAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7), (12,14,16)))
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (2)
    testCalls._2(0) should be ("1I-4")
    testCalls._2(1) should be ("1I-13")
  }

  "A MAFFT" should "find a 2-site longer paired insertions correctly" in {
    // position:   012344445678901234
    val readRef = "AAAA---AAAAAAA-AAA"
    val readFwd = "AAAAAAAAAAAAAAAAAA"
    val readRev = "AAAAAAAAAAAAAAAAAA"
    val cutSites = CutSites.fromIntervals(Array[Tuple3[Int,Int,Int]]((3,5,7), (11,14,16))) // lower the range here to capture the other event
    val testCalls = MAFFT2.callEdits(readRef,readFwd,readRev,cutSites)
    testCalls._2.length should be (2)
    testCalls._2(0) should be ("3I-4")
    testCalls._2(1) should be ("1I-11")
  } */

}
