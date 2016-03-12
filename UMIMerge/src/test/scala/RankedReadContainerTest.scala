package test.scala

import main.scala.utils.{RankedReadContainer, SortedReads, Utils}
import org.scalatest.junit.JUnitSuite

import scala.main.{ForwardReadOrientation, SequencingRead}
import org.scalatest._

/**
 * Created by aaronmck on 2/11/16.
 */
class RankedReadContainerTest extends FlatSpec with Matchers {
  val readName = "TestRead1"

  "A RankedReadContainer" should "sort two reads by their average quality score correctly" in {
    val read1 = SequencingRead.readFromNameAndSeq("test1","AAAAA","H")
    val read1R = SequencingRead.readFromNameAndSeq("test1","AAAAA","H")

    val read2 = SequencingRead.readFromNameAndSeq("test1","AAAAA","A")
    val read2R = SequencingRead.readFromNameAndSeq("test1","AAAAA","A")
    val rankedContainer1 = SortedReads(read1,read1R)
    val rankedContainer2 = SortedReads(read2,read2R)

    ((rankedContainer1 compare rankedContainer2) < 0) should be (true)
  }

  "A RankedReadContainer" should "handle ranking reads correctly" in {
    val read1 = SequencingRead.readFromNameAndSeq("test1","TTTTT","H")
    val read1R = SequencingRead.readFromNameAndSeq("test1","TTTTT","H")

    val read2 = SequencingRead.readFromNameAndSeq("test1","AAAAA","A")
    val read2R = SequencingRead.readFromNameAndSeq("test1","AAAAA","A")

    val container = new RankedReadContainer("test",2)

    val rankedContainer1 = SortedReads(read1,read1R)
    val rankedContainer2 = SortedReads(read2,read2R)

    container.addRead(rankedContainer1,true,true)
    container.addRead(rankedContainer2,true,true)
    container.pQ.size should be (2)

    container.addRead(rankedContainer1,true,true)
    container.addRead(rankedContainer1,true,true)
    container.pQ.size should be (2)

    container.pQ.dequeue().read1.bases should be ("TTTTT")
    container.pQ.dequeue().read1.bases should be ("TTTTT")
    container.pQ.size should be (0)
  }
}
