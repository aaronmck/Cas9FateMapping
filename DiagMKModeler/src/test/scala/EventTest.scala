package test.scala

import main.scala.Event
import org.scalatest.{Matchers, FlatSpec}

/**
 * Test out various conditions for edit distances between two sites -- just to make sure we're doing it right
 */
class EventTest extends FlatSpec with Matchers {
  val readName = "TestRead1"

  "An Event" should "find the edit distance between spanning event sites correctly" in {
    val evt1 = Event("test1","test1",1,Array[String]("4D-37","79D-65","79D-65","79D-65","79D-65","NONE","NONE","NONE","NONE","NONE"))
    val evt2 = Event("test2","test2",1,Array[String]("41D-20","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))
    Event.editDistance(evt1,evt2) should be(2)
  }

  "An Event" should "find a simple distance correctly" in {
    // these edit strings don't make biological sense, just here for testing
    val evt1 = Event("test1","test1",1,Array[String]("41D-20","4D-37","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))
    val evt2 = Event("test2","test2",1,Array[String]("41D-20","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))
    Event.editDistance(evt1,evt2) should be(1)
  }
}