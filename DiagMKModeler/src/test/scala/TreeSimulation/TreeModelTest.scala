package test.scala.TreeSimulation

import main.scala.Event
import main.scala.TreeSimulation.TreeModel

import org.scalatest.{Matchers, FlatSpec}

/**
 * Created by aaronmck on 10/31/15.
 */
class TreeModelTest extends FlatSpec with Matchers {

  "An TreeModel" should "split the model into two trees successfully" in {
    val evt1 = Event("test1","test1",1,Array[String]("4D-37","79D-65","79D-65","79D-65","79D-65","NONE","NONE","NONE","NONE","NONE"))
    val evt2 = Event("test2","test2",1,Array[String]("41D-20","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))

    val eventArray = Array[Event](evt1,evt2)
    val treeMod = TreeModel(eventArray,10)

    for (i <- 0 until 100) {
      treeMod.splitOnEvent()
    }

  }
}
