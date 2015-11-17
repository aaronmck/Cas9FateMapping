package test.scala


import main.scala.Event
import org.scalatest.FlatSpec
import org.scalatest.FlatSpec
import org.scalatest.Matchers

import org.scalatest.{Matchers, FlatSpec}

/**
 * Created by aaronmck on 11/15/15.
 */
class DistanceMetricTest extends FlatSpec with Matchers {
  "A log2 distance matrix" should "calculate distances correctly" in {
    val evt1 = Event("test1","test1",1,Array[String]("4D-37","79D-65","79D-65","79D-65","79D-65","NONE","NONE","NONE","NONE","NONE"))
    val evt2 = Event("test2","test2",1,Array[String]("41D-20","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))


  }
}
