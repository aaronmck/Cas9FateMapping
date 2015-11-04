package main.scala

import scala.collection.mutable._

// ------------------------------------------------------------------------------------
// our event helper class - stores a single line of our UMI stats table
// ------------------------------------------------------------------------------------
case class Event(name: String, sample: String, numberOfReads: Int, eventStrings: Array[String], knownParent: String = "UNKNOWN") {

  // make a map that contains the actual events and their positions as a way to remove duplicate - multisite events
  val eventMap = new HashMap[String, Array[Int]]()
  eventStrings.zipWithIndex.foreach { case (evt, index) => eventMap(evt) = eventMap.getOrElse(evt, Array[Int]()) :+ index }
  val containsNone = if (eventMap contains "NONE") 1 else 0

  // store our size for quick reference
  val uniqueSize = eventMap.size - containsNone

  // if and when we find a parent record it
  var probParent: Option[Event] = None

  val size = eventStrings.size

  // generate some output
  def toOutputString(): String = {
    val oth = if (probParent.isDefined) probParent.get.eventStrings.mkString("-") + "\t" + probParent.get.sample + "\t" + knownParent else "NONE\tNONE\tNONE"
    name + "\t" + sample + "\t" + numberOfReads + "\t" + eventStrings.mkString(",") + "\t" + oth
  }

  // given an array of events that have an edit distance one less than ours, find a likely parent
  def oneEdit(events: Array[Event]): Boolean = {
    events.foreach { evt => {
      val dst = Event.editDistance(this, evt)
      //println(eventStrings.mkString(",") + "\t" + dst + "\t" + evt.eventStrings.mkString(","))
      if (dst <= 1) {
        probParent = Some(evt)
        return true
      }
    }
    }
    return false
  }
}

//val evt1 = Event("test1","test1",1,Array[String]("4D-37","79D-65","79D-65","79D-65","79D-65","NONE","NONE","NONE","NONE","NONE"))
//val evt2 = Event("test2","test2",1,Array[String]("41D-20","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE","NONE"))

object Event {
  val NONE = "NONE"

  def headerString(): String = "name\tsample\tnumReads\teventString\tparentEventString\tparent\tknownParent"


  def editDistance(evt1: Event, evt2: Event): Int = {
    evt1.eventMap.map { case (event1, positions1) => {
      //println(event1 + "\t" + (evt2.eventMap contains event1))
      if (evt2.eventMap contains event1) {
        val positions2 = evt2.eventMap(event1)
        positions1.map { ps => if (positions2.indexOf(ps) < 0) 1 else 0 }.sum
      } else {
        1
      }
    }
    }.sum
  }

  // calculate the basic distance between two points
  def probJointEvents(event1: String, event2: String,
                      occuranceProportion: HashMap[String, Double],
                      epsilon: Double = 0.001): Double = {

    // both not edited or both the same edit
    if ((event1 == "NONE" && event2 == "NONE") || (event1 != "NONE" && event2 == event1))
      1 - epsilon
    // event 1 is NONE, event2 is not
    else if (event1 == "NONE" && event2 != event1)
      1.0 - occuranceProportion(event2)
    else if (event2 == "NONE" && event2 != event1)
      1.0 - occuranceProportion(event1)
    else // the events should have a distance inversely proportional to occurrence
      epsilon
  }
}
