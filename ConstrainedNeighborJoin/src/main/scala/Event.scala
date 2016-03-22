package main.scala

import java.io.PrintWriter

import scala.util.Random

// ------------------------------------------------------------------------------------------------------------------------
//
// Helper class -- contains most things we know about one of our events
//
// ------------------------------------------------------------------------------------------------------------------------
case class Event(name: String, sample: String, supportCount: Int, eventStrings: Array[String], inferred: Boolean = false) extends IndexedNode {

  var left: Option[Event] = None
  var right: Option[Event] = None
  var branchLeft = 0.0
  var branchRight = 0.0
  var id = -1
  var support_count = supportCount

  def setID(idVal: Int) = {id = idVal}
  def getID() = id
  def getSupport() = support_count
  def getSample() = sample
  def addSupport(additionalCount: Int) { support_count += additionalCount}
  def countNonWT(): Int = eventStrings.filter{evt => evt.contains("\\+")}.size

  def merge(otherNode: IndexedNode, branchLeft: Double, branchRight: Double, newID: Int, constrained: Boolean): IndexedNode = otherNode match {
    case other: Event => Event.merge(this,otherNode.asInstanceOf[Event],branchLeft,branchRight,newID,constrained)
    case _ => throw new IllegalStateException("Trying to merge a non-Event of IndexedNode with an Event")
  }

  def compatible(otherNode: IndexedNode): Boolean = otherNode match {
    case other: Event => Event.compatible(this, otherNode.asInstanceOf[Event])
    case _ => throw new IllegalStateException("Trying to compare a non-Event of IndexedNode with an Event")
  }

  def getEventStrings() = eventStrings

  /**
   * make a fancy string from this event for debug printing
   * @return a string representation
   */
  def toFancyString(): String = name + "\t" + sample + "\t" + support_count + "\t" +
    eventStrings.mkString("\t") + "\t" + left.isDefined + "\t" + right.isDefined

  /**
   * create the Newick string representation of this node
   * @param distance the distance from this node to it's parent
   * @param annotationOutput where to write the string to (to keep it easy)
   * @return the string representation
   */
  def newickString(distance: Double, fullDistance: Double, annotationOutput: PrintWriter): String = {
    val leftStr =  if (left.isDefined)  left.get.newickString(branchLeft, distance+branchLeft, annotationOutput)   else ""
    val rightStr = if (right.isDefined) right.get.newickString(branchRight, distance+branchRight,annotationOutput) else ""

    var newDist = if (distance.isNaN) 1.0 else distance
    val ret = id + ":" + newDist

    var newSample = if (sample startsWith "M") "_" else sample
    annotationOutput.write(id + "\t" + name + "\t" + newSample + "\t" + fullDistance + "\t" + support_count + "\t" + eventStrings.mkString("-") + "\n")

    //println(ret + "(" + leftStr + "," + rightStr + ")")
    val rnd = new Random()

    if (leftStr != "")
      if (rnd.nextDouble < 0.5)
        "(" + leftStr + "," + rightStr + ")" + ret
      else
        "(" + rightStr + "," + leftStr + ")" + ret
    else
      ret
  }
}


object Event {
  def scoreInternal(subNode: Event, event2: Event): Int = {
    val sm = subNode.eventStrings.zip(event2.eventStrings).map{case(evt1,evt2) =>
      if (evt1 == "NONE" || evt2 == "NONE" || evt2 == evt1)
        0
      else
        1
    }.sum
    /*if (subNode.left.isDefined)
      sm += Event.scoreInternal(subNode.left.get,event2)
    if (subNode.right.isDefined)
      sm += Event.scoreInternal(subNode.right.get,event2)*/
    return sm
  }

  /**
   * merge down this event with the passed in event.  This 'merge' means a union: any
   * events they share are kept, and NONEs they share are kept, and any only events where
   * one is NONE and other has a value become NONE in the merge.
   * @param event1 the first event
   * @param event2 the other event
   * @param branchL the branch to the left node (this node)
   * @param branchR the branch to the right node (the event2 distance)
   * @param mergeID the id to assign the return 'merged' node
   * @return a merged event
   */
  def merge(event1: Event, event2: Event, branchL: Double, branchR: Double, mergeID: Int, constrained: Boolean): Event = {
    if (event1.id == event2.id) {
      throw new IllegalStateException("Don't ask to merge two events with the same ID: " + event1.id)
    }

    val newEvtStr = event1.eventStrings.zip(event2.eventStrings).map{case(evt1,evt2) => {
      if (evt1 != "NONE" && evt2 != "NONE" && evt2 != evt1) {
        if (constrained) {
          println()
          println(Event.compatible(event1, event2))
          println(event1.eventStrings.mkString("\t"))
          println(event2.eventStrings.mkString("\t"))
          throw new IllegalStateException("UNABLE TO MERGE " + event1.toFancyString + " and " + event2.toFancyString)
        } else {
          "NONE"
        }
      }

      if ((evt1 == "NONE" && evt2 != "NONE") || (evt2 == "NONE" && evt1 != "NONE") ) {
        "NONE"
      } else {
        if (evt1 == evt2)
          evt1
        else
          "NONE"
      }
    }}

    val ret = Event("MERGED" + mergeID,
      "M_L" + event1.sample + "_" + event2.sample + "J",
      event1.support_count + event2.support_count,
      newEvtStr,
      true)

    ret.left = Some(event1)
    ret.right = Some(event2)
    ret.branchLeft = branchL
    ret.branchRight = branchR

    return ret
  }

  /**
   * are two events compatible? they're compatible if for all sites they either:
   * 1) share the same edits
   * 2) either or both have a NONE
   * they must satisfy these criteria at every position to be compatible
   *
   * @param event2 another event to compare to
   * @return true if they're compatible, false if not
   */
  def compatible(event1: Event, event2: Event): Boolean = {
    return true /*
    // now check our compatibility
    var sm = 0
    sm += Event.scoreInternal(event1, event2)
    sm += Event.scoreInternal(event2, event1)

    // our compatible condition and do the reverse in case event2 has merged nodes
    sm == 0*/
  }

}