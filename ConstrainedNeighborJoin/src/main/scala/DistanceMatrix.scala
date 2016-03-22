package main.scala

import java.io.{PrintWriter, File}

import scala.collection.mutable
import scala.collection.mutable._

/**
 * create a matrix that tracks the distances and compatibility between pairs
 * of events we've recorded
 *
 * @param events an array of all of our events
 * @param distanceCalculator a calculator that finds the distance between two events
 */
class DistanceMatrix(events: Array[IndexedNode], distanceCalculator: DistanceMetric[IndexedNode]) {

  // should we print debug info?
  val debug = false

  // our instance variables
  val idToEvent =     new HashMap[Int,IndexedNode]()
  val compatibility = new HashMap[Int,HashMap[Int,Boolean]]()
  val distances =     new HashMap[Int,HashMap[Int,Double]]()
  val meanDistances = new HashMap[Int,Double]()
  val minDistance =   new HashMap[Int,Int]()

  // we keep track of every event's ID
  var lastID = 0

  // a retained list of removed IDs
  var removedIDs = new HashMap[Int,Boolean]()

  /**
   * assign an ID to each event in order.  This ID mapping gets assigned by this distance matrix
   */
  println("Assigning event ids...")
  events.foreach{event => {
    event.setID(lastID)
    lastID += 1
  }}
  println("Assigned " + (lastID - 1) + " events")

  /**
   * first setup a mapping of IDs to events.  This mapping MUST stay unique, which we check
   */
  println("Recording events IDs and checking for collisions...")
  events.foreach{event => {
    if (idToEvent contains event.getID())
      throw new IllegalStateException("Unable to add duplicate ID " + event.getID() + " with fancy string ")
    idToEvent(event.getID()) = event
  }}

  /**
   * calculate the compatibility between all known events
   *
   * the compatibility test is taking a loooot of time, cache the results, so we only
   * have to update it for new pairs
   */
  println("Calculating compatability...")
  events.foreach{event1 => {
    compatibility(event1.getID()) = HashMap[Int,Boolean]()
    events.foreach{event2 => {
      compatibility(event1.getID())(event2.getID()) = event1.compatible(event2)
    }}
  }}

  /**
   * find the distance for each pair of events
   */
  println("calculating initial distances...")
  events.foreach { event => {
    // setup the sum, min distance, and distance map
    var sum = 0.0
    distances(event.getID()) = new HashMap[Int,Double]()
    var thisMinDist = Double.MaxValue

    events.zipWithIndex.foreach { case (event2, index) => {
      if (event.getID() != event2.getID()) {
        val newDist = distanceCalculator.distance(event, event2)
        distances(event.getID())(event2.getID()) = newDist
        sum += newDist

        if (newDist < thisMinDist) {
          minDistance(event.getID()) = event2.getID()
          thisMinDist = newDist
        }
      }
    }}
    meanDistances(event.getID()) = sum / (math.max(3,events.length.toDouble) - 2.0)
  }}

  /**
   * find the minimum value for each pair of events
   */
  println("calculating initial optimization functions...")
  events.foreach { event1 => {
    minDistance(event1.getID()) = -1
    var thisMinDist = Double.MaxValue

    events.zipWithIndex.foreach { case (event2, index) => {
      if (event1.getID() != event2.getID()) {
        val minimizationDistance = distances(event1.getID())(event2.getID()) - meanDistances(event1.getID()) - meanDistances(event2.getID())

        if (minimizationDistance < thisMinDist) {
          minDistance(event1.getID()) = event2.getID()
          thisMinDist = minimizationDistance
        }
      }
    }}
  }}

  /**
   * out basic distance metric -- the distance is the D_ij - Mu_i - Mu_j
   * @param id1 the first node's ID
   * @param id2 the second node's ID
   * @return a double representing the distance metric we're optimizing
   */
  def distanceFunction(id1: Int, id2: Int): Double = {
    distances(id1)(id2) - meanDistances(id1) - meanDistances(id2)
  }


  /**
   * recalculate the mean distance from one node to all others -- this must
   * be recomputed after each merging
   */
  def recalculateMeans(): Unit = {
    idToEvent.foreach{case(eventId, event) => {
      var sum = 0.0
      idToEvent.foreach{case(eventId2, event2) => {
        if (eventId != eventId2) {
          sum += distances(eventId)(eventId2)
        }
      }}
      meanDistances(eventId) = sum / (math.max(3,idToEvent.size.toDouble) - 2.0)
    }}
  }

  /**
   * add a new event to our internal registers -- this involves a fair amount of calculation for mean distances, etc.  It's a
   * bit of a imperative mess, but oh well.  Also worth noting:
   *
   * calculating distance here is different!
   *
   * @param evt the event to add
   */
  def addNewInternalEvent(evt: IndexedNode, mergeFrom1id: Int, mergeFrom2id: Int): Unit = {
    if (debug)
      println("addinging node " + lastID + " removing " + mergeFrom1id + " and " + mergeFrom2id + " " + (distances contains mergeFrom1id) + " " + (distances contains mergeFrom1id))

    // assign this new event
    if (evt.getID() >= 0)
      throw new IllegalStateException("Unable to incorporate new node with existing id " + evt.getID() + " fancy.string " + evt.toFancyString())

    evt.setID(lastID)
    lastID += 1
    //println("adding event " + evt.id)

    /**
     * add this event to the compatibility hashmaps
     */
    compatibility(evt.getID()) = new HashMap[Int,Boolean]()
    idToEvent.foreach{case(otherID, event2) => {
      val compat = evt.compatible(event2)
      compatibility(evt.getID())(otherID) = compat
      compatibility(otherID)(evt.getID()) = compat
    }}

    /**
     * calculate the distances between the new event and each previous event,
     * updating the mean distances as we go
     */
    // calculate the distances
    distances(evt.getID()) = new HashMap[Int,Double]()
    var sum = 0.0

    val oneToTwo = distances(mergeFrom1id)(mergeFrom2id)

    idToEvent.foreach{case(otherID, event2) => {

      if (otherID != mergeFrom1id && otherID != mergeFrom2id) {
        val oneToOther = distances(mergeFrom1id)(otherID)
        val twoToOther = distances(mergeFrom2id)(otherID)

        val dist = (oneToOther + twoToOther - oneToTwo)/2.0
        if (debug) {
          print("ADDING dist " + dist + " for one " + idToEvent(otherID).getEventStrings()(0))
          println(" distance between " + oneToTwo + " oneToOther " + oneToOther + " twoToOther " + twoToOther)
        }

        distances(evt.getID())(otherID) = dist
        distances(otherID)(evt.getID()) = dist
        sum += dist
      }
    }}
    meanDistances(evt.getID()) = sum / (math.max(3,idToEvent.size.toDouble) - 1.0)

    /**
     * now that we have the mean for this node, find it's (and other) min values
     */
    var thisMinDist = Double.MaxValue
    idToEvent.foreach{case(otherID, event2) => {
      val minimizationDistance = distanceCalculator.distance(evt,event2) - meanDistances(evt.getID()) - meanDistances(event2.getID())

      val compat = compatibility(evt.getID())(otherID) && compatibility(otherID)(evt.getID())

      if (minimizationDistance < thisMinDist && compat) {
        minDistance(evt.getID()) = otherID
        thisMinDist = minimizationDistance
      }

      val otherMin = minDistance(otherID)
     if (compat && (otherMin < 0 || minimizationDistance < distances(otherID)(otherMin))) {
        minDistance(otherID) = evt.getID()
      }
    }}
    if (!(minDistance contains evt.getID()))
      minDistance(evt.getID()) = -1

    // we add it to the mapping of events last, so that we can loop over this hashmap earlier
    idToEvent(evt.getID()) = evt

  }

  /**
   * remove an event from our all of our internal registers by ID
   *
   * @param idToRemove the ID to remove
   */
  def removeEventID(idToRemove: Int): Unit = {
    //if (idToRemove == 161)
    //println("removing event " + idToRemove)

    removedIDs(idToRemove) = true

    idToEvent.remove(idToRemove)
    compatibility.remove(idToRemove)
    compatibility.foreach{case(id,mapToComp) => {
      if (mapToComp contains idToRemove) {
        mapToComp.remove(idToRemove)
      }
    }}

    distances.remove(idToRemove)
    distances.foreach{case(id,mapToComp) => {
      if (mapToComp contains idToRemove) {
        mapToComp.remove(idToRemove)
      }
    }}

    meanDistances.remove(idToRemove)
    minDistance.remove(idToRemove)

    // now fix the min distances: if this was a min from another node,
    // we need to find that node a new min

    minDistance.keySet.foreach {case(fromNode) => {
      if (minDistance(fromNode) == idToRemove) {
        minDistance(fromNode) = findNewMin(fromNode)
      }
    }}
  }

  /**
   * find a new min distance from a indexed node
   * @param fromNode the node to find the min for
   * @return the new node index
   */
  def findNewMin(fromNode: Int): Int = {
    var newMin = Double.MaxValue
    var newIndex = -1
    distances(fromNode).foreach{case(otherNode,distance) => {
      val dst = distance - meanDistances(fromNode) - meanDistances(otherNode)
      if (fromNode != otherNode && dst < newMin) {
        newMin = dst
        newIndex = otherNode
      }
    }}
    newIndex
  }

  /**
   * find the next two closest events, and return this pair, if there isn't a valid pair return None
   * @return the pair of events that are closest: THIS is destructive, removing that pair
   */
  def findAndMergeNextClosestPair(useConstrained: Boolean): Option[IndexedNode] = {

    // first find the min distance
    var min = Double.MaxValue
    var minIndex = -1
    var minOtherIndex = -1

    minDistance.foreach { case(eventID, otherMinEvent) => {
      if (otherMinEvent >= 0) {
        val compat = if (useConstrained)
          compatibility(eventID)(otherMinEvent) && compatibility(otherMinEvent)(eventID)
        else
            true

        val dst = distances(eventID)(otherMinEvent) - meanDistances(eventID) - meanDistances(otherMinEvent)
        if (debug)
          println("DST: id = " + eventID + " other id = " + otherMinEvent + " dst = " + dst + " distance = " + distances(eventID)(otherMinEvent) + " mean1 " + meanDistances(eventID) + " mean2 "+  meanDistances(otherMinEvent))

        if (eventID != otherMinEvent && dst < min && compat) {
          min = dst
          minIndex = eventID
          minOtherIndex = otherMinEvent
        }
      }
    }}

    // if we couldn't find a matching event because no remaining events were compatible, return None
    if (minIndex < 0)
      return None

    if (debug)
      println("MIN is " + min)
    // now get the two events, drop them from our internal tracking variables, and return them
    val eventID1 = minIndex
    val eventID2 = minOtherIndex

    // print("Merging " + minIndex + " and " + minOtherIndex + " with distance " + min + " their true distance " + distances(minIndex)(minOtherIndex) + " ")
    // print("sample id " + idToEvent(minIndex).sample + " " + idToEvent(minOtherIndex).sample + " event1 ")
    //println(eventID1 + " event2 " + eventID2 + " in lookup? " + (idToEvent contains eventID1)

    // find the branch lengths
    val branchDisti = (0.5 * distances(eventID1)(eventID2)) + (0.5 * (meanDistances(eventID1) - meanDistances(eventID2)))
    val branchDistj = (0.5 * distances(eventID1)(eventID2)) + (0.5 * (meanDistances(eventID2) - meanDistances(eventID1)))

    if (debug)
      println("branchi = " + branchDisti + " branchDistj " + branchDistj)

    // now merge the two events down to a single event
    try {
      val newEvent = idToEvent(eventID1).merge(idToEvent(eventID2), branchDisti, branchDistj, lastID, useConstrained)
      lastID += 1

      // add the new events
      addNewInternalEvent(newEvent,eventID1,eventID2)

      // remove the old events
      removeEventID(eventID1)
      removeEventID(eventID2)

      // find new mean values
      recalculateMeans()

      if (debug)
        idToEvent.foreach{case(index,ode)  => {
          println(ode.getEventStrings()(0) + " => " + idToEvent.map{case(index2,ode) =>
            if (index != index2) ode.getEventStrings()(0) + "-" + distances(index)(index2) else ode.getEventStrings()(0) + "-" + 0.0}.mkString("\t"))
        }}

      return(Some(newEvent))
    } catch {
      case e: Exception => {
        println("eventID1 " + eventID1)
        println("eventID2 " + eventID2)
        println("does the idToEvent contain1 " + (idToEvent contains eventID1))
        println("does the idToEvent contain2 " + (idToEvent contains eventID2))
        //println("Failed merging " + eventID1 + " and " + eventID2 + " with compatibility " + compatibility(eventID1)(eventID2))
        throw(e)
      }
    }

  }

  /**
   * given a full set, merge down compatible, minimum distance events until we've come to the smallest possible set
   *
   * @return the set of minimum values
   */
  def minimizeSet(useConstraints: Boolean, printStatus: Boolean = true): Array[IndexedNode] = {
    while(findAndMergeNextClosestPair(useConstraints).isDefined) {
      if (printStatus)
        print(".")
    }
    val ret = idToEvent.map{tt => tt._2}.toArray

    if (printStatus)
      println("\nMinimize set of size " + ret.size)

    ret
  }

  /**
   * write the distance matrix for our file
   * @param outputFile the file to write to
   */
  def toDistanceFile(outputFile: File): Unit = {
    val output = new PrintWriter(outputFile.getAbsolutePath)

    output.write("othernode\t" + 0.until(lastID).map{case(id) => id}.mkString("\t") + "\n")
    0.until(lastID).map{case(id) => {
      output.write(id + "\t" + 0.until(lastID).map{case(thisid) => distances(id).getOrElse(thisid,0.0)}.mkString("\t") + "\n")
    }}
    output.close()
  }

  /**
   * write the annotation file out to disk
   * @param outputFile the output file to write
   */
  def toAnnotationFile(outputFile: File, stats: InputTable): Unit = {
    val output = new PrintWriter(outputFile.getAbsolutePath)

    output.write("taxa\tsample\tcount\teventString\tproportion\n")
    stats.getUniqueEvents().foreach(ft => {
      val sampleTotal = stats.getSampleCount()(ft.getSample()).toDouble
      output.write(ft.getID() + "\t" + ft.getSample() + "\t" +
        ft.getSupport() + "\t" + ft.getEventStrings().mkString("-") + "\t" + (ft.getSupport().toDouble /sampleTotal) + "\n")
    })

    output.close()
  }

}