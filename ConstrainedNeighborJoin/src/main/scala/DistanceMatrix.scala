package main.scala

import scala.collection.mutable._

/**
 * create a matrix that tracks the distances and compatibility between pairs
 * of events we've recorded
 *
 * @param events an array of all of our events
 * @param distanceCalculator a calculator that finds the distance between two events
 */
class DistanceMatrix(events: Array[Event], distanceCalculator: DistanceMetric) {

  // our instance variables
  val idToEvent =     new HashMap[Int,Event]()
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
    event.id = lastID
    lastID += 1
  }}
  println("Assigned " + (lastID - 1) + " events")

  /**
   * first setup a mapping of IDs to events.  This mapping MUST stay unique, which we check
   */
  println("Recording events IDs and checking for collisions...")
  events.foreach{event => {
    if (idToEvent contains event.id)
      throw new IllegalStateException("Unable to add duplicate ID " + event.id + " with fancy string " + event.toFancyString())
    idToEvent(event.id) = event
  }}

  /**
   * calculate the compatibility between all known events
   *
   * the compatibility test is taking a loooot of time, cache the results, so we only
   * have to update it for new pairs
   */
  println("Calculating compatability...")
  events.foreach{event1 => {
    compatibility(event1.id) = HashMap[Int,Boolean]()
    events.foreach{event2 => {
      compatibility(event1.id)(event2.id) = Event.compatible(event1,event2)
    }}
  }}

  /**
   * find the distance for each pair of events
   */
  println("calculating initial distances...")
  events.foreach { event => {
    // setup the sum, min distance, and distance map
    var sum = 0.0
    distances(event.id) = new HashMap[Int,Double]()
    var thisMinDist = Double.MaxValue

    events.zipWithIndex.foreach { case (event2, index) => {
      if (event.id != event2.id) {
        val newDist = distanceCalculator.distance(event, event2)
        distances(event.id)(event2.id) = newDist
        sum += newDist

        if (newDist < thisMinDist) {
          minDistance(event.id) = event2.id
          thisMinDist = newDist
        }
      }
    }}
    meanDistances(event.id) = sum / (events.length.toDouble - 2.0)
  }}

  /**
   * find the minimum value for each pair of events
   */
  println("calculating initial optimization functions...")
  events.foreach { event1 => {
    minDistance(event1.id) = -1
    var thisMinDist = Double.MaxValue

    events.zipWithIndex.foreach { case (event2, index) => {
      if (event1.id != event2.id) {
        val minimizationDistance = distances(event1.id)(event2.id) - meanDistances(event1.id) - meanDistances(event2.id)

        if (minimizationDistance < thisMinDist) {
          minDistance(event1.id) = event2.id
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
      meanDistances(eventId) = sum / (idToEvent.size.toDouble - 2.0)
    }}
  }

  /**
   * add a new event to our internal registers -- this involves a fair amount of calculation for mean distances, etc.  It's a
   * bit of a imperative mess, but oh well
   *
   * @param evt the event to add
   */
  def addNewEvent(evt: Event): Unit = {
    // assign this new event
    if (evt.id >= 0)
      throw new IllegalStateException("Unable to incorporate new node with existing id " + evt.id + " fancy.string " + evt.toFancyString())

    evt.id = lastID
    lastID += 1
    //println("adding event " + evt.id)

    /**
     * add this event to the compatibility hashmaps
     */
    compatibility(evt.id) = new HashMap[Int,Boolean]()
    idToEvent.foreach{case(otherID, event2) => {
      val compat = Event.compatible(evt,event2)
      compatibility(evt.id)(otherID) = compat
      compatibility(otherID)(evt.id) = compat
    }}

    /**
     * calculate the distances between the new event and each previous event,
     * updating the mean distances as we go
     */
    // calculate the distances
    distances(evt.id) = new HashMap[Int,Double]()
    var sum = 0.0

    idToEvent.foreach{case(otherID, event2) => {
      val dist = distanceCalculator.distance(evt,event2)

      distances(evt.id)(otherID) = dist
      distances(otherID)(evt.id) = dist
      sum += dist
    }}
    meanDistances(evt.id) = sum / (idToEvent.size.toDouble - 1.0)

    /**
     * now that we have the mean for this node, find it's (and other) min values
     */
    var thisMinDist = Double.MaxValue
    idToEvent.foreach{case(otherID, event2) => {
      val minimizationDistance = distanceCalculator.distance(evt,event2) - meanDistances(evt.id) - meanDistances(event2.id)

      val compat = compatibility(evt.id)(otherID) && compatibility(otherID)(evt.id)

      if (minimizationDistance < thisMinDist && compat) {
        minDistance(evt.id) = otherID
        thisMinDist = minimizationDistance
      }

      val otherMin = minDistance(otherID)
     if (compat && (otherMin < 0 || minimizationDistance < distances(otherID)(otherMin))) {
        minDistance(otherID) = evt.id
      }
    }}
    if (!(minDistance contains evt.id))
      minDistance(evt.id) = -1

    // we add it to the mapping of events last, so that we can loop over this hashmap earlier
    idToEvent(evt.id) = evt
    recalculateMeans()
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

    recalculateMeans()
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
  def findAndMergeNextClosestPair(): Option[Event] = {

    // first find the min distance
    var min = Double.MaxValue
    var minIndex = -1
    var minOtherIndex = -1

    minDistance.foreach { case(eventID, otherMinEvent) => {
      if (otherMinEvent >= 0) {
        val compat = compatibility(eventID)(otherMinEvent) && compatibility(otherMinEvent)(eventID)
        val dst = distances(eventID)(otherMinEvent) - meanDistances(eventID) - meanDistances(otherMinEvent)

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

    // now get the two events, drop them from our internal tracking variables, and return them
    val eventID1 = minIndex
    val eventID2 = minOtherIndex

    // print("Merging " + minIndex + " and " + minOtherIndex + " with distance " + min + " their true distance " + distances(minIndex)(minOtherIndex) + " ")
    // print("sample id " + idToEvent(minIndex).sample + " " + idToEvent(minOtherIndex).sample + " event1 ")
    //println(eventID1 + " event2 " + eventID2 + " in lookup? " + (idToEvent contains eventID1)

    // find the branch lengths
    val branchDisti = (0.5 * distances(eventID1)(eventID2)) + (0.5 * (meanDistances(eventID1) - meanDistances(eventID2)))
    val branchDistj = (0.5 * distances(eventID1)(eventID2)) + (0.5 * (meanDistances(eventID2) - meanDistances(eventID1)))

    // now merge the two events down to a single event
    try {
      val newEvent = Event.merge(idToEvent(eventID1), idToEvent(eventID2), branchDisti, branchDistj, lastID)
      lastID += 1

      // remove the old events
      removeEventID(eventID1)
      removeEventID(eventID2)

      // add the new events
      addNewEvent(newEvent)

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
  def minimizeSet(printStatus: Boolean = true): Array[Event] = {
    if (printStatus)
      print("<Minimize the final set of nodes>")
    while(findAndMergeNextClosestPair().isDefined) {
      if (printStatus)
        print(".")
    }
    val ret = idToEvent.map{tt => tt._2}.toArray

    if (printStatus)
      println("\nMinimize set of size " + ret.size)

    ret
  }
}