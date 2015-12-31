package simulator

import main.scala.{InputTable, EventCounter}

import scala.collection.mutable
import scala.util.Random

/**
 * Created by aaronmck on 12/14/15.
 */
class DrawProportionalEvent(events: InputTable) {
  val rand = new Random()


  // setup a container for events over sites
  val positionEventCounts = new mutable.HashMap[Int,mutable.HashMap[String,Int]]()
  val positionEventSize   = new mutable.HashMap[String,Int]()
  val positionTotalCounts = new mutable.HashMap[Int,Int]()
  var totalSum            = 0

  // load up the events and determine a map for site->count and a event->event_size
  events.getAllEvents().foreach{case(eventObj) => {
    val eventLengths = new mutable.HashMap[String,Int]()
    eventObj.getEventStrings().map{case(evt) => eventLengths(evt) = eventLengths.getOrElse(evt,0) + 1}
    eventLengths.foreach{case(evt,cnt) => positionEventSize(evt) = cnt}

    eventObj.getEventStrings().zipWithIndex.foreach {case(event,index) => {
      if (!(positionEventCounts contains index))
        positionEventCounts(index) = new mutable.HashMap[String,Int]

      positionEventCounts(index)(event) = positionEventCounts(index).getOrElse(event,0) + 1
      if (event != "NONE" && event != "UNKNOWN")
        positionTotalCounts(index) = positionTotalCounts.getOrElse(index,0) + 1
    }}
  }}

  totalSum = positionTotalCounts.values.sum

  /**
   * choose a new mutation
   * @param height the height of the tree
   * @return the event string, the position, the the length
   */
  def drawEvent(height: Int, availableSites: Array[Int], subtractProb: Boolean): Tuple3[String,Int,Int] = {

    // first choose a position for the edit by summing up the counts over the available sites
    // and choose randomly, but proportionally, to the total sums
    // ------------------------------------------------------------------------
    var siteTotals = Array[Int]()
    var tCount = 0
    availableSites.foreach{case(availSite) => {
      siteTotals :+= positionTotalCounts(availSite) + tCount
      tCount  += positionTotalCounts(availSite)
    }}

    val totalCount = siteTotals.sum

    val rPosition = rand.nextInt(totalCount)
    var position = -1
    siteTotals.zipWithIndex.foreach{case(total,index) =>
      if (total > rPosition)
        position = availableSites(index)
    }

    if (position < 0)
      throw new IllegalStateException("Unable to pick column")

    //  now that we have a position, choose and event from that position
    // ------------------------------------------------------------------------
    val siteEvents = positionEventCounts(position)

    var retEvent = "NO"
    val rPos = rand.nextInt(positionTotalCounts(position))
    val runningTotal = 0
    siteEvents.foreach{case(event,count) =>
      if ((count + runningTotal) > rPos)
        retEvent = event
    }

    //  now some math -- remove the proportional occurance from the
    // ------------------------------------------------------------------------
    var leavesGivenHeight = math.pow(2,height)
    val newRate = math.max(positionEventCounts(position)(retEvent) - leavesGivenHeight,1)
    positionEventCounts(position)(retEvent) = newRate.toInt
    positionTotalCounts(position) = positionTotalCounts(position) -
      math.min(positionEventCounts(position)(retEvent) - leavesGivenHeight,positionEventCounts(position)(retEvent) - 1).toInt
    return (retEvent,position,positionEventSize(retEvent))
  }

}
