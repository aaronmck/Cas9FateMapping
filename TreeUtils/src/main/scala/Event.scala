package main.scala

/**
  * Created by aaronmck on 4/27/16.
  */
case class Event(eventString: String) {
  val individualEvents = eventString.split(Event.seperator)
  val length = individualEvents.size
}

object Event {
  val seperator = "_"
  val complexEventSeperator = "&"
  val empty = "NONE"

  def intersection(events: Array[Event]): Option[Event] = {
    if (events.size == 0)
      return None

    var intersections = new Array[String](events(0).length)
    throw new IllegalArgumentException("CRAP")

  }
}
