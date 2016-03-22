package main.scala

import java.io.File

import scala.collection.mutable
import scala.collection.mutable.HashMap

/**
 * Created by aaronmck on 11/25/15.
 */
abstract class EventCounter(eventFile: InputTable) {
  def getEventCount(sample: String, event: String): Double
  def getAllSampleEvents(sample: String): HashMap[String,Double]
  def getAllEvents(): HashMap[String,Double]
}

case class NormalizedEventCounter(eventF: InputTable, numberOfCells: Double) extends EventCounter(eventF) {

  // a mapping of sample -> event -> count
  val sampleEventToCount = new HashMap[String, HashMap[String, Double]]()
  val eventToCountAll = new HashMap[String, Double]()

  eventF.getAllEvents().foreach { case(evt1) => {
    if (!(sampleEventToCount contains evt1.getSample()))
      sampleEventToCount(evt1.getSample()) = new HashMap[String,Double]()
    evt1.getEventStrings().foreach { evtString =>
      sampleEventToCount(evt1.getSample())(evtString) = sampleEventToCount(evt1.getSample()).getOrElse(evtString,0.0) + evt1.getSupport.toDouble
      eventToCountAll(evtString) = eventToCountAll.getOrElse(evtString,0.0) + evt1.getSupport.toDouble
    }
  }}

  // now for each sample we saw, normalize by the total coverage, and multiply by the suspected number of cells in the grown
  // organism
  val eventToCountNormalized = new HashMap[String, HashMap[String, Double]]()
  sampleEventToCount.foreach{case(sample,eventToCount) => {
    val sampleTotal = eventToCount.map{case(event,count) => count}.sum
    eventToCountNormalized(sample) = eventToCount.map{case(event,count) => (event,(count / sampleTotal) * numberOfCells)}
  }}

  val eventToCountAllNormalized = new HashMap[String, Double]()
  val allTotal = eventToCountAllNormalized.map{case(e,c) => c}.sum
  eventToCountAll.foreach{case(event,count) => {
    eventToCountAllNormalized(event) = (count / allTotal) * numberOfCells
  }}

  def getEventCount(sample: String,event: String ): Double = {
    if (!(eventToCountNormalized contains sample))
      eventToCountAllNormalized.getOrElse(event,0.001)
    else
      eventToCountNormalized(sample).getOrElse(event,0.01)
  }

  override def getAllSampleEvents(sample: String): mutable.HashMap[String, Double] = sampleEventToCount(sample)
  override def getAllEvents(): mutable.HashMap[String, Double] = eventToCountAll
}
