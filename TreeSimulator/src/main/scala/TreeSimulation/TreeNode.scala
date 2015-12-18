package main.scala.TreeSimulation

import java.io.PrintWriter
import akka.actor.Actor

import scala.collection.mutable._


/**
 * This class represents wat we know at each node in our evolutionary tree:
 * 1) the 'model': this is misnamed a bit, but it contains all the information we know about the current node structure
 * 2) the current depth of the tree at this point
 */
case class TreeNode(model: GlobalSiteModel,
                    depth: Int,
                    previousLogLikelihood: Double,
                    eventString: Array[String],
                     mu: Double) {

  def terminal(): Boolean = model.leaves.size < 1 || !model.remainingSamples

  var logLikelihood = previousLogLikelihood

  var newickString = eventString.mkString("_") + "_" + depth // + ":" + displayLL(previousLogLikelihood)
  var maxDepth = depth

  var leftNode : Option[TreeNode] = None
  var rightNode : Option[TreeNode] = None

  // get a distribution that represent the chance we mutate

  if (!terminal) {
  // draw the next event
    try {
      val site = model.siteLevelDist.drawFromDistribution()._1
      // println("Drawing from site " + site + " with counts " + model.eventDistributions(site).eventDistributions.size + " sizes are " +
      //  model.eventDistributions.map{dst => dst.eventDistributions.size}.mkString(","))

      val eventIndex = model.eventDistributions(site).drawFromDistribution()
      val event = model.eventDistributions(site).drawFromDistribution()

      // now split the tree into events that contain the event at the specific site (left) and those that don't (right)
      val leftEvents = GlobalSiteModel(model.leaves.filter { leaf => leaf.eventMap contains event._1 }, model.eventSize, model.alreadySpiltOn ++ Array[String](event._1))
      val rightEvents = GlobalSiteModel(model.leaves.filter { leaf => !(leaf.eventMap contains event._1) }, model.eventSize, model.alreadySpiltOn ++ Array[String](event._1))

      // println("spliting with left size " + leftEvents.leaves.size + " right events size " + rightEvents.leaves.size + " on event " + event)

      val leftString = eventString.clone()
      leftString(site) = event._1
      val rightString = eventString.clone()

      leftNode = Some(TreeNode(leftEvents, depth + 1, math.log10(event._2), leftString))
      rightNode = Some(TreeNode(rightEvents, depth + 1, math.log10(1), rightString))

      //currentLogLikelihood += leftNode.currentLogLikelihood
      //currentLogLikelihood += rightNode.currentLogLikelihood

      newickString = "(" + leftNode.get.newickString + "," + rightNode.get.newickString + ")" // ," + newickString + ")"
      logLikelihood += leftNode.get.logLikelihood
      logLikelihood += rightNode.get.logLikelihood
      maxDepth = math.max(leftNode.get.maxDepth, rightNode.get.maxDepth)

    } catch {
      case e: Exception => {
        println("Unable to split node at " + model.leaves.size)
        throw (e)
      }
    }
  }



  def displayLL(logLikelihood: Double) = math.abs(logLikelihood * -1.0).toString

  def writeAnnotations(eventToKnown: HashMap[String,String],
                       outputAnnotations: PrintWriter): Unit = {
    if (eventToKnown contains eventString.mkString("_"))
      outputAnnotations.write(eventString.mkString("_") + "_" + depth + "\t" + eventToKnown(eventString.mkString("_")) + "\n")
    else
      outputAnnotations.write(eventString.mkString("_") + "_" + depth + "\tinternal\n")
    if (leftNode.isDefined)
      leftNode.get.writeAnnotations(eventToKnown,outputAnnotations)
    if (rightNode.isDefined)
      rightNode.get.writeAnnotations(eventToKnown,outputAnnotations)
  }
}

