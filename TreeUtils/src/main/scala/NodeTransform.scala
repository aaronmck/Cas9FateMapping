package main.scala

import beast.evolution.tree.Node

import scala.annotation.tailrec
import collection.JavaConverters._

/**
  *
  */
object NodeTransform {

  /**
    *
    * @param rootNode the root node to start assigning names to
    */
  def renameInternalNodes(rootNode: Node): Unit = {
    val allNodes = rootNode.getAllChildNodes.asScala
    var lastInternalNumber = 0

    allNodes.foreach{nd => if (nd.getID == "") {
      nd.setID("internal_" + lastInternalNumber.toString)
      lastInternalNumber += 1
    }}
  }

  /*
  def mergeTargetsOneSide(events):
    events_to_targets = [set() for x in range(0,number_of_targets)]
    for evtInd,hmid in enumerate(events):
        for index, event in enumerate(hmid.split(split_token)):
            if evtInd == 0:
                events_to_targets[index].add(event)
            else:
                events_to_targets[index] = events_to_targets[index].intersection(set([event]))
    return events_to_targets
   */

}
