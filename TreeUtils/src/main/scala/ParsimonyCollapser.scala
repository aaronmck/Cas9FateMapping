package main.scala

import scala.collection.mutable

/**
  * Created by aaronmck on 5/2/16.
  */
object ParsimonyCollapser {
  /**
    * @param node the rich node
    */
  def checkCollapseNodes(node: RichNode) {
    var newChildren = Array[RichNode]()

    node.children.foreach{case(child) => {
      if (child.children.size > 0 && child.parsimonyGenotypeDistance(node) == 0) {
        // reconnect all the children-of-the-child nodes
        child.children.foreach(chd => newChildren :+= chd)
      } else {
        newChildren :+= child
      }
    }}

    // the order here is very important -- update child list, have them update their children lists, and finally update annotations recursively
    node.children = newChildren
    node.children.foreach{chd => checkCollapseNodes(chd)}
    node.resetChildrenAnnotations()
  } 
}
