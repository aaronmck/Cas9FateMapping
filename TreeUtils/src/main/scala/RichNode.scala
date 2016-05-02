package main.scala

import beast.evolution.tree.Node
import beast.util.TreeParser
import collection.JavaConverters._
import scala.annotation.tailrec
import scala.collection.mutable

/**
  * This class takes a root BEAST node and creates a descriptive depth-first tree.  This richer node-type can
  * store all of the annotation data we're interested in, and need for the JSON output
  */
case class RichNode(originalNd: Node,
                    annotations: AnnotationsManager,
                    parent: Option[RichNode],
                    numberOfTargets: Int = 10) {

  val name = originalNd.getID

  val originalNode = originalNd

  val myAnnotations = annotations.annotationMapping.get(name)

  // setup a store for our organ/taxa proportions, which we'll update with counts from the
  // children nodes
  val taxaProportions = new mutable.HashMap[String,Double]()

  // do we have taxa counts to fill in? internal nodes won't have these values but leaves will
  if (myAnnotations.isDefined) {
    taxaProportions(myAnnotations.get.taxa) = myAnnotations.get.proportion
  }

  // explicitly pull out our event string
  val eventString = if (myAnnotations.isDefined) Some(myAnnotations.get.event) else None

  // now store each of our children
  var children = Array[RichNode]()

  // store the events on each child branch we see
  var childrenEvents = Array[Array[String]]()
  var parsimonyEvents = Array[String]()

  originalNd.getChildren.asScala.foreach{nd => {
    val newChild = RichNode(nd, annotations, Some(this))

    // get the aggregate children events
    childrenEvents :+= newChild.childrenEvents.flatMap(x => x)

    // add to our existing taxa proportions
    newChild.taxaProportions.foreach{
      case(taxa,proportion) => taxaProportions(taxa) = taxaProportions.getOrElse(taxa,0.0) + proportion
    }

    // finally add our child to the array of children
    children :+= newChild
  }}

}

object RichNode {
  def toRichTree(inputTree: TreeParser, annotationsManager: AnnotationsManager): RichNode = {
    return RichNode(inputTree.getRoot, annotationsManager, None)
  }

  /**
    * given the results of the parsimony run, figure out the genotypes of each node in turn.
    * this means accumulating events from the root outwards, as single nodes only have the changes
    * compared to the previous node
    *
    * @param rootNode the root node, which we assume is all NONE in camin-sokal parsimony
    * @param parser the results from the parsimony run
    */
  def applyParsimonyGenotypes(rootNode: RichNode, parser: MixParser, numberOfTargets: Int = 10): Unit = {
    // assign the root node to the default NONEs
    rootNode.parsimonyEvents = Array[String]("NONE" * numberOfTargets)

    // now lookup each link between a subnode and the root, and assign it's genotypes recursively
    rootNode.children.foreach{newChild => recAssignGentoypes(rootNode,newChild,parser)}
  }

  /**
    * recursively walk down the tree assigning genotypes to each of the progeny nodes as we go
 *
    * @param parent the parent of this node
    * @param child this node, the child
    * @param parser the parser which contains all of the info about genotypes, etc
    */
  def recAssignGentoypes(parent: RichNode, child: RichNode, parser: MixParser): Unit = {
    // copy the parents genotype over to the child
    child.parsimonyEvents = parent.parsimonyEvents.clone()

    // find the link from out parent node -- there should only be one edge leading to this node ever
    val link = parser.lookupTos(child.name)

    // make a list of the events that we're adding
    link.chars.zipWithIndex.map{case(change,index) => change match {
      case ('.') => {/* do nothing */ }
      case ('1') => {
        val event = parser.numberToEvent(index + 1) // our first position is an edit, not NONE
        parser.eventToSites(event).foreach{site => {
          // check to make sure we're not conflicting and overwriting an already edited site
          if (child.parsimonyEvents(site) != "NONE")
            println("WARNING: Conflict at site " + site + " for parent " + parent.name + " for child " + child.name + " event " + event)
          child.parsimonyEvents(site) = event
        }}
      }
    }}

    // now for each of the children of this node, recursively assign genotypes
    child.children.foreach{newChild => recAssignGentoypes(child,newChild,parser)}
  }


}
