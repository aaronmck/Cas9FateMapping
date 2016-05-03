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

  var name = originalNd.getID

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
  val eventString: Option[Array[String]] = if (myAnnotations.isDefined) Some(myAnnotations.get.event.split(annotations.eventSeperator)) else None

  // now store each of our children
  var children = Array[RichNode]()

  // store the events on each child branch we see
  var childrenEvents = Array[String]()
  val parsimonyEvents = Array.fill(numberOfTargets)("NONE")

  originalNd.getChildren.asScala.foreach{nd => {
    val newChild = RichNode(nd, annotations, Some(this))

    // get the aggregate children events
    newChild.eventString.foreach(mp => childrenEvents :+= mp.mkString(annotations.eventSeperator))
    newChild.childrenEvents.foreach{mp => childrenEvents :+= mp}

    // add to our existing taxa proportions
    newChild.taxaProportions.foreach{
      case(taxa,proportion) => taxaProportions(taxa) = taxaProportions.getOrElse(taxa,0.0) + proportion
    }

    // finally add our child to the array of children
    children :+= newChild
  }}

  def isCollapsible(): Option[Tuple3[String,Int,RichNode]] = {
    if (childrenEvents.size > 1) {
      if (childrenEvents.toSet.size == 1) {
        return Some(new Tuple3[String,Int,RichNode](childrenEvents.toSet.mkString(","),childrenEvents.size,this))
      }
    }
    return None
  }
}

/**
  * We use the static methods here to transform nodes after they're created.
  */
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
    // lookup each link between a subnode and the root, and assign it's genotypes recursively
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
    parent.parsimonyEvents.zipWithIndex.foreach{case(edit,index) => child.parsimonyEvents(index) = edit}

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

  def recAssignNames(node: RichNode, parser: MixParser): String = {
    // if we have a leaf -- where there are no children -- assign the name
    if (node.children.size == 0) {
      val edge = parser.lookupTos(node.name)
      return edge.from
    } else {
      val names = node.children.map{case(nd) => recAssignNames(nd,parser)}.toSet.toList
      if (names.size != 1)
        throw new IllegalStateException("Unable to assign the name for node with children " + names.mkString(","))
      node.name = names(0)
      val edge = parser.lookupTos(names(0))
      return edge.from
    }
  }

  /**
    * check that our nodes are assigned consistent node identities between the parsimony and known annotations
    *
    * @param node the node
    * @param parser the mix parser
    */
  def recCheckNodeConsistency(node: RichNode, parser: MixParser): Unit = {
    // if we have a leaf -- where there are no children -- assign the name
    if (node.children.size == 0 && node.eventString.isDefined) {
      val differences = node.eventString.get.zip(node.parsimonyEvents).map{case(evt1,evt2) => if (evt1 == evt2) 0 else 1}.sum
      if (differences > 0) {
        println("FAIL " + node.eventString.get.mkString(",") + " - " + node.parsimonyEvents.mkString(","))
      }
    } else {
      node.children.map{case(nd) => recCheckNodeConsistency(nd,parser)}.toSet.toList
    }
  }


}
