package main.scala

import beast.evolution.tree.Node
import beast.util.TreeParser
import collection.JavaConverters._
import scala.StringBuilder
import scala.annotation.tailrec
import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer

/**
  * This class takes a root BEAST node and creates a descriptive depth-first tree.  This richer node-type can
  * store all of the annotation data we're interested in, and need for the JSON output
  */
case class RichNode(originalNd: Node,
                    annotations: AnnotationsManager,
                    parent: Option[RichNode],
                    numberOfTargets: Int = 10) extends Ordered[RichNode] {

  // store the orginal node for later if needed
  val originalNode = originalNd


  var name = originalNd.getID
  // handle some basic annotations about the ndoe
  val myAnnotations = annotations.annotationMapping.get(name)

  var distToParent = originalNd.getHeight - (if (parent.isDefined) parent.get.originalNd.getHeight else 0)
  var height = originalNd.getHeight

  var color = "black"
  var sampleName = "UNKNOWN"
  var count = 0
  val freeAnnotations = new mutable.HashMap[String, String]()

  // setup a store for our organ/taxa proportions, which we'll update with counts from the
  // children nodes
  var sampleProportions = new mutable.HashMap[String, Double]()

  // do we have taxa counts to fill in? internal nodes won't have these values but leaves will
  if (myAnnotations.isDefined) {
    sampleName = myAnnotations.get.sample
    color = annotations.cladeMapping(sampleName).color
    count = myAnnotations.get.count
    sampleProportions(sampleName) = myAnnotations.get.proportion
  }

  // explicitly pull out our event string
  val eventString: Option[Array[String]] = if (myAnnotations.isDefined) Some(myAnnotations.get.event.split(annotations.eventSeperator)) else None

  // our parsimony events
  val parsimonyEvents = Array.fill(numberOfTargets)("NONE")

  // now store each of our children and store the events on each child branch we see
  var children = Array[RichNode]()
  var childrenEvents = Array[String]()

  originalNd.getChildren.asScala.foreach { nd => {
    children :+= RichNode(nd, annotations, Some(this))
  }
  }

  // ******************************************************************************************************
  // our member methods
  // ******************************************************************************************************

  // when we change children, we have to change annotations
  def resetChildrenAnnotations(): Unit = {
    childrenEvents = Array[String]()
    sampleProportions = new mutable.HashMap[String, Double]()
    if (myAnnotations.isDefined) {
      sampleProportions(myAnnotations.get.sample) = myAnnotations.get.proportion
    }

    children.foreach { newChild => {
      // get the aggregate children events
      newChild.eventString.foreach(mp => childrenEvents :+= mp.mkString(annotations.eventSeperator))
      newChild.childrenEvents.foreach { mp => childrenEvents :+= mp }

      // add to our existing taxa proportions
      newChild.sampleProportions.foreach {
        case (sample, proportion) => sampleProportions(sample) = sampleProportions.getOrElse(sample, 0.0) + proportion
      }

    }
    }
  }

  /**
    * find the shared edits in nodes underneath this node
    *
    * @return the string of edits in common in our children, if nothing is shared we use '*'
    */
  def sharedEdits(): Array[String] = {
    val storage = Array.fill[mutable.HashSet[String]](numberOfTargets)(mutable.HashSet[String]())
    childrenEvents.foreach { event =>
      event.split(annotations.eventSeperator).zipWithIndex.foreach { case (edit, index) => storage(index) += edit }
    }
    storage.map { st => if (st.size == 1) st.iterator.next else "*" }
  }

  // figure out the shared events
  def parsimonyGenotypeDistance(otherNode: RichNode): Int = {
    parsimonyEvents.zip(otherNode.parsimonyEvents).map { case (evt1, evt2) => if (evt1 == evt2) 0 else 1 }.sum
  }

  // are we wild-type or not?
  def getConsistency(): String = {
    if (this.parsimonyEvents.map { mp => if (mp == "NONE") 0 else 1 }.sum == 0)
      return "WT"
    return "NOTWT"
  }

  // some recursive counting functions
  def countSubNodes(): Int = 1 + children.map { chd => chd.countSubNodes() }.sum

  def countSubUMIs(): Int = {
    val myCount = if (myAnnotations.isDefined) myAnnotations.get.count else 0
    myCount + children.map { chd => chd.countSubUMIs() }.sum
  }

  def countSubProportions(): Double = {
    val myCount = if (myAnnotations.isDefined) myAnnotations.get.proportion else 0.0
    myCount + children.map { chd => chd.countSubUMIs() }.sum
  }

  /**
    * we order RichNodes based on their parsimony allele string
    *
    * @param that the other node to compare to
    **/
  def compare(that: RichNode) = {
    val ret = this.parsimonyEvents.zip(that.parsimonyEvents).map { case (thisEvent, thatEvent) => {
      if (thisEvent == thatEvent)
        0
      else if (thisEvent == "NONE")
        -1
      else if (thatEvent == "NONE")
        1
      else if (thisEvent.compareTo(thatEvent) > 0)
        1
      else
        -1
    }
    }.sum
    //println(this.parsimonyEvents.mkString("_") + " " + that.parsimonyEvents.mkString("_") + ret)
    ret
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
    *
    * @param node           the RichNode to recurse on
    * @param ourBranchJustOrgan the color to set for the branches from this node to the parent
    */
  def assignBranchColors(node: RichNode, ourBranchJustOrgan: String = "false"): Unit = {
    node.freeAnnotations("justOrganSplit") = ourBranchJustOrgan

    // if we have children, check if their node all share the same edit, if so
    // propigate the 'gray' color to them
    if (node.children.size > 0) {
      val counts = node.children.map { child => child.parsimonyEvents.mkString("_") }
      node.children.map { child => {
        val newColor = counts.count(_ == child.parsimonyEvents.mkString("_")) match {
          case x if x > 1 => "true"
          case _ => "false"
        }
        assignBranchColors(child, newColor)
      }}
    }
  }

  /**
    * given the results of the parsimony run, figure out the genotypes of each node in turn.
    * this means accumulating events from the root outwards, as single nodes only have the changes
    * compared to the previous node
    *
    * @param rootNode the root node, which we assume is all NONE in camin-sokal parsimony
    * @param parser   the results from the parsimony run
    */
  def applyParsimonyGenotypes(rootNode: RichNode, parser: MixParser, numberOfTargets: Int = 10): Unit = {
    // lookup each link between a subnode and the root, and assign it's genotypes recursively
    rootNode.children.foreach { newChild => recAssignGentoypes(rootNode, newChild, parser) }
  }

  /**
    * recursively walk down the tree assigning genotypes to each of the progeny nodes as we go.  We have to do
    * this as the parsimony output is recursive; they define each nodes identity in terms of it's parents
    * genotype with edit's overlaid
    *
    * @param parent the parent of this node
    * @param child  this node, the child
    * @param parser the parser which contains all of the info about genotypes, etc
    */
  def recAssignGentoypes(parent: RichNode, child: RichNode, parser: MixParser): Unit = {
    // copy the parents genotype over to the child
    parent.parsimonyEvents.zipWithIndex.foreach { case (edit, index) => child.parsimonyEvents(index) = edit }

    // find the link from out parent node -- there should only be one edge leading to this node ever
    val link = parser.lookupTos(child.name)

    // make a list of the events that we're adding
    link.chars.zipWithIndex.foreach { case (change, index) => change match {
      case ('.') => {
        /* do nothing */
      }
      case ('1') => {
        val event = parser.numberToEvent(index + 1) // our first position is an edit, not NONE
        parser.eventToSites(event).foreach { site => {
          // check to make sure we're not conflicting and overwriting an already edited site
          if (child.parsimonyEvents(site) != "NONE")
            println("WARNING: Conflict at site " + site + " for parent " + parent.name + " for child " + child.name + " event " + event)
          child.parsimonyEvents(site) = event
        }
        }
      }
    }
    }

    // now for each of the children of this node, recursively assign genotypes
    child.children.foreach { newChild => recAssignGentoypes(child, newChild, parser) }
  }

  /**
    * our internal names come off the tree without names, even though the genotypes are assigned a name.  What we have
    * to do here is traverse to the leaves, which do have names, and walk backwards assigning parent names from known
    * leaf child->parent output
    *
    * @param node   the node to assign a name to
    * @param parser the mix output parser with all of the parent and child relationships
    * @return the string of the top node (lazy, this is used to recursively fill in names on the way
    *         back up)
    */
  def recAssignNames(node: RichNode, parser: MixParser): String = {
    // if we have a leaf -- where there are no children -- assign the name
    if (node.children.size == 0) {
      val edge = parser.lookupTos(node.name)
      return edge.from
    } else {
      val names = node.children.map { case (nd) => recAssignNames(nd, parser) }.toSet.toList
      if (names.size != 1)
        throw new IllegalStateException("Unable to assign the name for node with children " + names.mkString(","))
      node.name = names(0)
      val edge = parser.lookupTos(names(0))
      return edge.from
    }
  }

  /**
    * sort the nodes by their allele string
    *
    * @param node the node to start sorting on
    */
  def reorderChildrenByAlleleString(node: RichNode): Unit = {

      if (node.children.size > 0)
        print(node.children(0).count + " -> ")
      scala.util.Sorting.quickSort(node.children)
      if (node.children.size > 0)
        println(node.children(0).count)
      node.children.foreach { case (nd) => reorderChildrenByAlleleString(nd) }
  }

  /**
    * apply a function to all nodes in depth first order, storing the results in the annotations
    *
    * @param node the node to start sorting on
    */
  def applyFunction(node: RichNode, func: RichNode => Tuple2[String,String]): Unit = {
    val res = func(node)
    node.freeAnnotations(res._1) = res._2
    node.children.foreach { case (nd) => applyFunction(nd,func) }
  }

  /**
    * check that our nodes are assigned consistent node identities between the parsimony and known annotations
    *
    * @param node   the node
    * @param parser the mix parser
    */
  def recCheckNodeConsistency(node: RichNode, parser: MixParser): Unit = {
    // if we have a leaf -- where there are no children -- assign the name
    if (node.children.size == 0 && node.eventString.isDefined) {
      val differences = node.eventString.get.zip(node.parsimonyEvents).map { case (evt1, evt2) => if (evt1 == evt2) 0 else 1 }.sum
      if (differences > 0) {
        println("FAIL " + node.eventString.get.mkString(",") + " - " + node.parsimonyEvents.mkString(","))
      }
    } else {
      node.children.map { case (nd) => recCheckNodeConsistency(nd, parser) }.toSet.toList
    }
  }

  /**
    * find the depth of the deepest node on the tree
    *
    * @param node the node
    */
  def maxHeight(node: RichNode): Double = {
    if (node.children.size == 0)
      node.height
    else
      Math.max(node.height, node.children.foldLeft(0.0)((a, b) => Math.max(a, maxHeight(b))))
  }


  /**
    * recurrsively output the current tree as a JSON string
    *
    * @param node         the node to start at
    * @param parent       it's parent node, optional (for recussion)
    * @param distToRoot   the distance to the root node, as we convert to a zero or one system
    * @param noParentName if we don't have a parent, what should we call the previous node
    * @return a JSON string representation
    */
  def toJSONOutput(node: RichNode, parent: Option[RichNode], distToRoot: Double, noParentName: String = "null"): String = {
    val outputString = new ArrayBuffer[String]()
    outputString += RichNode.toJSON("name", node.name)
    outputString += RichNode.toJSON("parent", if (parent.isDefined) parent.get.name else noParentName)
    outputString += RichNode.toJSON("length", 1)
    outputString += RichNode.toJSON("rootDist", distToRoot)
    outputString += RichNode.toJSON("cladeTotal", node.countSubNodes())
    outputString += RichNode.toJSON("totatSubNodes", node.countSubProportions())
    outputString += RichNode.toJSON("color", node.color)
    outputString += RichNode.toJSON("sample", node.sampleName)
    node.freeAnnotations.foreach { case (key, value) =>
      outputString += RichNode.toJSON(key, value)
    }
    val sampleTot = if (node.annotations.sampleTotals contains node.sampleName) node.annotations.sampleTotals(node.sampleName) else 0
    outputString += RichNode.toJSON("organCountsMax", sampleTot)
    outputString += RichNode.toJSON("cladeTotal", node.count)
    outputString += RichNode.toJSON("max_organ_prop", if (sampleTot.toDouble > 0) node.count.toDouble / sampleTot.toDouble else 0.0)
    outputString += RichNode.toJSON("event", node.parsimonyEvents.mkString(node.annotations.eventSeperator))
    outputString += RichNode.toJSON("commonEvent", node.sharedEdits.mkString(node.annotations.eventSeperator))
    outputString += "\"organProportions\": { \n"
    outputString += node.sampleProportions.map { case (sample, prop) => RichNode.toJSON(sample, prop, terminator = "") }.mkString(",\n")
    outputString += "\n},\n"
    if (node.children.size > 0) {
      outputString += "\"children\": [{\n"
      outputString += node.children.map { child => RichNode.toJSONOutput(child, Some(node), distToRoot + 1.0) }.mkString("},\n{\n")
      outputString += "}],\n"
    }
    outputString += RichNode.toJSON("consistency", node.getConsistency, terminator = "\n")
    outputString.toArray.mkString("")
  }


  def toJSON(name: String, simple: Any, terminator: String = ",\n"): String = simple match {
    case x: String => "\"" + name + "\" : \"" + x + "\"" + terminator
    case x: Int => "\"" + name + "\" : " + x + terminator
    case x: Double => "\"" + name + "\" : " + x + terminator
    case x: Float => "\"" + name + "\" : " + x + terminator
    case _ => throw new IllegalStateException("We don't know what to do ")
  }
}