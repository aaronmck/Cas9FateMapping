package main.scala

import java.io.File

import beast.evolution.tree.Node
import beast.util.TreeParser

import scala.io.Source

/**
  * given all the parts we need to make a rich JSON version
  * of a parsimony tree from the MIX program
  */
class ParsimonyProcessor(mixTrees: File, mixOutput: File, annotations: File, sampleToClade: File, eventsToNumbers: File) {

  // find the best tree from the mix output
  val bestTreeContainer = BestTree(mixTrees)

  // parse out the annotations from the mix (PHYLIP) output
  val mixParser = new MixParser(mixOutput.getAbsolutePath, eventsToNumbers.getAbsolutePath, bestTreeContainer.maxIndex)

  // load our tree
  val treeParser = new TreeParser(bestTreeContainer.bestTreeString, false, true, true, 1)

  // load up any annotations we have
  val annotationMapping = new AnnotationsManager(annotations, sampleToClade)

  // traverse the nodes and add names to any internal nodes without names
  val rootNode = RichNode(treeParser.getRoot,annotationMapping, None)

  // reassign the names
  val rootName = RichNode.recAssignNames(rootNode, mixParser)

  // now apply the parsimony results to the root of the tree (recursively walking down the nodes)
  RichNode.applyParsimonyGenotypes(rootNode, mixParser)

  // check that the nodes we assigned are consistent
  RichNode.recCheckNodeConsistency(rootNode, mixParser)

  // make a collapser of nodes
  val collapser = new ParsimonyCollapser(rootNode,mixParser)
}
