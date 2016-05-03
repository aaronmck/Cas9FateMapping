package main.scala

import java.io.{PrintWriter, File}

import beast.evolution.tree.Node
import beast.util.TreeParser

import scala.io.Source

/**
  * given all the parts we need to make a rich JSON version
  * of a parsimony tree from the MIX program, make the enhanced tree
  * and then collapse out parsimony nodes that are
  */
class ParsimonyProcessor(mixTrees: File,
                         mixOutput: File,
                         annotations: File,
                         sampleToClade: File,
                         eventsToNumbers: File,
                         outputFile: File) {

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

  // count nodes before
  println("nodes " + rootNode.countSubNodes())

  // make a collapser of nodes
  ParsimonyCollapser.checkCollapseNodes(rootNode)

  // the updated numbers
  println("nodes " + rootNode.countSubNodes())

  // now output the adjusted tree
  val output = new PrintWriter(outputFile.getAbsolutePath)
  output.write("[{\n")
  val jsonString = RichNode.toJSONOutput(rootNode, None)
  output.write("}]\n")
  output.close()

}
