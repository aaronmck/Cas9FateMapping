package main.scala


import beast.util._
import beast.evolution.tree._
import scala.collection.JavaConversions._
import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._
import scala.util.Random

/**
 * created by aaronmck on 2/13/14
 *
 * Copyright (c) 2014, aaronmck
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2.  Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */
case class TreeConfig(inputTree: File = new File(Main.NOTAREALFILENAME),
                  outputSubTrees: File = new File(Main.NOTAREALFILENAME),
                  mappingFile: File = new File(Main.NOTAREALFILENAME))



object Main extends App {
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[TreeConfig]("UMIMerge") {
    head("TreeUtils", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("inputTree") required() valueName ("<file>") action { (x, c) => c.copy(inputTree = x) } text ("first read file ")
    opt[File]("outputSubTrees") required() valueName ("<file>") action { (x, c) => c.copy(outputSubTrees = x) } text ("second reads file")
    opt[File]("mappingFile") required() valueName ("<file>") action { (x, c) => c.copy(mappingFile = x) } text ("mapping node names to ")


    // some general command-line setup stuff
    note("processes reads with UMIs into merged reads\n")
    help("help") text ("prints the usage information you see here")
  }

  // *********************************** Run *******************************************************
  parser.parse(args, TreeConfig()) map {
    config: TreeConfig => {

      // load up the annotation manager
      println("Loading annotation manager from your annotation file...")
      val annotManager = AnnotationManager(config.mappingFile.getAbsolutePath)

      println("finding X% pure subtrees from your data...")
      findPureBranches(config.inputTree.getAbsolutePath,0.95, 0,annotManager, config.outputSubTrees.getAbsolutePath, 10000 )

      /*
      val cutHeight = 60.0

      println(root.getHeight)
      getAverageLeafDistanceToRoot(root)

      val output  = new PrintWriter(config.outputSubTrees)
      getSubtreesAtSetHeight(root,output,cutHeight)
      //getSubtreesAtSetDivisions(root,output,5)
      output.close()
      */
    }
  } getOrElse {
    println("Unable to parse the command line arguments you passed in, please check that your parameters are correct")
  }


  /**
   * manages our annotations.  Put this all in one place so we can make a scrambled version of it later
    *
    * @param mappingFileName the input file
   */
  case class AnnotationManager(mappingFileName: String) {

    // store leaf node names to their group and read
    val nodeNameToGroup = new HashMap[String,String]()
    val nodeNameToReads = new HashMap[String,Int]()

    // store the relationship between the species / group name to the node IDs
    val nodeSamplesToReads = new HashMap[String,Int]()

    // find the total number of reads for each embryo
    val tagToReadCount = new HashMap[String,Int]()

    val mappingFile = Source.fromFile(mappingFileName).getLines()
    val header = mappingFile.next().split("\t").zipWithIndex.map{case(tk,index) => (tk,index)}.toMap

    if (!(header contains "taxa")) throw new IllegalStateException("Unable to find taxa in your header")
    if (!(header contains "sample")) throw new IllegalStateException("Unable to find sample in your header")
    if (!(header contains "count")) throw new IllegalStateException("Unable to find count in your header")
    if (!(header contains "eventString")) throw new IllegalStateException("Unable to find eventString in your header")

    mappingFile.foreach{line => {
      val sp = line.split("\t")
      nodeNameToGroup(sp(0)) = sp(header("sample"))
      nodeNameToReads(sp(0)) = sp(header("count")).toInt
      nodeSamplesToReads(sp(header("sample"))) = nodeSamplesToReads.getOrElse(sp(header("sample")),0) + nodeNameToReads(sp(0))
    }}

    nodeSamplesToReads.foreach{case(tag,reads) => {
      println(tag + "\t" + reads)
      tagToReadCount(tag) = tagToReadCount.getOrElse(tag,0) + reads
    }}
  }


  /**
   * find the average distance you have to travel when going from the root node to any leaf
    *
    * @param rootNode the root node to consider
   */
  def getAverageLeafDistanceToRoot(rootNode: Node): Unit = {
    var heights = ArrayBuilder.make[Double]
    var oldNodes = new ArrayBuffer[Node]()
    oldNodes += rootNode

    rootNode.getAllLeafNodes.foreach { child => {
      //println(child.getID + "\t" + child.getHeight + "\t" + child.getLength)
      var nd = child
      var height = 0.0
      var edges = 0
      while (!nd.isRoot) {
        height += nd.getLength
        edges += 1
        nd = nd.getParent
      }
      heights += height
    }
    }
    println(heights.result.sum / heights.result.size)
  }


  /**
   *
   * @param treeFile the input tree
   * @param propCutoff the proportion of leaves that need to be the same label
   * @param minNodes the min number of nodes required to find a valid subtree
   * @param annotManager a mapping of the node labels to the class
   * @param output a place to print successful subtrees
   */
  def findPureBranches(treeFile: String,
                       propCutoff: Double,
                       minNodes: Int,
                       annotManager: AnnotationManager,
                       output: String,
                       numberOfSimulatedTrees: Int): Unit = {

    val isLabeled = true

    val outputFile = new PrintWriter(output)

    outputFile.write("isSimulation\tmaxProp\ttarget.reads\ttarget.nodes\ttotalEmbryoReadProp\tpropCutoff\tmaxLabel\tcounts\tnewick\n")

    val treeParser = new TreeParser(Source.fromFile(treeFile).getLines().mkString(""), false, true, isLabeled, 1)
    purityOverTree(treeParser.getRoot,propCutoff,minNodes,annotManager,outputFile,false)

    // let's do 100 simulations
    /*(0 until numberOfSimulatedTrees)foreach {index => {
      val treeParser = new TreeParser(Source.fromFile(treeFile).getLines().mkString(""), false, true, isLabeled, 1)

      val root = treeParser.getRoot
      val newRoot = permuteTreeLabels(root)
      purityOverTree(newRoot,propCutoff,minNodes,annotManager,outputFile,true)
    }}
    */
    outputFile.close()
  }


  /**
   * find the purity of the subtree
    *
    * @param root the root node
   * @param propCutoff the cutoff we need to call a subtree pure
   * @param minNodes the minimum number of nodes to report for a pure subtree
   * @param annotManager the annotations
   * @param outputFile the output file we're writing to
   * @param isSimulationTree is this a simulation or a real tree?
   */
  def purityOverTree(root: Node,
                     propCutoff: Double,
                     minNodes: Int,
                     annotManager: AnnotationManager,
                     outputFile: PrintWriter,
                     isSimulationTree: Boolean) = {

    var stillConsideringNodes = new ArrayBuffer[Node]()
    stillConsideringNodes += root

    while (stillConsideringNodes.size > 0) {
      var newNodesToReplace = new ArrayBuffer[Node]()

      stillConsideringNodes.toArray.foreach { nd => {
        val purity = determineSubtreePurity(nd, propCutoff, minNodes, annotManager, outputFile, isSimulationTree)

        if (!(purity isDefined))
          nd.getChildren.foreach { child => newNodesToReplace += child }

      }
      }
      stillConsideringNodes = newNodesToReplace
    }

  }

  /**
   * given a subtree, figure out how 'pure' it is -- how many of the nodes underneath
   * have the same label in our nodeNameToGroup tabl
    *
    * @param node the node to traverse over
   * @param propCutoff the proportion of leaves that need to be the same label
   * @param minNodes the min number of nodes required to find a valid subtree
   * @param annotManager all of our annotations for this tree
   * @param output a place to print successful subtrees
   * @return which subtree label was the max, the size of the leaf nodes below this node, and the proportion of the max label
   */
  def determineSubtreePurity (node: Node,
                              propCutoff: Double,
                              minNodes: Int,
                              annotManager: AnnotationManager,
                              output: PrintWriter,
                              isSimulationTree: Boolean) : Option[Tuple3[String, Int, Double]] = {

    val tagCounts = new HashMap[String,Int]()
    val tagReads = new HashMap[String,Int]()

    var sum = 0

    // println(node.getAllLeafNodes.size)

    node.getAllLeafNodes.foreach{leaf =>
      tagCounts(annotManager.nodeNameToGroup(leaf.getID)) = tagCounts.getOrElse(annotManager.nodeNameToGroup(leaf.getID),0) + 1
      tagReads(annotManager.nodeNameToGroup(leaf.getID)) =  tagReads.getOrElse(annotManager.nodeNameToGroup(leaf.getID),0) + annotManager.nodeNameToReads(leaf.getID)
      sum += 1
    }

    //tagCounts.foreach{case(key,value) => println(key + "\t" + value)}

    var maxProp = 0.0
    var maxLabel = ""
    tagCounts.foreach{case(group,count) => {
      val prop = count.toDouble/sum.toDouble
      if (prop > maxProp) {
        maxProp = prop
        maxLabel = group
      }
    }}

    if (maxProp >= propCutoff) {
      //println("Found maxlabel")
      val totalEmbryoReadProp = tagReads(maxLabel).toDouble/annotManager.tagToReadCount(maxLabel).toDouble
      output.write(isSimulationTree + "\t" + maxProp + "\t" + tagReads(maxLabel) + "\t" + tagCounts(maxLabel) + "\t" + totalEmbryoReadProp + "\t")
      output.write(maxProp + "\t" + maxLabel + "\t" + tagCounts.map{case(key,value) => key + ":" + value}.mkString(","))
      output.write("\t" + node.toNewick + "\n")
    }

    if (maxProp < propCutoff)
      return None

    return Some(maxLabel,node.getAllLeafNodes.size,maxProp)
  }


  /**
   * split out a subtree as a set depth
    *
    * @param rootNode the node to start at
   * @param outputNewick the output file to write to
   * @param divisions how many cell divisions should we split at?
   */
  def getSubtreesAtSetDivisions(rootNode: Node, outputNewick: PrintWriter, divisions: Int): Unit = {
    var stillConsideringNodes = new ArrayBuffer[Node]()
    stillConsideringNodes += rootNode
    var depth = 0
    while (depth < divisions) {
      println("round...")
      depth += 1
      var newNodesToReplace = new ArrayBuffer[Node]()
      stillConsideringNodes.toArray.foreach { nd => {
        nd.getChildren.foreach { child => newNodesToReplace += child }
      }
      }
      stillConsideringNodes = newNodesToReplace
    }
    stillConsideringNodes.foreach { child => outputNewick.write(child.toNewick + "\n") }
  }

  /**
   * find all the set trees as a specific height
    *
    * @param rootNode the node to start from
   * @param outputNewick the output file
   * @param height the height to cut at
   */
  def getSubtreesAtSetHeight(rootNode: Node, outputNewick: PrintWriter, height: Double): Unit = {
    var stillConsideringNodes = new ArrayBuffer[Node]()
    stillConsideringNodes += rootNode

    while (stillConsideringNodes.size > 0) {
      println("round...")
      var newNodesToReplace = new ArrayBuffer[Node]()
      stillConsideringNodes.toArray.foreach{nd => {
        if (nd.getHeight < height && nd.getChildren.size >= 1) {
          outputNewick.write(nd.toNewick + "\n")
          println("MIN HEIght reached ")
        } else if (nd.getLength < 0) {
          println("NEG BRANCH HIT")
          nd.getChildren.foreach { child => outputNewick.write(child.toNewick + "\n") }

        } else
          nd.getChildren.foreach{child => newNodesToReplace += child}
      }}
      stillConsideringNodes = newNodesToReplace
    }
  }


  /**
   * permute a tree's labels at the leaves.  I do this by iterating through each leaf, picking another random
   * leaf among all leaves, and swapping labels
    *
    * @param root the root node
   * @return the same root node, with leaves transformed
   */
  def permuteTreeLabels(root: Node): Node = {
    val rand = new Random()

    // for each leaf in the tree, pick a random other leaf, and swap their labels
    val leaves = root.getAllLeafNodes

    leaves.foreach(leaf => {
      val randomLeaf = rand.nextInt(leaves.size())
      val tmpName = leaves(randomLeaf).getID
      leaves(randomLeaf).setID(leaf.getID())
      leaf.setID(tmpName)
    })
    return root
  }

}
