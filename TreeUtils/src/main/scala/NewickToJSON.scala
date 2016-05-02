package main.scala

import beast.util._
import beast.evolution.tree._
import scala.collection.JavaConversions._
import scala.io._
import java.io._
import scala.collection.mutable.{ArrayBuffer, HashMap}
import scala.sys.process._
import java.util.zip._
import scala.util.Random
/*
case class NewickConfig(inputTree: File = new File(NewickToJSON.NOTAREALFILENAME),
                        outputJSON: File = new File(NewickToJSON.NOTAREALFILENAME),
                        annotations: Seq[String] = Seq(),
                        jars: Seq[File] = Seq())
*/

/**
  * Created by aaronmck on 2/1/16.
  */
object NewickToJSON { //extends App {
  /*
  val NOTAREALFILENAME = "/0192348102jr10234712930h8j19p0hjf129-348h512935"
  // please don't make a file with this name
  val NOTAREALFILE = new File(NOTAREALFILENAME)

  // parse the command line arguments
  val parser = new scopt.OptionParser[NewickConfig]("NewickToJSON") {
    head("NewickToJSON", "1.0")

    // *********************************** Inputs *******************************************************
    opt[File]("inputTree") required() valueName ("<file>") action { (x, c) => c.copy(inputTree = x) } text ("newick tree")
    opt[File]("outputJSON") required() valueName ("<file>") action { (x, c) => c.copy(outputJSON = x) } text ("the output json file")
    opt[Seq[String]]("annotation") valueName ("<jar1>,<jar2>...") action { (x, c) => c.copy(annotations = x) } text ("the annotation input files")
    //opt[Seq[File]]('j', "jars") valueName("<jar1>,<jar2>...") action { (x,c) =>
    //  c.copy(jars = x) } text("jars to include")

    // some general command-line setup stuff
    note("processes a newick file into a JSON richly annotated version\n")
    help("help") text ("prints the usage information you see here")
  }

  // *********************************** Run *******************************************************
  parser.parse(args, NewickConfig()) match {
    case Some(newickConfig) => {
      treeToJSON(newickConfig.inputTree.getAbsolutePath, newickConfig.annotations.toList, newickConfig.outputJSON.getAbsolutePath)
    }
    case None => {
      println("Unable to parse the command line arguments you passed in, please check that your parameters are correct")
    }
  }

  /**
    *
    * @param treeFile the input tree
    * @param annotationFiles an array of annotation files to richly annotate our tree nodes with
    * @param outputFile the output destination to write to
    */
  def treeToJSON(treeFile: String,
                 annotationFiles: List[String],
                 outputFile: String): Unit = {

    val nodeNameCol = "taxa"

    // load up the tree file
    // ------------------------------------------------------------------------------------------------------------------------
    val treeParser = new TreeParser(Source.fromFile(treeFile).getLines().mkString(""), false, true, true, 1)

    // now make a mapping from each of the input file sample annotations, to the annotation name, to the value
    val nodeToAnnotationToValue = new HashMap[String,HashMap[String,String]]()

    // process the annotations files, filling in the mapping from sample to annotation to value
    // ------------------------------------------------------------------------------------------------------------------------
    annotationFiles.foreach{annotationFile => {
      val lines = Source.fromFile(annotationFile).getLines()
      val header = lines.next().split("\t")
      val splitHeader = header.map{tk => (tk,header.indexOf(tk))}.toMap

      var loadedAnnotations = 0
      lines.foreach{line => {
        val sp = line.split("\t")
        val sample = sp(splitHeader(nodeNameCol))
        if (!nodeToAnnotationToValue.contains(sample))
          nodeToAnnotationToValue(sample) = new HashMap[String,String]()

        splitHeader.foreach{case(headerToken,index) => {
          if (headerToken != nodeNameCol)
            nodeToAnnotationToValue(sample)(headerToken) = sp(index)
        }}
        loadedAnnotations += 1
      }}
      println("Loaded " + loadedAnnotations + " annotations from " + annotationFile)
    }}

    // walk the node tree
    def recurseNode(node: Node): D3Node = {
      val ret = new D3Node(node.getID, node.getLength, nodeToAnnotationToValue)
      // translate a list into a pair of nodes
      //println(node.getChildren.size)
      ret.children = node.getChildren.map { case(child) =>
        recurseNode(child)
      }.toList
      return ret
    }

    // we store all of the nodes as we make them
    var nodes = recurseNode(treeParser.getRoot)

    val outputFl = new PrintWriter(outputFile)
    outputFl.write(nodes.toJSONString(0))
    outputFl.close()

  }


}


// out container for d3 nodes
class D3Node(nm: String,
             ln: Double,
             masterAnnotations: HashMap[String,HashMap[String,String]]) {

  val name = nm
  val len = ln
  var annotations = new HashMap[String,String]
  var children = List[D3Node]()

  if (masterAnnotations contains name)
    annotations = masterAnnotations(name)

  var eventArray = Array[Int]()
  if (annotations contains "eventString") {
    val spl = annotations("eventString").split("-")

    val events = new ArrayBuffer[Event]()

    spl.zipWithIndex.foreach{case(tg,index) => {
        tg.split("\\&").foreach{
          subevt => {
            events += Event.toEvent(subevt,index)
          }
        }
      }}
    val hmid = HMID(events.toArray)

    eventArray = hmid.eventToPerBase(120, 380)
  }


  // why is finding a decent JSON serialization library so hard?
  // ---------------------------------------------------------------------------------
  def toJSONString(tabStop: Int): String = {
    val tabs = (0 until tabStop).map{_ => "\t"}.mkString("")
    var ret = tabs + "{\n" + tabs + "\"name\" : \"" + name + "\","
    ret += tabs + "\n" + tabs + "\"branchlength\" : \"" + len + "\""

    if (annotations.size > 0)
      ret += ",\n" + annotations.map{case(key,value) => tabs + "\t\"" + key + "\" : \"" + value + "\""}.mkString(",\n")

    if (eventArray.size > 0)
      ret += ",\n" + tabs + "\"events\": [" + eventArray.mkString(",") + "]"

    if (children.size > 0) {
      ret += ",\n" + tabs + "\"children\" : [\n"
      ret += children.map{case(child) => child.toJSONString(Math.min(10,tabStop + 1))}.mkString(",\n")
      ret += "\n" + tabs + "]"
    }


    ret += "}"
    ret
  }
*/
}