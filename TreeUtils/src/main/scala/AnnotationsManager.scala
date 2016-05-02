package main.scala

import java.io.File

import scala.collection.mutable.HashMap
import scala.io.Source

/**
  * load up all the annotation files, and make a set of mappings
  * for each node name to a set of annotations.  On request,
  * annotate each node with it's appropriate info
  */
class AnnotationsManager(annotations: File, sampleToClade: File) {
  // taxa    sample  count   proportion      event
  // store leaf node names to their group and read

  val mappingFile = Source.fromFile(annotations.getAbsolutePath).getLines()
  val cladeFile = Source.fromFile(annotations.getAbsolutePath).getLines()
  val header = mappingFile.next().split("\t").zipWithIndex.map{case(tk,index) => (tk,index)}.toMap
  val cladeHeader = cladeFile.next().split("\t").zipWithIndex.map{case(tk,index) => (tk,index)}.toMap

  val annotationMapping = new HashMap[String,AnnotationEntry]()
  val cladeMapping = new HashMap[String,CladeEntry]()

  // check the annotation header
  if (!(header contains "taxa")) throw new IllegalStateException("Unable to find taxa in your header")
  if (!(header contains "sample")) throw new IllegalStateException("Unable to find sample in your header")
  if (!(header contains "count")) throw new IllegalStateException("Unable to find count in your header")
  if (!(header contains "proportion")) throw new IllegalStateException("Unable to find count in your header")
  if (!(header contains "event")) throw new IllegalStateException("Unable to find eventString in your header")

  // check the clade header
  if (!(cladeHeader contains "sample")) throw new IllegalStateException("Unable to find sample in your header")
  if (!(cladeHeader contains "clade")) throw new IllegalStateException("Unable to find clade in your header")
  if (!(cladeHeader contains "color")) throw new IllegalStateException("Unable to find color in your header")

  // map the annotation header
  mappingFile.foreach{line => {
    val sp = line.split("\t")
    annotationMapping(sp(header("sample"))) = AnnotationEntry(sp(header("taxa")),
      sp(header("sample")),
      sp(header("count")).toInt,
      sp(header("proportion")).toFloat,
      sp(header("event")))
  }}

  // now get the mapping for sample to clade and color
  cladeFile.foreach{line => {
    val sp = line.split("\t")
    cladeMapping(sp(header("sample"))) = CladeEntry(sp(header("sample")),
      sp(header("clade")),
      sp(header("color")))
  }}
}

// some containers
case class AnnotationEntry(taxa: String, sample: String, count: Int, proportion: Double, event:String)
case class CladeEntry(sample: String, clade: String, color:String)
