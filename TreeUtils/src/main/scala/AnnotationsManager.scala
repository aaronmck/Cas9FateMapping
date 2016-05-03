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
  val seperator = "\t"
  val eventSeperator = "_"

  println("loading the annotations from " + annotations.getAbsolutePath)
  println("loading the sample to clade information from " + sampleToClade.getAbsolutePath)

  val mappingFile = Source.fromFile(annotations.getAbsolutePath).getLines()
  val mappingFilesHeader = mappingFile.next()
  println("mapping files header " + mappingFilesHeader)
  val header = mappingFilesHeader.split("\t").zipWithIndex.map{case(tk,index) => (tk,index)}.toMap

  val cladeFile = Source.fromFile(sampleToClade.getAbsolutePath).getLines()
  val cladeHeaderString = cladeFile.next()
  println("cladeHeaderString files header " + cladeHeaderString)
  val cladeHeader = cladeHeaderString.split("\t").zipWithIndex.map{case(tk,index) => {println(tk); (tk,index)}}.toMap
  cladeHeader.foreach{case(id,value) => println(id + "\t" + value)}

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


  val annotationMapping = new HashMap[String,AnnotationEntry]()
  val cladeMapping = new HashMap[String,CladeEntry]()

  // map the annotation header
  mappingFile.foreach{line => {
    val sp = line.split(seperator)
    // println(sp(header("taxa")) + " - " + sp(header("event")))
    annotationMapping(sp(header("taxa"))) = AnnotationEntry(sp(header("taxa")),
      sp(header("sample")),
      sp(header("count")).toInt,
      sp(header("proportion")).toFloat,
      sp(header("event")))
  }}

  // now get the mapping for sample to clade and color
  cladeFile.foreach{line => {
    val sp = line.split(seperator)
    cladeMapping(sp(header("sample"))) = CladeEntry(sp(cladeHeader("sample")),
      sp(cladeHeader("clade")),
      sp(cladeHeader("color")))
  }}
}

// some containers
case class AnnotationEntry(taxa: String, sample: String, count: Int, proportion: Double, event:String)
case class CladeEntry(sample: String, clade: String, color:String)
