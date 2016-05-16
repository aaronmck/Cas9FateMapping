package main.scala

import java.io.File

import scala.collection.mutable.HashMap
import scala.io.Source

/**
  * load up all the annotation files, and make a set of mappings
  * for each node name to a set of annotations.  On request,
  * annotate each node with it's appropriate info
  */
class AnnotationsManager(annotations: File, sampleToClade: File, cladeIdentities: Option[File]) {
  // annotationMapping
  // taxa    sample  count   proportion      event
  // N0      17_Brain        3501    0.14872557349192864     1D+141_NONE_NONE_1D+220_10D+237_1I+273+T_10D+289_67D+322_67D+322_67D+322
  //
  // cladeMapping
  // sample  clade   color
  // 17_Blood        Blood   #FF0000
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
  val sampleTotals = new HashMap[String,Int]()

  // map the annotation header
  mappingFile.foreach{line => {
    val sp = line.split(seperator)
    sampleTotals(sp(header("sample"))) = sampleTotals.getOrElse(sp(header("sample")),0) + sp(header("count")).toInt
    annotationMapping(sp(header("taxa"))) = AnnotationEntry(sp(header("taxa")),
      sp(header("sample")),
      sp(header("count")).toInt,
      sp(header("proportion")).toFloat,
      sp(header("event")))
  }}

  // now get the mapping for sample to clade and color
  cladeFile.foreach{line => {
    val sp = line.split(seperator)
    println("cladeMapping " + cladeHeader("sample") + " + " + sp(cladeHeader("sample")))
    cladeMapping(sp(cladeHeader("sample"))) = CladeEntry(sp(cladeHeader("sample")),
      sp(cladeHeader("clade")),
      sp(cladeHeader("color")))
  }}

  var eventDefinitionsToColors : Option[HashMap[String,Array[String]]] = None

  // *******************
  // deal with the optional clade assignment color matching
  if (cladeIdentities.isDefined) {
    val inputEvtDefs = Source.fromFile(cladeIdentities.get).getLines()
    val evtDefHeader = inputEvtDefs.next()
    if (evtDefHeader != "clade_name\tclade_event\tclade_color")
      throw new IllegalArgumentException("Unable to find the correct header on the clade identities file: we require clade_name<tab>clade_event<tab>clade_color")

    // now process the events definitions into colors
    eventDefinitionsToColors = Some(new HashMap[String,Array[String]]())

    inputEvtDefs.foreach{line => {
      val sp = line.split("\t")
      (eventDefinitionsToColors.get)(sp(2)) = sp(1).split(eventSeperator)
    }}
  }

  /**
    * lookup the clade color for this event; if there isn't one return black our default
    *
    * @param node the event node
    * @return a color string
    */
  def setNodeColor(node: RichNode, parentNode: Option[RichNode]): Tuple2[String,String] = {
    if (!eventDefinitionsToColors.isDefined)
      return ("nodecolor","black")

    var assigned_colors = Array[String]()
    eventDefinitionsToColors.get.foreach{case(color,arrayOfEvents) => {
      val containCount = arrayOfEvents.map{case(chkEvt) => if (node.parsimonyEvents contains chkEvt) 1 else 0}.sum

      var parentContains = 0
      if (parentNode.isDefined)
        parentContains = arrayOfEvents.map{case(chkEvt) => if (parentNode.get.parsimonyEvents contains chkEvt) 1 else 0}.sum

      // the second part of this expression is to deal with Jamie's clade choices
      if (containCount == arrayOfEvents.size && (node.children.size > 0 || parentContains == arrayOfEvents.size)) {
        assigned_colors :+= color
      }
    }}

    // check for conflict, only assign if it's one color
    if (assigned_colors.size == 1)
      return ("nodecolor",assigned_colors(0))
    return ("nodecolor","black")
  }

}

// some containers
case class AnnotationEntry(taxa: String, sample: String, count: Int, proportion: Double, event:String)
case class CladeEntry(sample: String, clade: String, color:String)
