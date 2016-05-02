package main.scala

import scala.collection.mutable


/**
  * the genotype
  *
  * @param name the name of the node
  * @param numberToEvent
  * @param eventToSites
  * @param edge
  * @param parentGenomes
  * @param numberofTargets
  */
case class GenotypeIdentity(name: String,
                            numberToEvent: mutable.HashMap[Int,String],
                            eventToSites: mutable.HashMap[String,Array[Int]],
                            edge: Edge,
                            parentGenomes: Array[String],
                            numberofTargets: Int = 10) {

  var genotypes = parentGenomes.clone()

  edge.chars.zipWithIndex.foreach{case(ch,index) => {
    ch match {
      case '1' => {
        val event = numberToEvent(index + 1)
        val positions = eventToSites(event)
        positions.foreach{pos => {
          if (genotypes(pos) != "NONE")
            throw new IllegalStateException("Unable to reassign genotype for edge " + edge.to + " and from " + edge.from)
          genotypes(pos)  = event
        }}
      }
      case _ => {}
    }
  }}

}