package main.scala
import scala.io._
import java.io._
import scala.math._


/**
  * Given a set of tree's output from MIX, find the best one (by proportion field of the output)
  */
case class BestTree(treeFile: File) {

  //get all the trees as a single line, MIX uses ';' to seperate them
  val one_line = Source.fromFile(treeFile).getLines().mkString("")

  // find the best tree
  var maxIndex = 0
  var bestScore = 0.0
  val trees = one_line.split(";").zipWithIndex.map{case(tree,index) => {
    val treeString = tree.split("\\[")
    val score = treeString(1).stripSuffix("]").toDouble
    if (score >= bestScore) {
      bestScore = score
      maxIndex = index
    }
    (treeString(1).stripSuffix("]").toDouble,treeString(0))
  }}.toArray

  val bestTreeString = trees(maxIndex)._2
}
