package main.scala

import scala.collection.mutable

/**
  * Created by aaronmck on 5/2/16.
  */
class ParsimonyCollapser(rootNode: RichNode, parser: MixParser) {

  var bestCandidate = mutable.HashMap[String,Int]()
  var bestCandidateNode = mutable.HashMap[String,RichNode]()

  ParsimonyCollapser.checkCollapseNodes(rootNode,parser).foreach{case(tpl) =>
    if ((bestCandidate contains tpl._1) && (bestCandidate(tpl._1) < tpl._2)) {
      bestCandidate(tpl._1) = tpl._2
      bestCandidateNode(tpl._1) = tpl._3
    } else if (!(bestCandidate contains tpl._1)) {
      bestCandidate(tpl._1) = tpl._2
      bestCandidateNode(tpl._1) = tpl._3
    }
  }

  println("Best candidates list is:")
  var sum = 0
  bestCandidate.foreach{case(event,count) => {
    println(event + "\t" + count)
    sum += count
  }}
  println("sum = " + sum + " total " + ParsimonyCollapser.getNodeCount(rootNode))
}

object ParsimonyCollapser {
  /**
    *
    * @param node the rich node
    * @param parser the parsimony results
    */
  def checkCollapseNodes(node: RichNode, parser: MixParser): List[Tuple3[String,Int, RichNode]] = {
    (node.isCollapsible() ++ node.children.flatMap{case(nd) => checkCollapseNodes(nd,parser)}).toList
  }

  /**
    *
    * @param node the rich node
    */
  def getNodeCount(node: RichNode): Int = {
    1 + node.children.map{case(nd) => getNodeCount(nd)}.sum
  }
}
