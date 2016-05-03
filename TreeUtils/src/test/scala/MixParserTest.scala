package test.scala

import main.scala.{Edge, MixParser}
import org.scalatest.{Matchers, FlatSpec}

import scala.collection.mutable

/**
  * Created by aaronmck on 4/29/16.
  */
class MixParserTest extends FlatSpec with Matchers {

  val mixOutput17 = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2016_04_29_Tree_Data/data/parsimony_fish_17/fish17.parsimony.output"
  val eventToNumber17 = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2016_04_29_Tree_Data/data/parsimony_fish_17/17.eventToNumber"

  val treeToUse = 99 // last tree
  val parser = new MixParser(mixOutput17, eventToNumber17, treeToUse)
  val treeGenotypes = parser.activeTree.getOrElse(Array[Edge]())

  "MixParser" should "read the correct number of trees out" in {
    parser.currentTreeNumber should be (100)
  }

  "MixParser" should "get the correct number of nodes out" in {

    println(treeGenotypes(0).toFancyString)
    treeGenotypes.size should be (1001)
  }

  "MixParser" should "find the first node correctly" in {
    treeGenotypes(0).from should be ("root")

  }
  "MixParser" should "find the last node correctly" in {
    treeGenotypes(treeGenotypes.size - 1).from should be ("21")

  }

}
