package utils

import java.io.File

import scala.collection.mutable
import scala.io.Source

/**
 * Created by aaronmck on 11/17/15.
 */
case class CutSites(cutsiteFile: File, windowSize: Int) {
  val cutSites = new mutable.HashMap[Int,Int]()
  val startSites = new mutable.HashMap[Int,Int]()
  var fullSites = Array[Tuple3[String,Int,Int]]()
  var windows = Array[Tuple3[Int,Int,Int]]()

  Source.fromFile(cutsiteFile).getLines().drop(1).zipWithIndex.foreach{case(line,index) => {
    val sp = line.split("\t")

    cutSites(index) = sp(2).toInt
    startSites(index) = sp(1).toInt
    fullSites :+= (sp(0),sp(1).toInt,sp(2).toInt)
    windows :+= (sp(2).toInt - windowSize,sp(2).toInt,sp(2).toInt + windowSize)
  }}

  var size = windows.size
}
