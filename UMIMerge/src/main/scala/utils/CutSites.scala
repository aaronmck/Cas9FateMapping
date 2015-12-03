package utils

import java.io.File

import scala.collection.mutable
import scala.io.Source

/**
 * Created by aaronmck on 11/17/15.
 */
class CutSites {
  val cutSites = new mutable.HashMap[Int, Int]()
  val startSites = new mutable.HashMap[Int, Int]()
  var fullSites = Array[Tuple3[String, Int, Int]]()
  var windows = Array[Tuple3[Int, Int, Int]]()
  var size = windows.size
}

object CutSites {

  def fromFile(cutsiteFile: File, windowSize: Int): CutSites = {

    val cut = new CutSites()

    Source.fromFile(cutsiteFile).getLines().drop(1).zipWithIndex.foreach { case (line, index) => {
      val sp = line.split("\t")

      cut.cutSites(index) = sp(2).toInt
      cut.startSites(index) = sp(1).toInt
      cut.fullSites :+=(sp(0), sp(1).toInt, sp(2).toInt)
      cut.windows :+=(sp(2).toInt - windowSize, sp(2).toInt, sp(2).toInt + windowSize)
      //println((sp(2).toInt - windowSize) + "\t" +  sp(2).toInt + "\t" + (sp(2).toInt + windowSize))
    }
    }

    cut.size = cut.windows.size
    return (cut)
  }

  def fromIntervals(triples: Array[Tuple3[Int, Int, Int]]): CutSites = {
    val cut = new CutSites()
    triples.zipWithIndex.foreach { case (trip, idx) => {
      cut.cutSites(idx) = trip._2
      cut.startSites(idx) = trip._2 // not right but ok for testing
      cut.fullSites :+=("UKNOWN", trip._1, trip._3)
      cut.windows :+=(trip._1, trip._2, trip._3)
    }
    }
    cut.size = triples.size
    return (cut)
  }
}
