package utils

import java.io.File

import scala.collection.mutable
import scala.io.Source

/**
 * Stores all the information about the cut sites in this target HMID cassette
 */
class CutSites {
  // mapping of the cutsite index to the base position of the cutsite
  val cutSites = new mutable.HashMap[Int, Int]()

  // mapping cutsite index to the starting base position of the target
  val startSites = new mutable.HashMap[Int, Int]()

  // the full target sites -- from the start of the target to the start of the PAM
  var fullSites = Array[Tuple3[String, Int, Int]]()

  // the window over the cutsite to consider -- if an indel overlaps this regions we count it as an edit
  var windows = Array[Tuple3[Int, Int, Int]]()

  // our total size of targets in the cassette
  var size = windows.size
}

object CutSites {
  // a constant for the distance between the cutsite and the end of the pam, used to define the full target window
  val cutsiteToPamDistance = 6

  /**
   * load a cutsite object from a CSV file on disk
   * @param cutsiteFile the input file
   * @param windowSize how large of a window to include
   * @return the cutsite object
   */
  def fromFile(cutsiteFile: File, windowSize: Int): CutSites = {

    val cut = new CutSites()

    Source.fromFile(cutsiteFile).getLines().drop(1).zipWithIndex.foreach { case (line, index) => {
      val sp = line.split("\t")

      val name = sp(0)
      val start = sp(1).toInt - 1
      val cutSite = sp(2).toInt - 1

      cut.cutSites(index) = cutSite
      cut.startSites(index) = start
      cut.fullSites :+=(name, start, cutSite + cutsiteToPamDistance)
      cut.windows :+=(cutSite - windowSize, cutSite, cutSite + windowSize)
    }
    }

    cut.size = cut.windows.size
    return (cut)
  }

  /**
   * create a cutsite object from a set of intervals, mainly used in testing
   * @param triples the cut site windows to use as a downstream-cut-upstream window
   * @return a cutsite object for the triples
   */
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
