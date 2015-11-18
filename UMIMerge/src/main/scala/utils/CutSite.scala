package utils

import java.io.File

import scala.io.Source

/**
 * Created by aaronmck on 11/9/15.
 */
case class CutSite(name: String, startOfSite: Int, cutSite: Int) {

}

object CutSite {
  def toCutSiteArray(input: File): Array[CutSite] = {
    Source.fromFile(input).getLines().drop(1).map{lt => CutSite(lt.split("\t")(0),lt.split("\t")(1).toInt,lt.split("\t")(2).toInt)}.toArray
  }
}