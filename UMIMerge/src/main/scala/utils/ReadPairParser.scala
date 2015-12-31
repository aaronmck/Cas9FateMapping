package main.scala.utils

import java.io.File

import scala.collection.mutable.ArrayBuilder
import scala.io.Source
import scala.main.SequencingRead

/**
 * Created by aaronmck on 12/6/15.
 */
class ReadPairParser(readFile: File) extends Iterator[RefReadPair] {

  val reads = new PeekFilelineIterator(Source.fromFile(readFile).getLines)

  var returnReads: Option[RefReadPair] = getNextReadPair()

  def getNextReadPair(): Option[RefReadPair] = {

    // store our reference and reads
    var refName = ""
    var readName = ""
    var refString = ArrayBuilder.make[String]
    var readString = ArrayBuilder.make[String]
    var inRead = false
    var started = false
    // the read file

    while(reads.hasNext) {
      var line = reads.peek
      if ((line startsWith ">") && inRead) {
        return Some(RefReadPair(SequencingRead.readFromNameAndSeq(refName, refString.result().mkString("")),
          SequencingRead.readFromNameAndSeq(readName, readString.result().mkString(""))))
      }

      else if (line startsWith ">") {
        if (refName == "") {
          refName = line.stripPrefix(">")
        } else {
          inRead = true
          readName = line.stripPrefix(">")
        }
      } else if (inRead) {
        readString += line
      } else {
        refString += line
      }

      reads.next()
    }
    return None
  }

  override def hasNext: Boolean = returnReads isDefined

  override def next(): RefReadPair = {
    val ret = returnReads.get
    returnReads = getNextReadPair()
    return ret
  }
}


class PeekFilelineIterator(iter: Iterator[String]) extends Iterator[String] {

  var nextString: Option[String] = if (iter.hasNext) Some(iter.next()) else None

  override def hasNext: Boolean = nextString isDefined

  override def next(): String = {
    val ret = nextString.get
    nextString = if (iter.hasNext) Some(iter.next()) else None
    return ret
  }

  def peek: String = nextString.get
}

case class RefReadPair(reference: SequencingRead, read: SequencingRead)
