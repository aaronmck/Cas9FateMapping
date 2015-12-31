package main.scala.utils

import java.io.File

/**
 * Created by aaronmck on 12/6/15.
 */
class UnmergedReadParser(readFiles: File) extends Iterator[ReadPair] {

  val parser = new ReadPairParser(readFiles)

  // the two reads we're returning
  var read1 : Option[RefReadPair] = None
  var read2 : Option[RefReadPair] = None

  nextReads()

  // fill the reads from our iterator
  def nextReads() {
    if (parser.hasNext)
      read1 = Some(parser.next())
    if (parser.hasNext)
      read2 = Some(parser.next())

    while (read1.isDefined && read2.isDefined && read1.get.read.name.split(" ")(0) != read2.get.read.name.split(" ")(0)) {
      read1 = read2
      if (parser.hasNext)
        read2 = Some(parser.next())
      else
        read2 = None
    }
  }

  override def hasNext: Boolean = read1.isDefined && read2.isDefined

  override def next(): ReadPair = {
    val ret = ReadPair(read1.get,read2.get)
    nextReads()
    return ret
  }


}

// the pair of reference and read from the fasta file
case class ReadPair(pair1: RefReadPair, pair2: RefReadPair)