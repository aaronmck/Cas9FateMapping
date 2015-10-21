import aaronmck.github.io.seq.dna.alignment._
import scala.io._
import java.io._
import java.util.zip.GZIPInputStream
import scala.math.sqrt
import scala.util.matching.Regex
import scala.collection._
import aaronmck.github.io.seq.util.Util
import aaronmck.github.io.seq.util.Cigar

/**
  * Calculate the concensous bases per position, and compare this to the reference bases
  * 
  * args:
  * 1: the reference file
  * 2: the input reads, compressed as a fastq
  * 3: the output file
  * 
  */


def baseToInt(base: Char): Int = base match {
  case 'A' => 0
  case 'C' => 1
  case 'G' => 2
  case 'T' => 3
  case _ => 4
}
def intToBase(pos: Int): Char = pos match {
  case 0 => 'A'
  case 1 => 'C'
  case 2 => 'G'
  case 3 => 'T'
  case _ => 'N'
}

val reference = Source.fromFile(args(0)).getLines().drop(1).map{ln => ln}.mkString("").toUpperCase


val mp = mutable.HashMap[Int,Array[Int]]()


val inp = args(1)
Source.fromInputStream(Util.gis(inp)).getLines().grouped(4).foreach{read => {
  read(1).toUpperCase.zipWithIndex.foreach{case(base,index) => {
    var cur = mp.getOrElse(index,Array[Int](0,0,0,0,0))
    cur(baseToInt(base)) += 1
    mp(index) = cur
  }}
}}

val output = new PrintWriter(args(2))
output.write("pos\tref\trefProp\tA\tC\tG\tT\tN\tpropHighest\tconcBase\n")
for (i <- 0 until 350) {
  val hits = mp.getOrElse(i,Array[Int](0,0,0,0,0))
  val propHighest = hits.max.toDouble / hits.sum.toDouble
  val maxBase = intToBase(hits.indexOf(hits.max))
  val refBase = if (i < reference.length) reference(i) else 'N'
  val propRef = hits(baseToInt(refBase)).toDouble / hits.sum.toDouble
  output.write(i + "\t" + refBase + "\t" + propRef + "\t" + hits.mkString("\t") + "\t" + propHighest + "\t" + maxBase + "\n")
}
output.close()
