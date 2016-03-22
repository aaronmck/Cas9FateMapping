import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._

// read in compressed input streams with scala source commands
def gis(s: String) = new GZIPInputStream(new BufferedInputStream(new FileInputStream(s)))

// complement a base as a character
  def compBase(b: Char): Char = b match {
    case 'A' => 'T'
    case 'C' => 'G'
    case 'G' => 'C'
    case 'T' => 'A'
    case 'a' => 't'
    case 'c' => 'g'
    case 'g' => 'c'
    case 't' => 'a'
    case _ => 'N'
  }

// reverse complement a string of DNA bases
def reverseComplement(str: String) = str.map{t => compBase(t)}.reverse.mkString("")

// forward and reverse reads are compressed
if (args(0) endsWith ".gz") {
  val forwardReads = Source.fromInputStream(gis(args(0))).getLines().grouped(4)
  val reverseReads = Source.fromInputStream(gis(args(1))).getLines().grouped(4)

  val interleavedFile = new PrintWriter(args(2))

  forwardReads.zip(reverseReads).foreach{case(read1,read2) =>
    interleavedFile.write(read1(0) + "\n")
    interleavedFile.write(read1(1) + "\n")
    interleavedFile.write(read1(2) + "\n")
    interleavedFile.write(read1(3) + "\n")

    interleavedFile.write(read2(0) + "\n")
    interleavedFile.write(reverseComplement(read2(1)) + "\n")
    interleavedFile.write(read2(2) + "\n")
    interleavedFile.write(read2(3) + "\n")
  }

  interleavedFile.close()
} else {
  val forwardReads = Source.fromFile(args(0)).getLines().grouped(4)
  val reverseReads = Source.fromFile(args(1)).getLines().grouped(4)

  val interleavedFile = new PrintWriter(args(2))

  forwardReads.zip(reverseReads).foreach{case(read1,read2) =>
    interleavedFile.write(read1(0) + "\n")
    interleavedFile.write(read1(1) + "\n")
    interleavedFile.write(read1(2) + "\n")
    interleavedFile.write(read1(3) + "\n")

    interleavedFile.write(read2(0) + "\n")
    interleavedFile.write(reverseComplement(read2(1)) + "\n")
    interleavedFile.write(read2(2) + "\n")
    interleavedFile.write(read2(3) + "\n")
  }

  interleavedFile.close()
}
