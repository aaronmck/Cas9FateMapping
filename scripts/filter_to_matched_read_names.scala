import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._

// read in compressed input streams with scala source commands
def gis(s: String) = new GZIPInputStream(new BufferedInputStream(new FileInputStream(s)))

// forward and reverse reads
val forwardReads = Source.fromInputStream(Utils.gis(args(0))).getLines().grouped(4)
val reverseReads = Source.fromInputStream(Utils.gis(args(1))).getLines().grouped(4)

val interleavedFile = new PrintWriter(args(2))

fowardReads.zip(reverseReads).foreach{case(read1,read2) =>
  interleavedFile.write(read1(0) + "\n")
  interleavedFile.write(read1(1) + "\n")
  interleavedFile.write(read1(2) + "\n")
  interleavedFile.write(read1(3) + "\n")

  interleavedFile.write(read2(0) + "\n")
  interleavedFile.write(read2(1) + "\n")
  interleavedFile.write(read2(2) + "\n")
  interleavedFile.write(read2(3) + "\n")
}

interleavedFile.close()
