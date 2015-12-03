import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._


// pull out reads that match a certain UMI sequence from paired-end fastq files
// read in compressed input streams with scala source commands
def gis(s: String) = new GZIPInputStream(new BufferedInputStream(new FileInputStream(s)))


val lines = Source.fromInputStream(gis(args(0))).getLines().grouped(4)
val f1 = "ATCTCGAGCTCAA"
val f2 = "GAATTCTCGACCT"

lines.foreach{reads => {
  val id1 = (reads(1) indexOf f1)
  val id2 = (reads(1) indexOf f2)
  if (id1 >=0 && id2 >= 0)
    println(id1 + " " + id2 + "\t" + reads(1).slice(id2,id2+60))
}}
