import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._


// pull out reads that match a certain UMI sequence from paired-end fastq files
// read in compressed input streams with scala source commands
def gis(s: String) = new GZIPInputStream(new BufferedInputStream(new FileInputStream(s)))

val umi = args(0)

val barcode1file = Source.fromInputStream(gis(args(1))).getLines().grouped(4)
val barcode2file = Source.fromInputStream(gis(args(2))).getLines().grouped(4)

val output1 = new PrintWriter(args(3))
val output2 = new PrintWriter(args(4))

barcode1file.zip(barcode2file).foreach{case(read1,read2) => {
  if (read1(1) startsWith args(0)) {
    output1.write(read1(0) + "\n" + read1(1) + "\n" + read1(2) + "\n" + read1(3) + "\n")
    output2.write(read2(0) + "\n" + read2(1) + "\n" + read2(2) + "\n" + read2(3) + "\n")
  }
}}
output1.close()
output2.close()
