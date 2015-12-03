import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._

// read in compressed input streams with scala source commands
def gis(s: String) = new GZIPInputStream(new BufferedInputStream(new FileInputStream(s)))

val barcode1file = Source.fromInputStream(gis(args(0))).getLines().grouped(4)
val barcode2file = Source.fromInputStream(gis(args(1))).getLines().grouped(4)

val barcodeCombinations = new HashMap[String,HashMap[String,Int]]()

var lines = 0
barcode1file.zip(barcode2file).foreach{case(bar1,bar2) => {
  lines += 1
  val bar1Hash = barcodeCombinations.getOrElse(bar1(1),new HashMap[String,Int]())
  bar1Hash(bar2(1)) = bar1Hash.getOrElse(bar2(1),0) + 1

  barcodeCombinations(bar1(1)) = bar1Hash
}}

println(lines)
val output = new PrintWriter(args(2))
output.write("barcode1\tbarcode2\tcount\n")
barcodeCombinations.foreach{case(barcode1,hmap) =>
  hmap.foreach{case(barcode2,count) =>
    output.write(barcode1 + "\t" + barcode2 + "\t" + count + "\n")
  }}

output.close()
