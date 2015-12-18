import scala.io._
import java.io._
import scala.collection.mutable._
import scala.util.matching.Regex

val statsFiles = args(0).split(",")
val output = new PrintWriter(args(1))
val output2 = new PrintWriter(args(2))

output2.write("sample\tusedUMI\ttotalUMI\n")

val minBasesMatchingPostAlignment = 0.90
val passString = "PASS"
//----------------------------------------------------------------------------------------------------------------
// now for each file in the tearsheet do the following: load the barcodes from the sam file,
// load and add those to the appropriate row of the stats file
//----------------------------------------------------------------------------------------------------------------
var outputHeader = false

statsFiles.foreach{statFileLine => {
  /*val sp = line.split("\t")
  val outputDir = sp(3)
  val sample = sp(0)

  val statsFile = new File(outputDir + "/" + sample + "/" + sample + ".stats")*/
  val statsFile = new File(statsFileLine)

  if (statsFile.exists) {
    val openStatsFile = Source.fromFile(statsFile).getLines()
    val header = openStatsFile.next()
    val headerTokens = new HashMap[String,Int]()
    header.split("\t").zipWithIndex.map{case(tk,ind) => headerTokens(tk) = ind}

    if (!outputHeader) {
      output.write("sample\t" + header + "\n")
      outputHeader = true
    }

    var usedLines = 0
    var totalLines = 0
    openStatsFile.foreach{line => {
      val sp = line.split("\t")
      val fwdMatch = sp(headerTokens("fwdBasesMatching")).toDouble
      val revMatch = sp(headerTokens("revBasesMatching")).toDouble
      val passFail = sp(headerTokens("fail.reason"))

      if (fwdMatch > minBasesMatchingPostAlignment /*&& revMatch > minBasesMatchingPostAlignment*/ && passFail == passString && !(line contains "&")) {
        output.write(sample + "\t" + line + "\n")
        usedLines += 1
      }
      totalLines += 1
    }}
    println("processed " + statsFile + " with " + usedLines + " UMIs used of a total " + totalLines + " UMIs")
    output2.write(sample + "\t" + usedLines + "\t" + totalLines + "\n")

  } else {
    throw new IllegalStateException("Unable to find file " + statsFileLine)
  }
}}
output2.close()
output.close()
