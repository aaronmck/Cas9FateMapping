import scala.io._
import java.io._
import scala.collection.mutable._
import scala.util.matching.Regex

val statsFiles = args(0).split(",")
val output = new PrintWriter(args(1))
val output2 = new PrintWriter(args(2))

val minBasesMatchingPostAlignment = 0.90
val passString = "PASS"
var totalMaxTargets = 12

output2.write("sample\tusedReads\tuniqueCigars\ttotalReads\n")
output.write("readName\ttype\tmergedReadLen\tread1len\tread2len\tmatchRate1\talignedBases1\tmatchRate2\talignedBases2\thasForwardPrimer\thasReversePrimer\tkeep\t")

statsFiles.foreach{statsFileLine => {
  val statsFile = new File(statsFileLine)

  if (statsFile.exists) {
    val sampleSplit = statsFile.getName.split("\\.")
    val sample = sampleSplit.slice(0,sampleSplit.length -1).mkString(".")
    val openStatsFile = Source.fromFile(statsFile).getLines()
    val header = openStatsFile.next()
    val headerTokens = new HashMap[String,Int]()
    header.split("\t").zipWithIndex.map{case(tk,ind) => headerTokens(tk) = ind}

    var maxTarget = 0
    (1 until (100+1)).foreach{case(ind) => if (headerTokens contains ("target" + ind)) maxTarget = ind}
    if (maxTarget > totalMaxTargets)
      totalMaxTargets = maxTarget
  }
}}

println("Max targets in dataset: " + totalMaxTargets)
output.write((1 until (totalMaxTargets+1)).map{case(ind) => "target" + ind}.mkString("\t") + "\n")

statsFiles.foreach{statsFileLine => {
  val statsFile = new File(statsFileLine)

  if (statsFile.exists) {
    val sampleSplit = statsFile.getName.split("\\.")
    val sample = sampleSplit.slice(0,sampleSplit.length -1).mkString(".")
    val openStatsFile = Source.fromFile(statsFile).getLines()
    val header = openStatsFile.next()
    val headerTokens = new HashMap[String,Int]()
    header.split("\t").zipWithIndex.map{case(tk,ind) => headerTokens(tk) = ind}

    var maxTarget = 0
    (1 until (totalMaxTargets+1)).foreach{case(ind) => if (headerTokens contains ("target" + ind)) maxTarget = ind}

    val uniqueCigar = new HashMap[String,Boolean]()

    var usedLines = 0
    var totalLines = 0
    openStatsFile.foreach{line => {
      val sp = line.split("\t")
      val fwdMatch = sp(headerTokens("matchRate1")).toDouble
      val passFail = sp(headerTokens("keep"))

      if (fwdMatch > minBasesMatchingPostAlignment /*&& revMatch > minBasesMatchingPostAlignment*/ && passFail == passString && !(line contains "WT_") && !(line contains "UNKNOWN")) {
        output.write(sample + "\t" + sp.mkString("\t") + "\t")
        output.write((maxTarget until totalMaxTargets).map{tk => "NA"}.mkString("\t") + "\n")

        val mergedCigar = (1 until (maxTarget +1)).map{case(ind) => sp(headerTokens("target" + ind))}.mkString("-")
        uniqueCigar(mergedCigar) = true

        usedLines += 1
      }
      totalLines += 1
    }}
    println("processed " + statsFile + " with " + usedLines + " reads used of a total " + totalLines + " reads")
    output2.write(sample + "\t" + usedLines + "\t" + uniqueCigar.size + "\t" + totalLines + "\n")

  } else {
    throw new IllegalStateException("Unable to find file " + statsFileLine)
  }
}}
output2.close()
output.close()
