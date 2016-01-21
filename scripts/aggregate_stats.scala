import scala.io._
import java.io._
import scala.collection.mutable._
import scala.util.matching.Regex

// command line parameters
val statsFiles = args(0).split(",")
val output     = new PrintWriter(args(1))
val output2    = new PrintWriter(args(2))
var isUMI      = args(3) match {
  case("UMI") => true
  case("NOUMI") => false
  case(_) => throw new IllegalStateException("Unable to parse proper UMI arg from: " + args(3))
}


// columns to keep
val columnsToKeep =
  if (isUMI) Array[String]("UMI","keptPCT","fail.reason","initF","filteredF","finalF","initR","filteredR","finalR","fwdBasesMatching","revBasesMatching")
  else Array[String]("readName","type","mergedReadLen","read1len","read2len","matchRate1","alignedBases1","matchRate2","alignedBases2","hasForwardPrimer","hasReversePrimer","keep")

output2.write("sample\tusedEntries\ttotalEntries\n")
output.write("sample\t" + columnsToKeep.mkString("\t") + "\t")

val minBasesMatchingPostAlignment = 0.90
val passString = "PASS"
var outputHeader = false

var totalMaxTargets = 0

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
output.write((1 until (totalMaxTargets+1)).map{case(ind) => "target" + ind}.mkString("\t") + "\tHMID\tisMerged\n")


statsFiles.foreach{statsFileLine => {
  val statsFile = new File(statsFileLine)

  if (statsFile.exists) {
    val sample = statsFile.getName.split("\\.")(0)
    val openStatsFile = Source.fromFile(statsFile).getLines()
    val header = openStatsFile.next()
    val headerTokens = new HashMap[String,Int]()
    header.split("\t").zipWithIndex.map{case(tk,ind) => headerTokens(tk) = ind}
    val colMappingsToKeep = headerTokens.keys.filter(key => columnsToKeep contains key).map{key => (key,headerTokens(key))}.toMap

    if (colMappingsToKeep.size != columnsToKeep.size)
      throw new IllegalStateException("Unable to find all of the required header keys in the input file " + statsFile)

    val sortedColumns = colMappingsToKeep.values.toArray
    scala.util.Sorting.quickSort(sortedColumns)

    var usedLines = 0
    var totalLines = 0
    openStatsFile.foreach{line => {
      val sp = line.split("\t")
      val fwdMatch = if (isUMI) sp(headerTokens("fwdBasesMatching")).toDouble else sp(headerTokens("matchRate1")).toDouble
      val isMerged = if (isUMI) sp(headerTokens("readR")).toUpperCase == "MERGED" else sp(headerTokens("type")).toUpperCase == "MERGED"
      val isMergedStr = if (isMerged) "merged" else "tworead"
      val passFail = if (isUMI) sp(headerTokens("fail.reason")) else sp(headerTokens("keep"))

      if (fwdMatch > minBasesMatchingPostAlignment && passFail == passString && !(line contains "WT_") && !(line contains "UNKNOWN")) {
        val targets = (1 until (totalMaxTargets +1)).map{case(ind) => if (headerTokens contains ("target" + ind)) sp(headerTokens("target" + ind)) else "NA"}
        val mergedCigar = targets.filter(tk => tk != "NA").mkString("-")
        val mergedTabCigar = targets.mkString("\t")
        
        output.write(sample + "\t" + sortedColumns.map{case(col) => sp(col)}.mkString("\t") + "\t" + mergedTabCigar + "\t" + mergedCigar + "\t" + isMergedStr + "\n")
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
