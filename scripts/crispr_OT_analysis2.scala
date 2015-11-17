import scala.io._
import java.io._
import scala.collection.mutable._
import scala.sys.process._
import java.util.zip._
import scala.util.matching.Regex

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
def reverseComplement(str: String) = str.map { t => compBase(t) }.reverse.mkString("")

val targetStart = 177
val targetEnd = 189

val editStart = 142
val editEnd = 153

// command line args:
// 1) control SAM file
// 2) edit SAM file
// 3) output file
val controlUmiToTarget = new HashMap[String,String]()
val controlUmiCount = new HashMap[String,Int]()
val controlUmiAgree = new HashMap[String,Int]()
val controlUmiDisagree = new HashMap[String,Int]()


val output = new PrintWriter(args(2))

// controls
Source.fromFile(args(0)).getLines.foreach { line => {
  if (!(line startsWith "@")) {
    val sp = line.split("\t")
    val seq = sp(9)
    val postUmi = seq indexOf "CTGTTAGCCATA"
    val pretarget = seq indexOf "AGCTTACATCGCT"

    if (postUmi > 0 & pretarget > 0) {
      val target = seq.slice(pretarget+13,pretarget+13+20)
      val umi = seq.slice(postUmi-12,postUmi)

      controlUmiCount(umi) = controlUmiCount.getOrElse(umi,0) + 1
      if (!(controlUmiToTarget contains umi)) {
        controlUmiToTarget(umi) = target
      } else {
        val oldTarget = controlUmiToTarget(umi)
        if (target == oldTarget)
          controlUmiAgree(umi) = controlUmiAgree.getOrElse(umi,0) + 1
        else
          controlUmiDisagree(umi) = controlUmiDisagree.getOrElse(umi,0) + 1
      }
    }
  }
}}

// cases
val caseUmiAgree = new HashMap[String,Int]()
val caseUmiDisagree = new HashMap[String,Int]()
val caseUmiMissing = new HashMap[String,Int]()
output.write("umi\ttarget\tcontrol\tpos.edit\tcontrol.ok\tmatching\tcontrol.count\tcontrol.agree\tcontrol.disagree\n")
Source.fromFile(args(1)).getLines.foreach { line => {
  if (!(line startsWith "@")) {
    val sp = line.split("\t")
    val seq = sp(9)
    val postUmi = seq indexOf "CTGTTAGCCATA"
    val pretarget = seq indexOf "AGCTTACATCGCT"

    if (postUmi > 0 & pretarget > 0) {
      val target = seq.slice(pretarget+13,pretarget+13+20)
      val umi = seq.slice(postUmi-12,postUmi)

      val exists = controlUmiToTarget contains umi
      val matching = (target == controlUmiToTarget.getOrElse(umi,"GARBAGE"))
      if (!exists)
        caseUmiMissing(umi) = caseUmiMissing.getOrElse(umi,0) + 1
      else if (matching)
        caseUmiAgree(umi) = caseUmiAgree.getOrElse(umi,0) + 1
      else
        caseUmiDisagree(umi) = caseUmiDisagree.getOrElse(umi,0) + 1

      val posEdit = if(!(target.endsWith("AGT"))) "EDIT" else "NONE"
      val controlOK = if(!(controlUmiToTarget.getOrElse(umi,"GARBAGE") endsWith ("AGT"))) "BAD" else "OK"
      output.write(umi + "\t" + target + "\t" + controlUmiToTarget.getOrElse(umi,"GARBAGE") + "\t" + posEdit + "\t" + controlOK + "\t" + matching + "\t" + controlUmiCount.getOrElse(umi,0) + "\t" + controlUmiAgree.getOrElse(umi,0) + "\t" + controlUmiDisagree.getOrElse(umi,0) + "\n")
    }
  }
}}
println(caseUmiAgree.size)
println(caseUmiDisagree.size)
println(caseUmiMissing.size)
output.close()
