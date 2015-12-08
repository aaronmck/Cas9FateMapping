import scala.collection.immutable._
import scala.io._
import java.io._

case class CRISPRCircle(bufferSize: Int, outputFile: String, cRISPRs: Array[String], maxMismatches: Int = 4) {
  val circBuffer = new Array[Char](bufferSize)
  var added = 0
  val output = new PrintWriter(outputFile)
  var currentPos = 0
  var contig = "UNKNOWN"

  def addLine(line: String) {line.foreach{chr => {addBase(chr)}}}

  def addBase(chr: Char) {
    added +=1
    circBuffer(currentPos) = chr
    if (added >= bufferSize) {

      //val curStr = toTarget()
      //val revStr = reverseCompString(curStr)

      //cRISPRs.foreach{crispr => {
      //  checkCRISPR(crispr, curStr, revStr, added)
      //}}
      cRISPRs.foreach{crispr => {
        val fw = checkCRISPRFw(crispr)
        if (fw.isDefined)
          output.write(contig + "\t" + (added - bufferSize) + "\t" + (added) + "\t" + crispr + "\t" + fw.get._1 + "\t+\t" + fw.get._2 + "\n")
        val rv = checkCRISPRRev(crispr)
        if (rv.isDefined)
          output.write(contig + "\t" + (added - bufferSize) + "\t" + (added) + "\t" + crispr + "\t" + rv.get._1 + "\t-\t" + rv.get._2 + "\n")
      }}
    }
    currentPos = (currentPos + 1) % bufferSize

  }

  def reset(cntig: String) {
    contig = cntig
    println("loading " + cntig)
    currentPos = 0
    added = 0
  }

  def checkCRISPR(crispr: String, fwStr: String, bkStr: String, pos: Int, maxChecked: Int = 20) {
    val fwMatch = crispr.zip(fwStr.slice(0,maxChecked)).map{case(fr,st) => if(fr == st) 0 else 1}.sum
    val bkMatch = crispr.zip(bkStr.slice(0,maxChecked)).map{case(fr,st) => if(fr == st) 0 else 1}.sum
    if (fwMatch <= maxMismatches)
      output.write(contig + "\t" + (added - bufferSize) + "\t" + (added) + "\t" + crispr + "\t" + fwStr + "\t+\t" + fwMatch + "\n")
    if (bkMatch <= maxMismatches)
      output.write(contig + "\t" + (added - bufferSize) + "\t" + (added) + "\t" + crispr + "\t" + bkStr + "\t-\t" + bkMatch + "\n")
  }


  def checkCRISPRFw(crispr: String, maxChecked: Int = 20): Option[Tuple2[String,Int]] = {
    //if (circBuffer(currentPos) != 'G')
    //  return None

    var pos = currentPos
    // println("pos = " + pos)
    var cPos = 0
    var str = new StringBuilder()
    var mismatch = 0

    do {
      pos = pos + 1
      if (pos >=bufferSize) pos = 0
      if (circBuffer(pos) == 'N')
        return None
      //println("adding " + circBuffer(pos))
      str.append(circBuffer(pos))
      if (circBuffer(pos) != crispr(cPos))
        mismatch += 1
      if (mismatch > maxMismatches)
        return None
      cPos += 1
    } while(currentPos != pos && (cPos) < maxChecked)

    while (pos != currentPos) {
      pos = pos + 1
      if (pos >=bufferSize) pos = 0
      str.append(circBuffer(pos))
    }

    //println("hit at\n" + str + " for\n" + crispr + " mis = " + mismatch)
    Some(Tuple2[String,Int](str.toString,mismatch))
  }

  def checkCRISPRRev(crispr: String, maxChecked: Int = 20): Option[Tuple2[String,Int]] = {
    //if (((currentPos+1) == circBuffer.size && circBuffer(0) != 'C') || ((currentPos + 1) != circBuffer.size && circBuffer(currentPos+1) != 'C'))
    //  return None

    var pos = currentPos + 1
    // println("pos = " + pos)
    var cPos = 0
    var str = new StringBuilder()
    var mismatch = 0

    do {
      pos = pos - 1
      if (pos < 0) pos = circBuffer.size - 1
      if (circBuffer(pos) == 'N')
        return None
      str.append(reverseComp(circBuffer(pos)))
      if (reverseComp(circBuffer(pos)) != crispr(cPos))
        mismatch += 1
      if (mismatch > maxMismatches)
        return None
      cPos += 1
    } while(cPos < maxChecked)

    for (i <- 0 until 3) {
      pos = pos - 1
      if (pos < 0) pos = circBuffer.size - 1
      str.append(reverseComp(circBuffer(pos)))
    }

    //println("hit at\n" + str + " for\n" + crispr + " mis = " + mismatch)
    Some(Tuple2[String,Int](str.toString,mismatch))
  }


  def reverseComp(c: Char): Char = if (c == 'A') 'T' else if (c == 'C') 'G' else if (c == 'G') 'C' else if (c == 'T') 'A' else c
  def reverseCompString(str: String): String = str.map{reverseComp(_)}.reverse.mkString

  def toTarget(): String = circBuffer.slice(currentPos % bufferSize,bufferSize).mkString + circBuffer.slice(0,currentPos % bufferSize).mkString
  def close() = output.close()
}

var hits = Array[String]()


val molly_hits = Source.fromFile(args(0)).getLines()
molly_hits.foreach{ln => {
  val sp = ln.split("\t")
  hits :+= sp(3)
}}
val cls = CRISPRCircle(23,args(1),hits) // "TCTTAAGCAGAACAAGGGCAGGG")


var lineCount = 0
Source.fromFile(args(2)).getLines().foreach{line => {
  if (line.startsWith(">"))
    cls.reset(line.split(" ")(0).slice(1,line.split(" ")(0).length))
  else {
    cls.addLine(line.toUpperCase)
    lineCount += 1

    if (lineCount % 100000 == 0)
      println("basecount = " + cls.added)
  }
}}
cls.close()
