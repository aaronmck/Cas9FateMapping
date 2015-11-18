package aligner

import scala.util.matching.Regex

/**
 * Created by aaronmck on 11/17/15.
 */
object Cigar {
  val cigarRegex = new Regex("(\\d+)([MIDSNHP=X])")

  // convert a phred scaled qual score to an int
  def qual33ToInt(qual: Char): Int = qual.toInt - 33

  // zip by the cigar string, returning a tuple of:
  // 1: read base
  // 2: reference base
  // 3: do they match? no if an indel
  // 4: position in the reference
  // 5: cigar character
  def zipCigar(baseString: String,referenceString: String,cigarString: String, quality: String) : Array[Tuple6[Char,Char,Int,Boolean,Int,Char]] = {
    var returnTuples = Array[Tuple6[Char,Char,Int,Boolean,Int,Char]]()
    var baseIndex = 0
    var refIndex = 0
    var qualPos = 0
    (cigarRegex findAllIn cigarString).matchData.foreach{matching => {
      (0 until matching.subgroups(0).toInt).foreach(ignored => {
        matching.subgroups(1) match {
          case "M" => {
            try {
              returnTuples :+= (baseString(baseIndex),referenceString(refIndex),qual33ToInt(quality(qualPos)),baseString(baseIndex) == referenceString(refIndex),refIndex,'M')
            } catch {
              case e: Exception => { println(cigarString); throw e}
            }
            baseIndex += 1
            qualPos += 1
            refIndex += 1
          }
          case "I" => {
            returnTuples :+= (baseString(baseIndex),'-',qual33ToInt(quality(qualPos)),false,refIndex,'I')
            baseIndex += 1
            qualPos += 1
          }
          case "D" => {
            returnTuples :+= ('-',referenceString(refIndex),0,false,refIndex,'D')
            refIndex += 1
          }
          case "S" => {
            returnTuples :+= ('+','+',qual33ToInt(quality(qualPos)),false,refIndex,'S')
            baseIndex += 1
            qualPos += 1
          }
          case _ => {
            println("Uknown cigar")
          }
        }
      })
    }}
    returnTuples
  }

  // find the number of CIGAR events as a proportion of our read length
  def cigarLength(cigar: String) : Int =
    (cigarRegex findAllIn cigar).matchData.map{ matching => {
      matching.subgroups(1) match {
        case "M" => matching.subgroups(0).toInt
        case "I" => matching.subgroups(0).toInt
        case _ => 0
      }}}.sum


  // find the number of CIGAR events as a proportion of our read length
  def cigarComplexity(cigar: String, readStrLen: Int) : Double = (cigarRegex findAllIn cigar).matchData.length.toDouble / readStrLen.toDouble
}
