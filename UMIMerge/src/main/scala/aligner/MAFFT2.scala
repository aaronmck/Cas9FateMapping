package aligner

import java.io.{PrintWriter, File}

import main.scala.Consensus
import main.scala.utils.Utils
import utils.CutSites

import scala.collection.mutable.ArrayBuffer
import scala.io.Source
import scala.main._
import scala.sys.process._

/**
 * Created by aaronmck on 11/17/15.
 */
case class MAFFT2(forwardRead: SequencingRead, reverseRead: SequencingRead, refString: String, cuts: CutSites, primers: List[String]) {

  // ----------------------------------------------------------------------
  // setup a MAFFT run and farm it out the machine
  // ----------------------------------------------------------------------
  val readCombined = Consensus.consensus(MAFFT2.alignTo(Array[SequencingRead](forwardRead,reverseRead),Some(refString)))
  val readCombined2 = Consensus.consensus(MAFFT2.alignTo(Array[SequencingRead](readCombined),Some(refString)))

  println(readCombined.bases)
  println(readCombined2.bases)

}

/**
 * a simple case class to hold alignments -- results we get back from parsing reads aligned with MAFFTv7
 * @param refPos the reference position for the start of the event
 * @param refBase the reference bases over the event
 * @param readBase the read bases over the event
 * @param cigarCharacter the cigar character for the event -- I, D, and M are valid
 */
case class Alignment(val refPos: Int, refBase: String, readBase: String, cigarCharacter: String) {
  def combine(next: Alignment): Array[Alignment] =
    if (next.cigarCharacter != cigarCharacter)
      return Array[Alignment](this,next)
    else
      return Array[Alignment](Alignment(this.refPos,this.refBase+next.refBase,this.readBase+next.readBase,cigarCharacter))
  
  def prettyPrint: String = refPos + ":" + refBase + ":" + readBase + ":" + cigarCharacter

  def toEditString: String = readBase.length + "" + cigarCharacter + "+" + refPos + {if (cigarCharacter == "I") {"+" + readBase} else {""}}
}

object MAFFT2 {

  /**
   * call edits over a matched reference and read string
   * @param reference the reference string
   * @param read the read string, THE SAME LENGTH as the reference,  i.e. out of an MSA program
   * @param minMatchOnEnd the minimum number of matches to end our calling, if we don't see an event like this we backtrack Is and Ds until we find an M of this size
   * @param debugInfo should we dump a ton of debug info
   * @return a list of alignments over the read/ref combo
   */
  def callEdits(reference: String, read: String, minMatchOnEnd: Int, debugInfo: Boolean = false): List[Alignment] = {
    var referencePos = 0
    var inRef = false

    var refToEvent = List[Alignment]()

    reference.zip(read).foreach{ case(refBase: Char,readBase: Char) =>
      if (debugInfo)
        print("BASES: " + refBase + "," + readBase + " ")
      (refBase,readBase) match {
        case('-',readB) if !inRef => {/* we might be in the situation where we haven't started the real alignment, take the offset */}
        case('-',readB) if  inRef => { // insertion

          if (refToEvent.isEmpty)   refToEvent :+= Alignment(referencePos,refBase.toString,readBase.toString,"I")
          else                      refToEvent =   refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos,refBase.toString,readBase.toString,"I"))

          if (debugInfo)
            println("1: " + refToEvent.map{st => st.prettyPrint}.mkString("<>") + " " + refToEvent.size)
        }
        case(refB,'-') if !inRef => {
          // deletion before read starts -- we haven't aligned yet
          referencePos += 1
        }
        case(refB,'-') => { // deletion
          inRef = true
          if (refToEvent.isEmpty)   refToEvent :+= Alignment(referencePos,refBase.toString,readBase.toString,"D")
          else                      refToEvent =   refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos,refBase.toString,readBase.toString,"D"))
          referencePos += 1
          if (debugInfo)
            println("2: " + refToEvent.map{st => st.prettyPrint}.mkString("<>") + " " + refToEvent.size)
        }
        case(refB,readB) => { // match / mismatch
          inRef = true
          if (refToEvent.isEmpty)   refToEvent :+= Alignment(referencePos,refBase.toString,readBase.toString,"M")
          else                      refToEvent =   refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos,refBase.toString,readBase.toString,"M"))
          referencePos += 1
          if (debugInfo)
            println("3: " + refToEvent.map{st => st.prettyPrint}.mkString("<>") + " " + refToEvent.size)
        }


      }
    }
    // get a bit aggressive here -- start from both ends -- strip off insertions and deletions until we hit a match or mismatch of at least 10 bases
    return filterEnds(refToEvent,minMatchOnEnd)
  }

  /**
   * filter the alignments, if we have poor matches on the ends that allow us to accept and indel, roll it back until we've matched enough bases
   * @param eventList the list of alignments over the read
   * @param minMatch the minimum number of match bases to 'anchor' the ends, otherwise strip the trash off
   * @return a filtered alignment set
   */
  def filterEnds(eventList: List[Alignment], minMatch: Int) : List[Alignment] = {

    // make this as clear as possible
    var firstIndex = -1
    for (i <- 0 until eventList.size)
      if (firstIndex < 0 && eventList(i).cigarCharacter == "M" && eventList(i).readBase.length >= minMatch)
        firstIndex = i

    var lastIndex = -1
    for (i <- (eventList.size -1).until(-1,-1))
      if (lastIndex < 0 && eventList(i).cigarCharacter == "M" && eventList(i).readBase.length >= minMatch)
        lastIndex = i

    //println(firstIndex+ " " + (lastIndex+1) + " " + "PRE: " + eventList.mkString("-") + " POST: " + eventList.slice(firstIndex,lastIndex+1).mkString("-"))
    return (eventList.slice(firstIndex,lastIndex+1))
  }


  /**
   * given a read and reference, align and call events at the cut-sites
   * @param reference ref string
   * @param fwdRead read string
   * @param revRead reverse read string
   * @param cutSites the cutsutes to consider
   * @param minMatchOnEnd the minimum number of matches on the ends to keep from peeling crappy indels off
   * @param debug should we dump a lot of debug info
   * @return the rate of matching for cigar "M" bases for both reads and the array of events over cutsites
   */
  def cutSiteEvents(reference: String, fwdRead: SequencingRead, revRead: SequencingRead, cutSites: CutSites, minMatchOnEnd: Int, debug: Boolean = false): Tuple3[Double,Double,Array[String]] = {
    val alignmentsF = MAFFT2.alignTo(Array[SequencingRead](fwdRead),Some(reference),false, debug)
    val alignmentsR = MAFFT2.alignTo(Array[SequencingRead](revRead),Some(reference),true,  debug)

    val events1 = MAFFT2.callEdits(alignmentsF(0).bases, alignmentsF(1).bases, minMatchOnEnd, debug)
    val events2 = MAFFT2.callEdits(alignmentsR(0).bases, alignmentsR(1).bases, minMatchOnEnd, debug)

    if (debug) {
      println("Events1 : " + events1.mkString(","))
      println("Events2 : " + events2.mkString(","))
    }

    val matchRate1 = percentMatch(alignmentsF(0).bases,alignmentsF(1).bases)
    val matchRate2 = percentMatch(alignmentsR(0).bases,alignmentsR(1).bases)

    return (matchRate1,matchRate2,combineTwoReadEdits(events1, events2, cutSites,debug))
  }

  /**
   * do two intervals overlap
   * @param pos1Start
   * @param pos1End
   * @param pos2Start
   * @param pos2End
   * @return
   */
  def overlap(pos1Start: Int, pos1End: Int, pos2Start: Int, pos2End: Int): Boolean = (pos1Start,pos1End, pos2Start, pos2End) match {
    case (pos1S,pos1E, pos2S, pos2E) if pos1S < pos2S && pos1E > pos2E => true
    case (pos1S,pos1E, pos2S, pos2E) if pos2S < pos1S && pos2E > pos1E => true
    case (pos1S,pos1E, pos2S, pos2E) if pos1S < pos2E && pos1E > pos2S => true
    case (pos1S,pos1E, pos2S, pos2E) if pos2S < pos1E && pos2E > pos1S=> true
    case _ => false
  }

  /**
   * combine the edits over two reads
   * @param edits1 the set of edits from read 1
   * @param edits2 the set of edits from read 2
   * @param cutSites the cutsites we consider
   * @param debug should we dump a lot of debugging info
   * @return an array of events over the target cut sites
   */
  def combineTwoReadEdits(edits1: List[Alignment],edits2: List[Alignment], cutSites: CutSites, debug: Boolean = false): Array[String] = {
    var ret = Array[String]()
    cutSites.windows.foreach{case(start,cut,end) => {
      var candidates = Array[Alignment]()

      edits1.foreach{edit =>
        if (edit.cigarCharacter != "M" && overlap(start,end,edit.refPos,edit.refPos+edit.refBase.length))
          candidates :+= edit
      }

      edits2.foreach{edit =>
        if (edit.cigarCharacter != "M" && overlap(start,end,edit.refPos,edit.refPos+edit.refBase.length))
          candidates :+= edit
      }

      if (debug)
        println("Site: " + start + "-" + end + ": " + candidates.mkString("\t") + "<<<")

      if (candidates.size == 0)
        ret :+= "NONE"
      else if (candidates.size == 1)
        ret :+= candidates(0).toEditString
      else if (candidates.size == 2) {
        // Do collision detection here -- do we have the same event?
        if (candidates(0).toEditString == candidates(1).toEditString)
          ret :+= candidates(0).toEditString
        else
          ret :+= candidates(0).toEditString + "&" + candidates(1).toEditString
      } else {
        throw new IllegalStateException("Unable to process sites with more than 2 edits!: " + candidates.map{cd => cd.toEditString}.mkString(":"))
      }
    }}

    return ret
  }

  /**
   * for non gap bases, what is our matching proportion?
   * @param ref the reference string
   * @param read the read string of the same length as the reference string
   * @return a proportion of bases that match
   */
  def percentMatch(ref: String, read: String): Double = {
    var bases = 0
    var matches = 0
    ref.zip(read).map{case(refBase,readBase) => {
      if (refBase != '-' && readBase != '-') {
        if (refBase == readBase)
          matches += 1
        bases += 1
      }
    }}
    matches.toDouble / bases.toDouble
  }

  /**
   * align two sequences
   * @param reads a list of reads, which we pull sequence out of
   * @param ref the reference sequence
   * @param reverseComplement should we reverse complement the reads?
   * @param debug dump a lot of info
   * @return a array of aligned sequencing reads, the reference will be the first
   */
  def alignTo(reads: Array[SequencingRead],ref: Option[String], reverseComplement: Boolean = false, debug: Boolean = false): Array[SequencingRead] ={
    val tmp = java.io.File.createTempFile("UMIMergerBWA", ".txt")
    val tmpWriter = new PrintWriter(tmp)
    var readDirections = Array[ReadDirection]()

    // write the reads / reference to the input file
    if (ref.isDefined)
      tmpWriter.write(">reference\n" + ref.get + "\n")
    reads.foreach { rd =>
      if (reverseComplement)
        tmpWriter.write(">" + rd.name + "\n" + Utils.reverseComplement(rd.bases.filter { bs => bs != '-' }.mkString("")) + "\n")
      else
        tmpWriter.write(">" + rd.name + "\n" + rd.bases.filter { bs => bs != '-' }.mkString("") + "\n")
    }
    tmpWriter.close()

    // make an array of SequenceReads to store the result in
    var ret = Array[SequencingRead]()
    /**
     * Run MAFFT and capture the output
     */
    val out = new StringBuilder
    val err = new StringBuilder
    val logger = ProcessLogger(
      (o: String) => out.append(o + "\n"),
      (e: String) => err.append(e + "\n"))

    if (debug)
      println("mafft --maxiterate 1000 --genafpair " + tmp)
    ("mafft --maxiterate 1000 --genafpair " + tmp) ! logger

    /**
     * now readback each read figure out the alignment
     */
    var readNames = Array[String]()
    var readStrings = Array[String]()
    var currentRead = ""

    out.toString().split("\n") foreach { line =>
      if (line startsWith ">") {
        if (currentRead != "")
          readStrings :+= currentRead
        currentRead = ""
        readNames :+= line.stripPrefix(">")
      } else {
        currentRead += line.toUpperCase
      }
    }
    if (currentRead != "")
      readStrings :+= currentRead

    if ((ref.isDefined && readStrings.length != reads.size + 1) || (!ref.isDefined && readStrings.length != reads.size)) {
      throw new IllegalStateException("DIDNT get all our reads back: " + readStrings.length + " instead of " + reads.size)
    }
    for (i <- 0 until readStrings.length) {
      if (debug)
        println(readStrings(i))
      ret :+= SequencingReadFromNameBases(readNames(i),readStrings(i))
    }

    tmp.delete()
    return ret
  }


  /**
   * create a sequencing read from a name a base string, just to save time in the alignment return
   * @param name the name of the read
   * @param bases the bases of the read
   * @return a sequencing read
   */
  def SequencingReadFromNameBases(name: String, bases: String): SequencingRead = {
    return SequencingRead(name,bases,"H"*bases.length,ConsensusRead,"UNKNOWN")
  }

}