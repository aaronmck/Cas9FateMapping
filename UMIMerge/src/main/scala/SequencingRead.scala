package scala.main

import main.scala.utils.Utils

/**
 * holds a sequencing read
 */
case class SequencingRead(name: String, bases: String, quals: String, readOrientation: ReadDirection, umi: String, cigar: Option[String] = None, position: Int = 0) {
  if (bases.length != quals.length && quals != "*")
    throw new IllegalArgumentException("Read " + name + " has an unequal base and qual string!\n" + bases + "\n" + quals + "\n")

  val intQuals = quals.map{qual => Utils.phredCharToQscore(qual)}.toArray
  def length = bases.length

  /**
   * find the last base in the string before a series of dashes
   * @return an integer of the last position in a read that's not a dash
   */
  def trueEnd(): Int = {
    ((bases.length - 1).until(-1,-1)).foreach{i => if (bases(i) != '-') return i}
    return -1
  }

  /**
   * slice a read at a position, producing a read with appropriate qual scores
   * @param fromPos the base to start at, inclusive
   * @param toPos the base to end at, inclusive
   * @return a new SequencingRead representing the sliced-down version
   */
  def slice(fromPos: Int, toPos: Int): SequencingRead = {
    if (fromPos < 0)
      throw new IllegalArgumentException("from position of " + fromPos + " < 0")
    if (toPos > bases.length)
      throw new IllegalArgumentException("to position of " + toPos + " > readlength of " + toPos)
    if (fromPos >= toPos)
      throw new IllegalArgumentException("from position of " + fromPos + " greater than or equal to toPos of " + toPos)

    SequencingRead(name, bases.slice(fromPos, toPos),quals.slice(fromPos, toPos), readOrientation, umi)
  }

  /**
   * filter the reads down by a specific combination of quality score drop over a window
   * @param windowSize the sliding window to filter over
   * @param minWindowQual the qual score we have to achive over a window to keep the rest of the read
   */
  def qualityThresholdRead(windowSize: Int = 5, minWindowQual: Double = 10): SequencingRead = {
    val cutPos = intQuals.toArray.sliding(windowSize).zipWithIndex.map{case(basesInWindow,index) => {
      if (basesInWindow.sum / windowSize.toDouble < minWindowQual) {
        //println("for read " + bases + " found a window with qual of " +  (basesInWindow.sum / windowSize.toDouble) + " at position " + index)
        index
      } else
        0
    }}.filter(x => x != 0).toArray
    //println("quals = " + (quals.zip(bases).zipWithIndex.map{case((ql,bs),index) => index + "--" + bs + "," + Utils.phredCharToQscore(ql)}.mkString(";")))
    //println("for read " + bases + " the cut pos is " + cutPos.mkString(","))
    if (cutPos.size > 0) {
      //println(cutPos(0))
      val cutMinusWindow = math.max(0, cutPos(0) - minWindowQual)
      SequencingRead(name, bases.slice(0, cutPos(0)), quals.slice(0, cutPos(0)), readOrientation, umi)
    } else
      SequencingRead(name, bases, quals, readOrientation, umi)
  }

  /**
   * find the 'distance' between this read and another
   * @param read2 the second sequencing read to consider
   * @return a double, the mismatched bases, normalized by the *longer* of the two read lengths
   */
  def distance(read2: SequencingRead): Double =
    (bases.zip(read2.bases).map{case(b1,b2) => if (b1 == b2) 0 else 1}.sum.toDouble + math.abs(bases.length - read2.bases.length)) /
      math.max(bases.length,read2.bases.length).toDouble

  /**
   *
   * @return the average quality score value
   */
  def averageQual(): Double = intQuals.sum.toDouble / intQuals.length.toDouble

  /**
   * does the read contain the name primer?
   * @param primer the primer sequence (please make sure the reverse complement is done for reverse reads)
   * @param window the window of bases added to the primer length when searching the beginning of the read
   * @return true if the read starts with the primer
   */
  def startsWithPrimer(primer: String, window: Int = 5): Boolean = (bases.slice(0,primer.length + 5) contains primer)

  /**
   *
   * @return a string representing the read in fastq format
   */
  def toFastqString(umi: String, rev: Boolean): String = {
    val plusOrMinus = readOrientation match {
      case ForwardReadOrientation => "+"
      case ReverseReadOrientation => "+" // nevermind, this is not imporant
      case ReferenceRead => "+" // ehh lets assume
      case ConsensusRead => "+" // ehh lets assume
    }
    if (!rev)
      "@" + umi + "_" + name + "\n" + bases + "\n" + plusOrMinus + "\n" + quals
    else
      "@" + umi + "_" + name + "\n" + Utils.reverseComplement(bases) + "\n" + plusOrMinus + "\n" + quals
  }


  /**
   * find the first and last non-deletion character in a read -- useful when we align to the reference
   * and we want to remove cruft like reading into adapters, etc
   * @return a tuple pair of the first non-dash base and last non-dash base (zero-indexed)
   */
  def firstAndLastActualBases(): Tuple2[Int,Int] = {
    var firstNonDash = 0
    var lastNonDash = bases.length - 1


    while(firstNonDash < bases.length & bases(firstNonDash) == '-')
      firstNonDash += 1
    while(lastNonDash > 0 & bases(lastNonDash) == '-')
      lastNonDash -= 1
    (firstNonDash,lastNonDash)
  }

  /**
   * generate a cigar string from the read aligned another sequence (the reference most likely).  The difference
   * here is that we consider a unique situation -- if we find mismatches attached on either end to an indel, we roll those
   * bases into the indel event.  This hopefully will lower the diversity of indels we have in our final alignment
   *
   * @param otherRead the read to align to
   * @return a tuple of cigar elements and their lengths
   */
  def callCigarsAgainstOtherRead(otherRead: SequencingRead): Array[CigarHit] = {

    // find the first and last real letter of the other read -- most likely the reference read
    val fwdFirstAndLast = this.firstAndLastActualBases()

    // lets not assume they're the same size, assert that they are
    if (otherRead.length != this.length)
      throw new IllegalStateException("Unable to compare two reads of different sizes")

    var events = Array[CigarEvent]()
    //(0 until bases.length).foreach {position => {
    (fwdFirstAndLast._1 until fwdFirstAndLast._2 + 1).foreach {position => {
      if (bases(position) != '-' && otherRead.bases(position) != '-' && bases(position) == otherRead.bases(position))
        events :+= Match
      else if (bases(position) != '-' && otherRead.bases(position) != '-')
        events :+= Mismatch
      else if (bases(position) == '-' && otherRead.bases(position) != '-')
        events :+= Deletion
      else if (bases(position) != '-' && otherRead.bases(position) == '-')
        events :+= Insertion
      //else
      //  throw new IllegalStateException("WHY!!! " + bases(position) + " and " + otherRead.bases(position))
    }}


    // now track through the edits and fold in events
    var lastEvent: CigarEvent = Unset
    var currentEventSize = 0
    var currentMismatchesToClaim = 0
    var returnEvents = Array[CigarHit]()
    var startingPos = fwdFirstAndLast._1

    // large state machine -- spelled out to be as clear as I can
    events.foreach { evt => {
      evt match {
        case Match => {
          lastEvent match {
            case Match => {
              // match - match, just extend
              currentEventSize += (1 + currentMismatchesToClaim)
            }
            case Mismatch => {
              // convert the mismatch into a match, we would of added it to a indel if it terminated an indel event
              lastEvent = Match
              currentEventSize += (1 + currentMismatchesToClaim)
            }
            case Unset => {
              lastEvent = Match
              currentEventSize += 1
            }
            case _ => {
              // match - indel - push an indel event and reset to match
              returnEvents :+= CigarHit(lastEvent, currentEventSize, startingPos)
              startingPos += currentEventSize
              lastEvent = Match
              currentEventSize = 1
            }
          }
          currentMismatchesToClaim = 0 // reclaim all mismatches
        }

        case Mismatch => {
          lastEvent match {
            case Mismatch => {
              // one more mismatch to figure out later
            }
            case Match => {
              // start a tally of mismatches
              currentMismatchesToClaim += 1
            }
            case Unset => {
              lastEvent = Mismatch
              currentMismatchesToClaim += 1
            }
            case _ => {
              // add this mismatch to the previous event
              currentEventSize += 1
            }
          }
        }
        case Insertion => {
          lastEvent match {
            case Insertion => {
              currentEventSize += 1
            }
            case Match => {
              returnEvents :+= CigarHit(lastEvent, currentEventSize, startingPos)
              startingPos += currentEventSize
              lastEvent = Insertion
              currentEventSize = (1 + currentMismatchesToClaim)
            }
            case Mismatch => {
              // claim this mismatch into the current indel
              lastEvent = Insertion
              currentEventSize += (1 + currentMismatchesToClaim)
            }
            case Deletion => {
              // maybe raise a flag here, this is suspect.....
              returnEvents :+=CigarHit(lastEvent, currentEventSize, startingPos)
              startingPos += currentEventSize
              lastEvent = Insertion
              currentEventSize = 1
            }
            case Unset => {
              lastEvent = Insertion
              currentEventSize += 1
            }
          }
          currentMismatchesToClaim = 0
        }
        case Deletion => {
          lastEvent match {
            case Deletion => {
              currentEventSize += 1
            }
            case Match => {
              returnEvents :+=CigarHit(lastEvent, currentEventSize, startingPos)
              startingPos += currentEventSize
              lastEvent = Deletion
              currentEventSize = (1 + currentMismatchesToClaim)
            }
            case Mismatch => {
              // convert the mismatch into a match
              lastEvent = Deletion
              currentEventSize += (1 + currentMismatchesToClaim)
            }
            case Insertion => {
              returnEvents :+=CigarHit(lastEvent, currentEventSize, startingPos)
              startingPos += currentEventSize
              currentEventSize = 1
            }
            case Unset => {
              lastEvent = Deletion
              currentEventSize += 1
            }
          }
          currentMismatchesToClaim = 0
        }
        case Unset => {
          // do nothing
          currentMismatchesToClaim = 0
        }
      }
    }
    }
    if (currentEventSize > 0)
      returnEvents :+= CigarHit(lastEvent, currentEventSize + currentMismatchesToClaim, startingPos)

    //println(returnEvents.map{ty => ty.event.toString}.mkString(","))
    return returnEvents
  }
}


object SequencingRead {
  /**
   * produce a reverse complement of a read, taking some care to correspond the bases and quals
   * @param sequencingRead the sequencing read to reverse
   * @return a sequencingread object that represents the reverse read
   */
  def reverseComplement(sequencingRead: SequencingRead): SequencingRead = {
    //name: String, bases: String, quals: String, forwardRead: Boolean, umi: String
    val newBases = Utils.reverseComplement(sequencingRead.bases.toUpperCase())
    val newQuals = sequencingRead.quals.reverse

    SequencingRead(sequencingRead.name,newBases,newQuals,sequencingRead.readOrientation,sequencingRead.umi)
  }

  /**
   * strip the insertions out of a read, most likely for output
   * @param sequencingRead the input read
   * @return a sequencing read representation of the read with insertions stripped out
   */
  def stripDownToJustBases(sequencingRead: SequencingRead): SequencingRead = {
    //name: String, bases: String, quals: String, readOrientation: ReadDirection, umi: String
    var bases = ""
    var quals = ""
    for (i <- 0 until sequencingRead.bases.length) {
      if (sequencingRead.bases(i) != '-') {
        bases += sequencingRead.bases(i)
        quals += sequencingRead.quals(i)
      }
    }

    SequencingRead(sequencingRead.name,bases,quals,sequencingRead.readOrientation,sequencingRead.umi)
  }
}

object SequencingReadQualOrder extends Ordering[SequencingRead] {
  def compare(a:SequencingRead, b:SequencingRead) = a.averageQual() compare b.averageQual()
}

/**
 * the cigar event
 * @param encoding the encodings string
 */
sealed abstract class CigarEvent(var encoding: String = "U")
case object Unset extends CigarEvent("U")
case object Insertion extends CigarEvent("I")
case object Deletion extends CigarEvent("D")
case object Match extends CigarEvent("M")
case object Mismatch extends CigarEvent("m") // this isn't a valid output state, just for bookkeeping

/**
 * enumerate out the possible read orientations: forward, reverse, consensus, and reference
 */
sealed trait ReadDirection 
case object ForwardReadOrientation extends ReadDirection
case object ReverseReadOrientation extends ReadDirection
case object ReferenceRead extends ReadDirection
case object ConsensusRead extends ReadDirection

case class CigarHit(event: CigarEvent, length: Int, position: Int)