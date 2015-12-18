case class Alignment(val refPos: Int, refBase: String, readBase: String, cigarCharacter: String) {
  def combine(next: Alignment): Array[Alignment] =
    if (next.cigarCharacter != cigarCharacter)
      return Array[Alignment](this, next)
    else
      return Array[Alignment](Alignment(this.refPos, this.refBase + next.refBase, this.readBase + next.readBase, cigarCharacter))

  def prettyPrint: String = refPos + ":" + refBase + ":" + readBase + ":" + cigarCharacter

  def toEditString: String = readBase.length + "" + cigarCharacter + "+" + refPos + {
    if (cigarCharacter == "I") {
      "+" + readBase
    } else {
      ""
    }
  }
}

def callEdits(reference: String, read: String, minMatchOnEnd: Int, debugInfo: Boolean = false): List[Alignment] = {
    var referencePos = 0
    var inRef = false

    var refToEvent = List[Alignment]()

    reference.zip(read).foreach { case (refBase: Char, readBase: Char) =>
      //if (debugInfo)
      //  println("BASES: " + refBase + "," + readBase + " ")
      (refBase, readBase) match {
        case ('-', readB) if !inRef => {
          /* we might be in the situation where we haven't started the real alignment, take the offset */
        }
        case ('-', readB) if inRef => {
          // insertion

          if (refToEvent.isEmpty) refToEvent :+= Alignment(referencePos, refBase.toString, readBase.toString, "I")
          else refToEvent = refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos, refBase.toString, readBase.toString, "I"))

          if (debugInfo)
            println("1: " + refToEvent.map { st => st.prettyPrint }.mkString("<>") + " " + refToEvent.size)
        }
        case (refB, '-') if !inRef => {
          // deletion before read starts -- we haven't aligned yet
          referencePos += 1
        }
        case (refB, '-') => {
          // deletion
          inRef = true
          if (refToEvent.isEmpty) refToEvent :+= Alignment(referencePos, refBase.toString, readBase.toString, "D")
          else refToEvent = refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos, refBase.toString, readBase.toString, "D"))
          referencePos += 1
          if (debugInfo)
            println("2: " + refToEvent.map { st => st.prettyPrint }.mkString("<>") + " " + refToEvent.size)
        }
        case (refB, readB) => {
          // match / mismatch
          inRef = true
          if (refToEvent.isEmpty) refToEvent :+= Alignment(referencePos, refBase.toString, readBase.toString, "M")
          else refToEvent = refToEvent.init ++ refToEvent.last.combine(Alignment(referencePos, refBase.toString, readBase.toString, "M"))
          referencePos += 1
          if (debugInfo)
            println("3: " + refToEvent.map { st => st.prettyPrint }.mkString("<>") + " " + refToEvent.size)
        }


      }
    }
    // get a bit aggressive here -- start from both ends -- strip off insertions and deletions until we hit a match or mismatch of at least 10 bases
    return filterEnds(refToEvent, minMatchOnEnd)
  }

  def editsToCutSiteCalls(edits: List[List[Alignment]], cutSites: CutSites, debug: Boolean = false): Tuple2[Boolean, Array[String]] = {
      var ret = Array[String]()
      var retCovered = Array[Boolean]()
      var collision = false

      cutSites.windows.foreach { case (start, cut, end) => {
        var candidates = Array[Alignment]()


        edits.foreach { editList =>
          editList.foreach { edit => {
            if ((edit.cigarCharacter == "D" && overlap(start, end, edit.refPos, edit.refPos + edit.refBase.length)) ||
              edit.cigarCharacter == "I" && overlap(start, end, edit.refPos, edit.refPos + 1))
              candidates :+= edit
          }
          }
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
          else {
            ret :+= candidates(0).toEditString + "&" + candidates(1).toEditString
            collision = true
          }
        } else {
          ret :+= candidates.map { t => t.toEditString }.mkString("&")
          collision = true
        }
      }
      }

      return (collision, ret)
    }

  def percentMatch(ref: String, read: String, minimumAlignedBases: Int = 50): Double = {
    var bases = 0
    var matches = 0
    ref.zip(read).foreach { case (refBase, readBase) => {
      if (refBase != '-' && readBase != '-') {
        if (refBase == readBase)
          matches += 1
        bases += 1
      }
    }
    }

    if (bases < minimumAlignedBases)
      return -1.0
    matches.toDouble / bases.toDouble
  }



val events1 = AlignmentManager.callEdits(alignmentsF(0).bases, alignmentsF(1).bases, minMatchOnEnd, false)
val events2 = AlignmentManager.callEdits(alignmentsR(0).bases, alignmentsR(1).bases, minMatchOnEnd, false)

val combined = editsToCutSiteCalls(List[List[Alignment]](events1, events2), cutSites, debug)

val matchRate1 = percentMatch(alignmentsF(0).bases, alignmentsF(1).bases)
val matchRate2 = percentMatch(alignmentsR(0).bases, alignmentsR(1).bases)
