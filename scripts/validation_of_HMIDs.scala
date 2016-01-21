import scala.io._
import java.io._
import scala.collection.mutable._
import scala.math._

val orginal_dir = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2015_12_18_Adult_target/data/pipeline_output/"
val validation_dir = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2016_01_15_adult_fish_2_and_3/data/pipeline_output/"
val validation_table = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/2016_01_17_validation_of_adult_fish/data/validation_dirs.txt"

val output_per_organ = new PrintWriter("per_organ_overlap.txt")
val output_per_hmid = new PrintWriter("per_hmid_per_organ.txt")

// validation results for each sample
case class ValidationPackage(fullPath: String, subDir: String) {
  var statsReads = 0
  val hmids = new HashMap[String,Int]()

  val inputFile = new File(fullPath + "/" + subDir + "/" + subDir.stripSuffix("/") + ".stats")
  val headerTokens = (1 until 11).map {tk => "target" + tk}.toArray

  var first = true
  if (inputFile.exists()) {
    val inputFl = Source.fromFile(inputFile).getLines()
    val header = inputFl.next().split("\t")
    val targetPositions = headerTokens.map{tk => header.indexOf(tk)}.toArray

    inputFl.foreach{line => {
      if ((line contains "PASS") && !(line contains "WT") && !(line contains "UNKNOWN")) {
        statsReads += 1
        val sp = line.split("\t")
        val target = targetPositions.map{pos => sp(pos).split("\\+")(0)}.mkString("-")
        if (first) {
          //println(target)
          first = false
        }
        hmids(target) = hmids.getOrElse(target,0) + 1
      }
    }}
  }
  //println("loaded " + hmids.size + " hmids from " + inputFile)
}

// val validation_pairs = new HashMap[String,Tuple3[ValidationPackage,ValidationPackage,ValidationPackage]]()

// read in each of the validation files and setup the validation package results
output_per_organ.write("fish\torgan\t")
output_per_organ.write("oLen\toIn1\toIn2\toInB\toInBO\t")
output_per_organ.write("oLenC\toIn1C\toIn2C\toInBC\toInBCO\t")

output_per_organ.write("va1len\tval1InO\tval1InVal2\tval1InB\tval1InBO\t")
output_per_organ.write("va1lenC\tval1InOC\tval1InVal2C\tval1InBC\tval1InBCO\t")

output_per_organ.write("val2Len\tval2InOrig\tval2InVal1\tval2InB\tcal2InBO\t")
output_per_organ.write("val2LenC\tval2InOrigC\tval2InVal1C\tval2InBC\tval2InBCO\n")

output_per_hmid.write("organ\tfish\torganfish\thmid\toriginal_count\tvalidation1_count\tvalidation2_count\torig_prop\tval1_prop\tval2_prop\n")

Source.fromFile(validation_table).getLines().drop(1).foreach{line => {
  val tokens = line.split("\t")
  val fish = tokens(0)
  val organ = tokens(1).replace(" ","_")
  val org_dir = tokens(2)
  val val1_dir = tokens(3)
  val val2_dir = tokens(4)

  val orig = ValidationPackage(orginal_dir, org_dir)
  val val1 = ValidationPackage(validation_dir, val1_dir)
  val val2 = ValidationPackage(validation_dir, val2_dir)

  var master_hmid_set = orig.hmids.keys.toSet ++ val1.hmids.keys.toSet ++ val2.hmids.keys.toSet
  master_hmid_set.foreach{hmid =>
    output_per_hmid.write(organ + "\t" + fish + "\t" +
      fish + organ + "\t" + hmid + "\t" + orig.hmids.getOrElse(hmid,0) +
      "\t" + val1.hmids.getOrElse(hmid,0) + "\t" + val2.hmids.getOrElse(hmid,0) + "\t" +
      (orig.hmids.getOrElse(hmid,0)/orig.statsReads.toDouble) + "\t" + (val1.hmids.getOrElse(hmid,0)/val1.statsReads.toDouble) +
      "\t" + (val2.hmids.getOrElse(hmid,0)/val2.statsReads.toDouble) + "\n")
      
  }

  // check what proportion of the original HMIDs are seen in the validation sets
  val origInVal1 = orig.hmids.map{case(hmid,count) => if (val1.hmids contains hmid) 1 else 0}.sum
  val origInVal2 = orig.hmids.map{case(hmid,count) => if (val2.hmids contains hmid) 1 else 0}.sum
  val origInValB = orig.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) && (val1.hmids contains hmid)) 1 else 0}.sum
  val origInValBO = orig.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) || (val1.hmids contains hmid)) 1 else 0}.sum

  val origInVal1C = orig.hmids.map{case(hmid,count) => if (val1.hmids contains hmid) count else 0}.sum
  val origInVal2C = orig.hmids.map{case(hmid,count) => if (val2.hmids contains hmid) count else 0}.sum
  val origInValBC = orig.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) && (val1.hmids contains hmid)) count else 0}.sum
  val origInValBCO = orig.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) || (val1.hmids contains hmid)) count else 0}.sum

  val val1InOrig = val1.hmids.map{case(hmid,count) => if (orig.hmids contains hmid) 1 else 0}.sum
  val val1InVal2 = val1.hmids.map{case(hmid,count) => if (val2.hmids contains hmid) 1 else 0}.sum
  val val1InB = val1.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) && (orig.hmids contains hmid)) 1 else 0}.sum
  val val1InBO = val1.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) || (orig.hmids contains hmid)) 1 else 0}.sum

  val val1InOrigC = val1.hmids.map{case(hmid,count) => if (orig.hmids contains hmid) count else 0}.sum
  val val1InVal2C = val1.hmids.map{case(hmid,count) => if (val2.hmids contains hmid) count else 0}.sum
  val val1InBC = val1.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) && (orig.hmids contains hmid)) count else 0}.sum
  val val1InBCO = val1.hmids.map{case(hmid,count) => if ((val2.hmids contains hmid) || (orig.hmids contains hmid)) count else 0}.sum

  val val2InOrig = val2.hmids.map{case(hmid,count) => if (orig.hmids contains hmid) 1 else 0}.sum
  val val2InVal1 = val2.hmids.map{case(hmid,count) => if (val1.hmids contains hmid) 1 else 0}.sum
  val val2InB = val2.hmids.map{case(hmid,count) => if ((val1.hmids contains hmid) && (orig.hmids contains hmid)) 1 else 0}.sum
  val val2InBO = val2.hmids.map{case(hmid,count) => if ((val1.hmids contains hmid) || (orig.hmids contains hmid)) 1 else 0}.sum

  val val2InOrigC = val2.hmids.map{case(hmid,count) => if (orig.hmids contains hmid) count else 0}.sum
  val val2InVal1C = val2.hmids.map{case(hmid,count) => if (val1.hmids contains hmid) count else 0}.sum
  val val2InBC = val2.hmids.map{case(hmid,count) => if ((val1.hmids contains hmid) && (orig.hmids contains hmid)) count else 0}.sum
  val val2InBCO = val2.hmids.map{case(hmid,count) => if ((val1.hmids contains hmid) || (orig.hmids contains hmid)) count else 0}.sum

  output_per_organ.write(fish + "\t" + organ.replace(" ","_") + "\t")
  output_per_organ.write(orig.hmids.size + "\t" + origInVal1 + "\t" + origInVal2 + "\t" + origInValB + "\t" + origInValBO + "\t")
  output_per_organ.write(orig.hmids.map{case(k,c) => c}.sum + "\t" + origInVal1C + "\t" + origInVal2C + "\t" + origInValBC + "\t" + origInValBCO + "\t")

  output_per_organ.write(val1.hmids.size + "\t" + val1InOrig + "\t" + val1InVal2 + "\t" + val1InB + "\t" + val1InBO + "\t")
  output_per_organ.write(val1.hmids.map{case(k,c) => c}.sum + "\t" + val1InOrigC + "\t" + val1InVal2C + "\t" + val1InBC + "\t" + val1InBCO + "\t")

  output_per_organ.write(val2.hmids.size + "\t" + val2InOrig + "\t" + val2InVal1 + "\t" + val2InB + "\t" + val2InBO + "\t")
  output_per_organ.write(val2.hmids.map{case(k,c) => c}.sum + "\t" + val2InOrigC + "\t" + val2InVal1C + "\t" + val2InBC + "\t" + val2InBCO + "\n")

}}


output_per_organ.close()
output_per_hmid.close()
