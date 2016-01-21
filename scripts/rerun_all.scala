import scala.sys.process._
import java.io._
import scala.io._
import java.util._
import scala.actors.Actor
import scala.actors.Actor._

// a list of directories to process
val baseDir = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup")
val directories = Array[File](
  new File("2015_05_10_Rerun_HEK_analysis/"),
  new File("2015_07_08_HEK4_exp3_redo/"),
  new File("2015_07_14_Target_3_analysis/"),
  new File("2015_07_26_barcoded_target1/"),
  new File("2015_10_06_Zebrafish_Initial_MiSeq/"),
  new File("2015_10_18_Full_TYR_sequencing/"),
  new File("2015_10_22_Cell_Culture_Lineage/"),
  new File("2015_10_31_OT_Deep_Sequencing/"),
  new File("2015_11_07_Dilution_Embryos/"),
  new File("2015_11_23_Deep_Sequence_Dilution/"),
  new File("2015_12_13_Adult_target1/"),
  new File("2015_12_18_Adult_target/"),
  new File("2015_12_29_V1_V2/"))

// This uses !! to get the whole result as a string
val dirContents = "ls".!!

val resetScript = "/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/scripts/reset_pipeline_output.py"

def setupCrisprPipeline(dir: File): String = {
  val base = "java -Xmx4g -jar /net/gs/vol1/home/aaronmck/tools/queue-protected/protected/gatk-queue-package-distribution/target/gatk-queue-package-distribution-3.5.jar"
  val script = " -S /net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/pipelines/CRISPR_analysis_PE_V2.scala "
  val tearsheet = " -i " + dir.getAbsolutePath + "/data/crispr_tearsheet.txt "
  val data = " --aggLocation " + dir.getAbsolutePath + "/data/pipeline_output/aggregate/"
  val expName = " --expName " + dir
  val runConditions = " -run -qsub -resMemReqParam mfree"
  val runscript = Array[String](base,script,tearsheet,data,expName,runConditions)
  return(runscript.mkString(" \\\n"))
}


class RunDirectory(dir: File) extends Actor {
  def act() {
    val tearsheet = new File(dir + File.separator + "data/crispr_tearsheet.txt")
    val run_script = new File(dir + File.separator + "bin/run_crispr_pipeline.sh")
    val output_dir = new File(dir + File.separator + "data/pipeline_output")
    val pipeline_output_out = new File(dir + File.separator + "data/output_of_run_" + Calendar.getInstance().getTime().toString.replace(" ","_") + "_out.txt")
    val pipeline_output_err = new File(dir + File.separator + "data/output_of_run_" + Calendar.getInstance().getTime().toString.replace(" ","_") + "_err.txt")

    val out = new StringBuilder
    val err = new StringBuilder

    val logger = ProcessLogger(
      (o: String) => out.append(o + "\n"),
      (e: String) => err.append(e + "\n"))

    if (tearsheet.exists) {
      //print("Reseting the output directory...")
      //val resetOutput = ("python " + resetScript + " --dir " + output_dir.getAbsolutePath).!
      //  println("result = " + resetOutput)
      val resetOutput = 0


      if (resetOutput == 0) {
        println("Rerunning directory " + dir)
        val runOutput = new PrintWriter(run_script)
        runOutput.write(setupCrisprPipeline(dir))
        runOutput.close()

        ("sh " + run_script) ! logger
        val output_results_err = new PrintWriter(pipeline_output_err)
        output_results_err.write(err.toString)
        output_results_err.close()

        val output_results_out = new PrintWriter(pipeline_output_out)
        output_results_out.write(out.toString)
        output_results_out.close()
      }
    }
    println("FINISHED with " + dir)
  }
}

// setup each of the run directories
directories.foreach{dir => {
  val newRun = new RunDirectory(dir)
  newRun.start
}}
