import scala.io._
import java.io._

import java.lang.ProcessBuilder
import scala.concurrent._
import scala.collection.JavaConversions._
import java.io.{File,FileInputStream,FileOutputStream}

def runMix(inputMixFile: String, weightsFile: String, dirToRunIn: String): Unit = {
  val mixprogram = List[String]("mix")
  val pb = new ProcessBuilder(mixprogram)
  pb.directory(new File(dirToRunIn))
  pb.redirectOutput(ProcessBuilder.Redirect.INHERIT);
  val process = pb.start()

  val stdin = process.getOutputStream();
  val writer = new BufferedWriter(new OutputStreamWriter(stdin))

  val a = Array(inputMixFile, "P","W","4","5","Y",weightsFile).foreach { s =>
    writer.write((s + "\n"))
    writer.flush()
  }

  process.waitFor()

  writer.close()
}

// copy the mix file
val src = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/bin/phylip-fast/phylip-3.696/exe/mix")
val dest = new File(args(2) + "/mix")
println(dest.getAbsolutePath)
dest.createNewFile
new FileOutputStream(dest) getChannel() transferFrom( new FileInputStream(src) getChannel, 0, Long.MaxValue )

// test run
runMix(args(0),args(1), args(2))
