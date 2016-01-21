import scala.io._
import java.io._

val inputFileList = Source.fromFile(args(0)).getLines()
val output = new PrintWriter(args(1))
output.write("fish\torgan\thmid\tindex\tcount\tproportion\n")
inputFileList.foreach{file => {
  val sample = file.split("\\/")(2)
  val fish = sample.slice(0,1)
  val organ = sample.slice(1,sample.length)
  val in = Source.fromFile(file).getLines().drop(1).foreach{line => {
    output.write(fish + "\t" + organ + "\t" + line + "\n")
  }}

}}
output.close()
