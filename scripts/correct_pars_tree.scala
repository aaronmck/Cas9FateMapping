import scala.io._
import java.io._
import scala.util.matching._
import scala.collection.mutable._

val fullStr = Source.fromFile(args(0)).getLines().mkString("")
val outputTree = new PrintWriter(args(1))

val lookup = new HashMap[String,Int]()
val tokenSplit = fullStr.split(":")
outputTree.write(tokenSplit.map{tk => {
  if (tk contains "_") {
    val sp = tk.split("_")
    val key = sp(1)
    lookup(key) = lookup.getOrElse(key,0) + 1
    sp(0) + "_" + key + lookup(key)
  } else {
    tk
  }
}}.mkString(":") + "\n")
outputTree.close()

val outputAnnotations = new PrintWriter(args(2))
outputAnnotations.write("taxa\tembryo\tumi\n")

val fullStr2 = Source.fromFile(args(1)).getLines().mkString("").split(":")
fullStr2.foreach{ln => {
  ln.split("\\(+").foreach{ln2 => {
    ln2.split(",").foreach{ln3 => {
      if (ln3 contains "_") {
        val sp = ln3.split("_")
        outputAnnotations.write(sp(0) + "_" + sp(1) + "\t" + sp(0) + "\t" + sp(0) + sp(1) + "\n")
      }
    }}
  }}
}}
outputAnnotations.close()
