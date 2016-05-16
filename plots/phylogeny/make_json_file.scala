import java.io._
import scala.io._

def getListOfFiles(dir: String):List[File] = {
  val d = new File(dir)
  if (d.exists && d.isDirectory) {
    d.listFiles.filter(_.isFile).toList
  } else {
    List[File]()
  }
}


def name_converter(inputName:String):String = {
  val sp = inputName.split("_")
  inputName match {
    case x if x startsWith "epi" => {
      if (sp(2) startsWith "0")
        "9hpf_embryo_V7_number_" + sp(1) + "_one_third_X_cas9" 
      else
        "9hpf_embryo_V7_number_" + sp(1) + "_one_X_cas9"
    }
    case x if x startsWith "embryos" => {
      "30hpf_embryo_V6_number_" + sp(2) + "_one_X_cas9"
    }
    case x if x startsWith "Dome" => {
      if (sp(2) startsWith "0")
        "4.3hpf_embryo_V7_number_" + sp(1) + "_one_third_X_cas9"
      else
        "4.3hpf_embryo_V7_number_" + sp(1) + "_one_X_cas9"
    }
    case x if x startsWith "30hr" => {
      if (sp(2) startsWith "0")
        "30hpf_embryo_V7_number_" + sp(1) + "_one_third_X_cas9"
      else
        "30hpf_embryo_V7_number_" + sp(1) + "_one_X_cas9"
    }
    case x if x startsWith "3d" => {
      if (sp(2) startsWith "0")
        "72hpf_embryo_V7_number_" + sp(1) + "_one_third_X_cas9"
      else
        "72hpf_embryo_V7_number_" + sp(1) + "_one_X_cas9"
    }
    case x if x startsWith "tree" => {
      if (sp(1) == "7B")
        "ADR1"
      else
        "ADR2"
    }
    case _ => inputName
  }
}

val outputFile = new PrintWriter("master_list.json")
outputFile.write("[\n")
var outputs = Array[String]()
val width = 800

getListOfFiles("./tree_data/").filter(fl => fl.getAbsolutePath endsWith ".json").foreach{jsonFile => {
  
  val nodeCount = Source.fromFile(jsonFile).getLines().filter(ln => ln contains "name").size
  val averageDepth = Source.fromFile(jsonFile).getLines().filter(ln => ln contains "rootDist").map{ln => ln.split(":")(1).split(",")(0).trim.toDouble}.toArray
  val avg = averageDepth.sum / averageDepth.size.toDouble
  var height = 5 * (nodeCount + 50)
  outputs :+= "{\"name\" :\"" + name_converter(jsonFile.getName) + "\", \"tree_file\" :\"" + jsonFile + "\",\"height\":" + height + ", \"barheight\": " + 5 + ", \"width\": " + width + "}"
  println(jsonFile + "\t" + avg)
}}
outputFile.write(outputs.mkString(",\n"))
outputFile.write("]\n")
outputFile.close()

// var adult_17_NJ = {tree_file:"tree_data/output_tree_adult_17_proportional.json",width:dims.width, height:3000, barheight: 5}
