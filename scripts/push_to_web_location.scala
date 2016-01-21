import java.io.{File,FileInputStream,FileOutputStream}
import scala.sys.process._
import java.io._

// "scala -J-Xmx1g " + toWebPublishScript + " " + webLocation + " " + perBase + " " + topR + " " + topReadC + " " + allReadC

val webLocation = new File(args(0))
val perbaseFile = new File(args(1))
val topReadsFile = new File(args(2))
val topReadCountsFile = new File(args(3))
val cutSiteFile = new File(args(4))
val allReadFile = new File(args(5))

val javaScriptFile = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/plots/read_plot/read_editing_mutlihistogram.js")
val htmlFile = new File("/net/shendure/vol10/projects/CRISPR.lineage/nobackup/codebase/plots/read_plot/read_editing_mutlihistogram.html")

// check that everything exists
if (!webLocation.exists)
  throw new IllegalArgumentException("Unable to find location webLocation " + webLocation)

if (!webLocation.exists)
  throw new IllegalArgumentException("Unable to find location perbaseFile " + perbaseFile)

if (!webLocation.exists)
  throw new IllegalArgumentException("Unable to find location topReadsFile " + topReadsFile)

if (!webLocation.exists)
  throw new IllegalArgumentException("Unable to find location topReadCountsFile " + topReadCountsFile)

if (!javaScriptFile.exists)
  throw new IllegalArgumentException("Unable to find location javaScriptFile " + javaScriptFile)

if (!htmlFile.exists)
  throw new IllegalArgumentException("Unable to find location htmlFile " + htmlFile)

if (!cutSiteFile.exists)
  throw new IllegalArgumentException("Unable to find location cutsiteFile " + htmlFile)

if (!allReadFile.exists)
  throw new IllegalArgumentException("Unable to find location allReadFile " + allReadFile)

// copy a file
//def copyFile(src: File, dest:File) =new FileOutputStream(dest) getChannel() transferFrom(new FileInputStream(src) getChannel, 0, Long.MaxValue )
def copyToDir(inputFile: File, outputDir: File): Boolean = (("cp " + inputFile.getAbsolutePath + " " + outputDir.getAbsolutePath).! == 0)

if (!copyToDir(perbaseFile,webLocation))        throw new IllegalArgumentException("unable to copy " + perbaseFile + " to " + webLocation)
if (!copyToDir(topReadsFile,webLocation))       throw new IllegalArgumentException("unable to copy " + topReadsFile + " to " + webLocation)
if (!copyToDir(topReadCountsFile,webLocation))  throw new IllegalArgumentException("unable to copy " + topReadCountsFile + " to " + webLocation)
if (!copyToDir(htmlFile,webLocation))           throw new IllegalArgumentException("unable to copy " + htmlFile + " to " + webLocation)
if (!copyToDir(javaScriptFile,webLocation))     throw new IllegalArgumentException("unable to copy " + javaScriptFile + " to " + webLocation)
if (!copyToDir(cutSiteFile,webLocation))        throw new IllegalArgumentException("unable to copy " + cutSiteFile + " to " + webLocation)
if (!copyToDir(allReadFile,webLocation))        throw new IllegalArgumentException("unable to copy " + allReadFile + " to " + webLocation)

// now in the web directory make a little javascript file that tells the main JS which files to use
val fileList = new PrintWriter(webLocation + File.separator + "JS_files.js")
fileList.write("var occurance_file = \"" + topReadCountsFile.getName() + "\"\n")
fileList.write("var top_read_melted_to_base = \"" + topReadsFile.getName() + "\"\n")
fileList.write("var per_base_histogram_data = \"" + perbaseFile.getName() + "\"\n")
fileList.write("var cut_site_file = \"" + cutSiteFile.getName() + "\"\n")
fileList.close()
