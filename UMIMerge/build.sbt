name := "UMIMerge"

version := "1.0"

scalaVersion := "2.11.7"

resolvers += Resolver.sonatypeRepo("public")

resolvers += "erichseifert.de" at "http://mvn.erichseifert.de/maven2"

unmanagedBase <<= baseDirectory { base => base / "project" }

libraryDependencies += "com.github.scopt" %% "scopt" % "3.2.0"

libraryDependencies += "org.scalatest" % "scalatest_2.11" % "2.2.4" % "test"

libraryDependencies += "org.scalatest" % "scalatest_2.11" % "2.2.4"

libraryDependencies += "org.apache.commons" % "commons-math3" % "3.5"

scalacOptions += "-target:jvm-1.7"


lazy val umi = (project in file(".")).
  settings(
      mainClass in (Compile, packageBin) := Some("main.scala.UMIProcessing"),
      mainClass in (Compile, run) := Some("main.scala.UMIProcessing")
  )


lazy val deepseq = (project in file(".")).
  settings(
    mainClass in (Compile, packageBin) := Some("main.scala.DeepSeq"),
      mainClass in (Compile, run) := Some("main.scala.DeepSeq"),
      assemblyJarName in assembly := "something.jar"
  )