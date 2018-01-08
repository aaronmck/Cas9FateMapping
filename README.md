An updated version of this codebase is available here: [Single-cell GESTALT codebase](https://github.com/aaronmck/SC_GESTALT)

# FateMappingCode
Collection of code from the [Shendure lab](https://github.com/shendurelab/) / [Schier lab](http://www.schierlab.fas.harvard.edu/) lineage tracing project using the CRISPR/Cas9 system. This is unsupported code, and stored here as an example of how we did our analysis in the paper, which may or may not be the way you'll want to do your own analysis.  The paper version of this repository is labeled v1, though some of the tree display code is newer than that (changes for our website).

# Abstract

Genome editing of synthetic target arrays for lineage tracing (GESTALT).

Multicellular systems develop from single cells through a lineage, but current lineage tracing approaches scale poorly to whole organisms. Here we use genome editing to progressively introduce and accumulate diverse mutations in a DNA barcode over multiple rounds of cell division. The barcode, an array of CRISPR/Cas9 target sites, records lineage relationships in the patterns of mutations shared between cells. In cell culture and zebrafish, we show that rates and patterns of editing are tunable, and that thousands of lineage-informative barcode alleles can be generated. We find that most cells in adult zebrafish organs derive from relatively few embryonic progenitors. Genome editing of synthetic target arrays for lineage tracing (GESTALT) will help generate large-scale maps of cell lineage in multicellular systems.

## Paper

The full paper can be found on [Science's website here](http://science.sciencemag.org/content/early/2016/05/25/science.aaf7907), and our pre-print is available on the [bioRxiv server here](http://biorxiv.org/content/early/2016/05/11/052712).

## Data

Data used in the paper is available in the [Dryad repository](http://datadryad.org/resource/doi:10.5061/dryad.478t9), which links out to GEO for the [sequencing data](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81713).

## Dependencies

There are lots (Sorry). For languages we rely upon Scala (Java/JVM), Python, Javascript/HTML, D3, and R generally (see the GitHub stats bar at the top).  For computational tools we used Trimmomatic, NEEDLEALL, MAFFT, GATK Queue, and Grid Engine to name a few.

## Directories

* ConstrainedNeighborJoin	- an old neighbor-joining approach we didn't use in the paper
* TreeSimulator	- again something that didn't make it into the paper, some intital work towards a simulator of events
* TreeUtils - the tool that takes Newick trees and annotations, and produces richly annotated JSON tree files for use in D3
* UMIMerge - code for merging and processing UMI tagged reads
* pipelines	- pipline code, dependent on the GATK Queue processing engine
* plots	- various plotting code, in either R, python, D3 (javascript + HTML), or Scala
* scripts - lots of scripts used throughout the processing and analysis
