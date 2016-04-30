# load libraries
library(colorspace)
library(gplots)

# point to the  directory that hosts the .allReadCounts files
# setwd("Directory Name")
setwd("C:/Users/James/Documents/Postdoc Work/Lineage Tracing/Clades/04012016_alignments")

# set fish ID number
fish <- "17"

# load .allReadCounts files
brain <- read.delim(paste(fish,"_Brain.allReadCounts.txt", sep=""),stringsAsFactors=F)
blood <- read.delim(paste(fish,"_Blood.allReadCounts.txt", sep=""),stringsAsFactors=F)
eye1 <- read.delim(paste(fish,"_Eye1.allReadCounts.txt", sep=""),stringsAsFactors=F)
eye2 <- read.delim(paste(fish,"_Eye2.allReadCounts.txt", sep=""),stringsAsFactors=F)
gills <- read.delim(paste(fish,"_Gills.allReadCounts.txt", sep=""),stringsAsFactors=F)
heartdis <- read.delim(paste(fish,"_Heart_diss.allReadCounts.txt", sep=""),stringsAsFactors=F)
heartchunk <- read.delim(paste(fish,"_Heart_chunk.allReadCounts.txt", sep=""),stringsAsFactors=F)
heartgfpm <- read.delim(paste(fish,"_Heart_GFP-.allReadCounts.txt", sep=""),stringsAsFactors=F)
heartgfpp <- read.delim(paste(fish,"_Heart_GFP+.allReadCounts.txt", sep=""),stringsAsFactors=F)
intestine <- read.delim(paste(fish,"_Intestine.allReadCounts.txt", sep=""),stringsAsFactors=F)
uppergi <- read.delim(paste(fish,"_Upper_GI.allReadCounts.txt", sep=""),stringsAsFactors=F)

# make all_organ data.frame and list
all <- rbind(blood,heartchunk,heartdis,heartgfpm,heartgfpp,uppergi,intestine,gills,brain,eye1,eye2)
all.list <-list(blood,heartchunk,heartdis,heartgfpm,heartgfpp,uppergi,intestine,gills,brain,eye1,eye2,all)

# throw away reads containing unknowns, should be obsolete now that we fixed this issue
all.list <- lapply(all.list, function(x) {
  x[grep("UNKNOWN", x$event, invert=TRUE),]
})

# ---------------
# remove blood HMIDs from all organs
# ---------------

# make a grep-compatible vector of all the blood HMIDs (deal with the '+' character)
blood.hmids <- gsub("\\+", "\\\\+", blood$event)

# remove the top blood alleles
num.alleles <- 5
blood.hmids.top <- paste(blood.hmids[1:num.alleles], collapse="|")
all.list.noblood <- lapply(all.list, function(x) {
  x[grep(blood.hmids.top, x$event, invert=TRUE),]
})

# ---------------
# shorten the list to remove blood data frame
# ---------------

all.list.noblood <- all.list.noblood[2:11]

# ---------------
# fix the proportions so that sum of proportions = 1 again, so that we can find top dominant alleles
# ---------------

# count the number of cells after removing blood
cells.by.organ <- unlist(lapply(all.list.noblood, function(x){sum(x[1:nrow(x),]$count)}))

# calculate new proportions
all.list.noblood.fixedprop <- lapply(all.list.noblood, function(this.organ){
  this.organ$newprop <- sapply(this.organ$count, function(this.hmid.count){this.hmid.count/sum(this.organ$count)})
})

# add new proportions back to list, rename
all.list.noblood <- Map(cbind, all.list.noblood, all.list.noblood.fixedprop)
for (x in 1:10){names(all.list.noblood[[x]])[5] <- "newprop"} # fix funky column name

# confirm that proportions now sum correctly
unlist(lapply(all.list.noblood, function(x){sum(x$prop)}))    # old proportions
unlist(lapply(all.list.noblood, function(x){sum(x$newprop)})) # new proportions

# ---------------
# subset for dominant alleles where abundance > proportion threshold 
# ---------------

# first, figure out how many elements per organ to include
prop.threshold <- 0.05
number.elements <- unlist(lapply(all.list.noblood, function(x){sum(x[1:nrow(x),]$newprop > prop.threshold)}))

# second, find top alleles and their abundance in every organ

# sorry this is a bunch of nested for loops *cringe* - note # comments for explanation
# it works, but I should switch to apply() at some point for clarity

top.alleles <- c(heartchunk=NA,heartdis=NA,heartgfpm=NA,heartgfpp=NA,uppergi=NA,intestine=NA,gills=NA,brain=NA,eye1=NA,eye2=NA)
top.alleles.events <- c(pattern=NA)

for (this.organ in 1:10){                                         # for every organ,
  for (this.abundant.hmid in 1:number.elements[this.organ]){      # for the most abundant HMIDs,
    pattern <- all.list.noblood[[this.organ]]$event[this.abundant.hmid]
    for (i in 1:10){                    # find abundance of that HMID in every organ  
      this.prop <- all.list.noblood[[i]]$newprop[grep(pattern, all.list.noblood[[i]]$event, fixed=TRUE)]
      if (identical(this.prop, numeric(0))){results$prop[i] <- 0}
      else {results$prop[i] <- this.prop[1]}
    }
    top.alleles <- rbind(top.alleles,results$prop[1:10]) # add those values to top.alleles
    top.alleles.events <- c(top.alleles.events, pattern)
  }
}

# flip vertically for plotting purposes
top.alleles.flip <- top.alleles[order(nrow(top.alleles):1),]

# ---------------
# set colors manually using colorspace palette, with parameters listed in comment
# ---------------

pal <- choose_palette() # set h1=10, C1=100, C2=0, L1=20, P1 = 2.7
cols2 <- rev(pal(1000))

# ---------------
# make a nice heatmap with legend, using colors defined above
# ---------------

heatmap.2(top.alleles.flip, Rowv=FALSE, Colv="Rowv",dendrogram="none", trace="none",
          col = cols2, 
          scale="row",density.info="none")

# ---------------
# output a cleaner version of top.alleles for browsing
# ---------------

top.alleles.x <- cbind(name=top.alleles.events, as.data.frame(top.alleles))
top.alleles.x <- top.alleles.x[-1,] # get rid of first row
write.csv(top.alleles.x, file="top.alleles.new.csv")
