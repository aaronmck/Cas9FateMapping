# ---------------
# load libraries
# ---------------

library(RColorBrewer)
library(clusterSim)
library(gplots)


# ---------------
# load data
# ---------------

# point to the  directory that hosts the .allReadCounts files
# setwd("Directory Name")
setwd("C:/Users/James/Documents/Postdoc Work/Lineage Tracing/Clades/2ndReAlignments")
setwd("C:/Users/James/Documents/Postdoc Work/Lineage Tracing/Clades/04012016_alignments")

# set fish ID number
fish <- "7B"
fish <- "17"

# initialize a couple vectors
tissues <- c("blood","heartchunk","heartdis","heartgfpm","heartgfpp","uppergi","intestine","gills","brain","eye1","eye2","all")
results <- data.frame(tissues,count = NA, prop = NA)

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
# remove blood alleles from all organs
# ---------------

# make a grep-compatible vector of all the blood alleles (deal with the '+' character)
blood.hmids <- gsub("\\+", "\\\\+", blood$event)

# remove the top blood alleles (default is 5)
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

# confirm that new proportions now sum correctly
unlist(lapply(all.list.noblood, function(x){sum(x$prop)}))    # old proportions
unlist(lapply(all.list.noblood, function(x){sum(x$newprop)})) # new proportions


# ---------------
# subset for dominant alleles where abundance > proportion threshold 
# ---------------

# first, figure out how many elements per organ to include
# using a threshold for what is a top/dominant allele, default is 5%

prop.threshold <- 0.05
number.elements <- unlist(lapply(all.list.noblood, function(x){sum(x[1:nrow(x),]$newprop > prop.threshold)}))

# second, find top alleles and their abundance in every organ

# sorry this is a bunch of nested for loops *cringe* - note # comments for explanation
# it works, but I should switch to apply() at some point for clarity

top.alleles <- c(heartchunk=NA,heartdis=NA,heartgfpm=NA,heartgfpp=NA,uppergi=NA,
                 intestine=NA,gills=NA,brain=NA,eye1=NA,eye2=NA)
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
top.alleles <- top.alleles[-1,] # remove top row of NAs


# ---------------
# set colors using RColorBrewer
# ---------------

cols <- colorRampPalette(brewer.pal(9,"Reds"))(1000)


# ---------------
# linear scale each allele to normalize for abundance differences between alleles
# ---------------

# normalize with type n8, x/max, so every row gets set to max of 1
top.alleles.scaled <- data.Normalization(top.alleles,type="n8",normalization="row")


# ---------------
# make a nice heatmap with key, using colors defined above
# ---------------

names <- c("piece of heart","DHCs","NCs","cardiomyocytes","intestinal bulb","post. intestine","gills","brain","left eye","right eye")
heatmap.2(t(top.alleles.scaled), col = cols, labRow = names, 
          Rowv=FALSE, Colv="Rowv",
          dendrogram="none", trace="none",
          margins = c(8,8), 
          key.title="", keysize=1.0, 
          key.xlab="scaled proportion",density.info="none",
          scale="none",
          lmat=rbind(3:4,2:1),lwid=c(1,4))


# ---------------
# output a cleaner version of top.alleles for browsing 
# ---------------

top.alleles.x <- cbind(name=top.alleles.events[-1], as.data.frame(top.alleles))
write.csv(top.alleles.x, file="dominant.alleles.csv")
