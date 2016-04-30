# point to the  directory that hosts the .allReadCounts files
setwd("Directory Name")

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

# set up results display
tissues <- c("blood","heartchunk","heartdis","heartgfpm","heartgfpp","uppergi","intestine","gills","brain","eye1","eye2","all")
results <- data.frame(tissues,count = NA, prop = NA)
pattern <- ""

# make a grep-able version of the event string, replace the example with your event
pattern <- gsub("\\+", "\\\\+", "24D+136_NONE_NONE_9D+217_12D+244_87D+265_87D+265_87D+265_87D+265_NONE")

# for every tissue, show proportion and counts of this edit(s)
for (i in 1:12){
  results$prop[i] <- sum(all.list[[i]]$count[grep(pattern, all.list[[i]]$event)])/sum(all.list[[i]]$count)
  results$count[i] <- sum(all.list[[i]]$count[grep(pattern, all.list[[i]]$event)])
  }
results

# number of different HMIDs within this clade 
length(unique(all[grep(pattern, all$event),]$event))
