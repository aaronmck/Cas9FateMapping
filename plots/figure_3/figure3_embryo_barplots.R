library(ggplot2)
library(ggthemes)
library(reshape2)

# developmental time course simulations
# dev_colors = c("#515251","#0086FF","#0AA600","#F26E12") <- colors from the cdfs, reverse list here due to ordering demands
dev_colors = c("#F26E12","#0AA600","#0086FF","#515251")

# shape the data down to usable columns and samples
embryo_uniques = read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_05_04_embryo_rerun/2016_05_04_embryo_rerun_supplement.txt")
embryo_uniques = embryo_uniques[embryo_uniques$sample != "3d_1b_0.3x" & embryo_uniques$sample != "3d_1b_1x",]
embryo_uniques = embryo_uniques[,complete.cases(t(embryo_uniques))]
# embryo_uniques = embryo_uniques[embryo_uniques$passHMIDs > 100 & embryo_uniques$uniqueHMIDs > 20,]
embryo_uniques = embryo_uniques[embryo_uniques$passHMIDs > 150,]

# pull out the rate information
split_to_rate = function(x) {
  rate = unlist(strsplit(as.character(x),"_"))[3]
  if (rate == "1x") {
    return("1X")
  } else {
    return("1/3X")
  }
}
split_to_hour = function(x) {
  sample = unlist(strsplit(as.character(x),"_"))[1]
  if (sample == "30hr") {
    return("30HR")
  }
  if (sample == "3d") {
    return("72HR")
  }
  if (sample == "epi90") {
    return("9HR")
  }
  if (sample == "Dome") {
    return("4.3HR")
  }
  return("UNKNOWN")
}
embryo_uniques$rate = sapply(embryo_uniques$sample,split_to_rate)
embryo_uniques$stage = sapply(embryo_uniques$sample,split_to_hour)

# now sort by hour, then by unique counts
embryo_uniques$stage = factor(as.character(embryo_uniques$stage),levels=c("4.3HR","9HR","30HR","72HR"))
embryo_uniques = embryo_uniques[order(embryo_uniques$stage,embryo_uniques$uniqueHMIDs),]
embryo_uniques$index = seq(1,nrow(embryo_uniques))
embryo_uniques$sample = factor(embryo_uniques$sample,levels=embryo_uniques$sample)

# plot by concentation and timepoint
gg = ggplot(embryo_uniques) + 
  geom_bar(aes(sample,uniqueHMIDs,fill=stage),stat="identity",width=.7) + 
  scale_fill_manual(values=dev_colors) + 
  theme_tufte() + xlab("") + ylab("Unique alleles") +
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) +
  theme(strip.text=element_text(family = "sans", size=0)) +
  facet_grid(. ~ rate,scales = "free", space = "free") + 
  geom_hline(aes(yintercept=1000),color="white") +
  geom_hline(aes(yintercept=2000),color="white") +
  geom_hline(aes(yintercept=3000),color="white") +
  geom_hline(aes(yintercept=4000),color="white") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.position="none")
ggsave(gg,file="~/Desktop/unique_HMIDs_by_embryo.pdf",width=8,height=3)

# now sort by hour, then by unique counts
##embryo_uniques$stage = factor(as.character(embryo_uniques$stage),levels=c("4.3HR","9HR","30HR","72HR"))
#embryo_uniques = embryo_uniques[order(embryo_uniques$stage,embryo_uniques$passHMIDs),]
#embryo_uniques$index = seq(1,nrow(embryo_uniques))
#embryo_uniques$sample = factor(embryo_uniques$sample,levels=embryo_uniques$sample)

# plot passing
gg2 = ggplot(embryo_uniques) + 
  geom_bar(aes(sample,passHMIDs,fill=stage),stat="identity",width=.7) + 
  scale_fill_manual(values=dev_colors) + 
  theme_tufte() + xlab("") + ylab("HMIDs") +
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) +
  theme(strip.text=element_text(family = "sans", size=12)) +
  facet_grid(. ~ rate,scales = "free", space = "free") + 
  geom_hline(aes(yintercept=10000),color="white") +
  geom_hline(aes(yintercept=20000),color="white") +
  geom_hline(aes(yintercept=30000),color="white") +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="none")
ggsave(gg2,file="~/Desktop/total_HMIDs_by_embryo.pdf",width=8,height=3)


# now sort by hour, then by unique counts
#embryo_uniques$stage = factor(as.character(embryo_uniques$stage),levels=c("4.3HR","9HR","30HR","72HR"))
#embryo_uniques = embryo_uniques[order(embryo_uniques$stage,embryo_uniques$edited.HMID.prop),]
#embryo_uniques$index = seq(1,nrow(embryo_uniques))
#embryo_uniques$sample = factor(embryo_uniques$sample,levels=embryo_uniques$sample)

# plot by concentation and timepoint
gg4 = ggplot(embryo_uniques) + 
  geom_bar(aes(sample,meanSitesEdited,fill=stage),stat="identity",width=.7) + 
  scale_fill_manual(values=dev_colors) + 
  theme_tufte() + xlab("embryo") + ylab("Edited sites") +
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) +
  theme(strip.text=element_text(family = "sans", size=0)) +
  facet_grid(. ~ rate,scales = "free", space = "free") +
  theme(axis.text.x = element_blank()) +
  geom_hline(aes(yintercept=2),color="white") +
  geom_hline(aes(yintercept=4),color="white") +
  geom_hline(aes(yintercept=6),color="white") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(legend.position="none") +
  ylim(c(0,10))


ggsave(gg4,file="~/Desktop/mean_editing_rate.pdf",width=8,height=3)

