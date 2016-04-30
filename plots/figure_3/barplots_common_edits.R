library(ggplot2)
library(ggthemes)
library(reshape2)


# shape the data down to usable columns and samples
embryo_uniques = read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_03_05_Embryos/embryos_sup_table.txt")
embryo_uniques = embryo_uniques[embryo_uniques$sample != "3d_1b_0.3x" & embryo_uniques$sample != "3d_1b_1x",]
embryo_uniques = embryo_uniques[,complete.cases(t(embryo_uniques))]
embryo_uniques = embryo_uniques[embryo_uniques$passHMIDs > 100 & embryo_uniques$uniqueHMIDs > 20,]

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

top_events = "/Users/aaronmck/google_drive/UW/shendure/CRISPR/2016_03_26_figure_3E/top_events_per_embryo.txt"

top_event_table = read.delim(top_events,stringsAsFactors=F)
top_event_table = top_event_table[order(top_event_table$propotions,decreasing = T),]
top_event_table$index = seq(1,nrow(top_event_table))

dev_colors = c("#F26E12","#0AA600","#0086FF","#515251")

top_event_table$sample[grep("_0",as.character(top_event_table$sample))] = 
  sapply(top_event_table$sample[grep("_0",as.character(top_event_table$sample))],
       function(x) {paste(x,".3x",sep="")})
top_event_table$rate = sapply(top_event_table$sample,split_to_rate)
top_event_table$stage = sapply(top_event_table$sample,split_to_hour)

idx <- match(top_event_table$sample,embryo_uniques$sample)
top_event_table$uniqueHMIDs <- embryo_uniques$Values[idx]

top_event_table$stage = factor(as.character(top_event_table$stage),levels=c("4.3HR","9HR","30HR","72HR"))

top_event_table = top_event_table[order(top_event_table$stage,top_event_table$propotions),]
top_event_table$sample = factor(top_event_table$sample,levels=embryo_uniques$sample)
top_event_table[top_event_table$type == "other",]$propotions = -1.0*top_event_table[top_event_table$type == "other",]$propotions

# plot by concentation and timepoint
prom_edits = ggplot(data=top_event_table) + 
  geom_bar(data=subset(top_event_table,type=="top"),aes(sample,propotions,fill=stage),stat="identity",width=.7) + 
  geom_bar(data=subset(top_event_table,type=="other"),aes(sample,propotions),fill="#A7A5A5",stat="identity",width=.7) + 
  scale_fill_manual(values=dev_colors) + 
  ylim(c(-0.25,0.60)) + scale_y_continuous(breaks=c(-.25, 0.0, 0.25, 0.50)) + 
  theme_tufte() + xlab("") + ylab("Dominant edit proportion") +
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) +
  theme(strip.text=element_text(family = "sans", size=12)) +
  facet_grid(. ~ rate,scales = "free", space = "free") + 
  # geom_hline(aes(yintercept=.125),color="white") +
  geom_hline(aes(yintercept=.25),color="white") +
  # geom_hline(aes(yintercept=-.125),color="white") +
  geom_hline(aes(yintercept=-.25),color="white") +
  geom_hline(aes(yintercept=0),color="black",size=1) +
  geom_hline(aes(yintercept=.50),color="black",linetype = "dashed") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  # theme(axis.text.x = element_blank()) +
  theme(legend.position="none")

ggsave(prom_edits,file="~/Desktop/proportion_of_events.pdf",width=8,height=4)
ggsave(prom_edits,file="~/Desktop/proportion_of_events.png",width=8,height=3)

