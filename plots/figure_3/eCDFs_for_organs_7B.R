# embryos <- read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_03_05_Embryos/embryos_merged.txt",stringsAsFactors = F)
sevenB <- read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_03_11_Adult_Fish_7_9_12_Reanalysis/merged_7B_allReads_march_16th_2016.txt",stringsAsFactors = F)

ggplot(embryos) + 
  geom_point(aes(array,running_prop,col=stage),size=1) + 
  theme_classic() + 
  scale_color_brewer(palette="Set2") + 
  facet_grid(.~ rate) + 
  xlab("HMIDs normalized to 1000 cells") +
  ylab("Cumulative density of HMIDs")

# make a master table
total_table = embryos[,c("stage","array","running_prop")]
colnames(total_table) <- c("sample","array","running_prop")                    
total_table = as.data.frame(rbind(total_table,sevenB[,c("sample","array","running_prop")]))


organ_cdfs = ggplot(sevenB) + 
  geom_line(aes(array,running_prop,col=organ),size=1) + 
  theme_tufte() + 
  scale_color_brewer(palette="Set2") +
  xlab("HMIDs normalized to 1000 cells") +
  ylab("Cumulative density of HMIDs") +
  theme(legend.text=element_text(family = "sans", size=20)) + 
  theme(legend.title=element_text(family = "sans", size=20)) + 
  theme(axis.text=element_text(family = "sans", size=20)) + 
  theme(axis.title=element_text(family = "sans", size=20))
ggsave(organ_cdfs,file="~/Desktop/organ_and_embryos_eCDF.png",width=6,height=4)