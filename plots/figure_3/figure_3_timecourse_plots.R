library(reshape2)
library(ggplot2)
library('ggthemes')

txt <- read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_03_05_Embryos/embryo_summary.txt",stringsAsFactors = F)

txt = txt[txt$passHMIDs > 100,]
txt$type <- sapply(txt$sample,function(x) {unlist(strsplit(x,"_"))[1]})
txt$rate <- sapply(txt$sample,function(x) {unlist(strsplit(x,"_"))[3]})
txt$type = factor(txt$type,levels=c("Dome","epi90","30hr","3d"))
txt = txt[order(txt$type,txt$meanSitesEdited),] # 
txt$index = seq(1,nrow(txt))

melt.txt = melt(txt,measure.vars = c("meanSitesEdited","meanCutSites"),id.vars = c("sample","type","index","rate"))

txt$type = as.character(txt$type)
txt$type[txt$type == "Dome"] = "4.3 hpf"
txt$type[txt$type == "epi90"] = "9 hpf"
txt$type[txt$type == "30hr"] = "30 hpf"
txt$type[txt$type == "3d"] = "72 hpf"

txt$type = factor(txt$type,levels=c("4.3 hpf","9 hpf","30 hpf","72 hpf"))

# color_set = c("#FF5F00","#0107FA")
# color_set = c("#999999", "#E69F00", "#56B4E9")
color_set = c("#999999", "#01A4AC","black")

# we dont have a ton of points, so maybe this isn't the best representation.  Turn on jitter to see
editing.by.embryo = ggplot(txt) + # geom_jitter(position=position_dodge(1)) + 
  scale_fill_manual(values=color_set,name="Concentration")  + theme_tufte() + 
  geom_jitter(aes(y=meanSitesEdited/10.0, x=type, fill=rate), colour="black",pch=21, size=3, position=position_dodge(1)) + 
  ylab("Edited proportion") + xlab("Developmental stage") + 
  theme(axis.text=element_text(family = "sans", size=15)) + 
  theme(axis.title=element_text(family = "sans", size=15)) + ylim(c(0,1)) +
  geom_smooth(method = "lm", se = FALSE)

editing.by.embryo = ggplot(txt) + # geom_jitter(position=position_dodge(1)) + 
  scale_fill_manual(values=color_set,name="Concentration")  + theme_tufte() + 
  geom_bar(aes(y=meanSitesEdited/10.0,x=sample, fill=rate),stat="identity") + 
  ylab("Edited proportion") + xlab("Developmental stage") + 
  theme(axis.text=element_text(family = "sans", size=15)) + 
  theme(axis.title=element_text(family = "sans", size=15)) + ylim(c(0,1)) +
  geom_smooth(method = "lm", se = FALSE)

unique.by.embryo = ggplot(txt) + # geom_jitter(position=position_dodge(1)) + 
  scale_fill_manual(values=color_set,name="Concentration")  + theme_tufte() + 
  geom_jitter(aes(y=uniqueHMIDs/passHMIDs, x=type, fill=rate), colour="black",pch=21, size=3, position=position_dodge(.3)) + 
  ylab("Unique HMID proportion") + xlab("Developmental stage") + 
  theme(axis.text=element_text(family = "sans", size=15)) + 
  theme(axis.title=element_text(family = "sans", size=15)) + ylim(c(0,1)) +
  geom_smooth(method = "lm", se = FALSE)

# ggsave(editing.by.embryo,file="embryo_editing_by_type.pdf",width=5,height=3)
ggsave(editing.by.embryo,file="embryo_editing_by_type.png",width=5,height=3,dpi = 300)
ggsave(unique.by.embryo,file="unique_editing_by_type.png",width=5,height=3,dpi = 300)

unique.by.embryo = ggplot(txt) + 
  geom_point(aes(x=passHMIDs, y=uniqueHMIDs, col=rate, shape=type), size=3) + 
  geom_rangeframe() + 
  theme_tufte() + 
  scale_shape_manual(values = c(15, 16, 17, 18),name="Developmental\nstage") + 
  scale_color_manual(values=color_set,name="Concentration") + 
  ylab("Unique HMIDs") + 
  xlab("Total captured HMIDs") + 
  geom_line(aes(x=passHMIDs, y=uniqueHMIDs/passHMIDs, col="black"), size=1) + 
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) 

unique.by.embryo = ggplot(txt) + 
  geom_bar(aes(x=passHMIDs, y=uniqueHMIDs, col=rate, shape=type), size=3) + 
  geom_rangeframe() + 
  theme_tufte() + 
  scale_shape_manual(values = c(15, 16, 17, 18),name="Developmental\nstage") + 
  scale_color_manual(values=color_set,name="Concentration") + 
  ylab("Unique HMIDs") + 
  xlab("Total captured HMIDs") + 
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12))

# ggsave(unique.by.embryo,file="unqiue_vs_captured_HMIDs_by_stage.pdf",width=5,height=3)
ggsave(unique.by.embryo,file="unqiue_vs_captured_HMIDs_by_stage.png",width=5,height=3)
