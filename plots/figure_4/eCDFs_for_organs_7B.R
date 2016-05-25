library(ggplot2)
library(ggthemes)
library(reshape2)
library(data.table)
library(grid)
library(gridExtra)

sevenB <- read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_04_08_Adult_Fish_7_9_12/merged_adult_7B_May_8th_2016_noblood_with_all.txt",stringsAsFactors = F)
sevenB$experiment = "adult"

# use this if we want to keep the median (we selected the lower-count median sample if there was a tie)
embryos_to_keep = c("3d_3_1x","30hr_3_0","epi90_2_1x","Dome_3_1x")
# if we want the most, use the following line
#embryos_to_keep = c("3d_1_1x","30hr_3_1x","epi90_12_1x","Dome_10_1x")

embryos <- read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_05_04_embryo_rerun/merged_v7_embryo_data_gt_100_May_5th_2016.txt",stringsAsFactors = F)
embryos = embryos[is.element(embryos$sample,embryos_to_keep),]


embryos <- embryos[,c("sample","stage","event","array","count","proportion","running_prop")]
embryos$experiment = "embryo"

# rename columns so that we can rbind it all together
colnames(embryos) <- c("sample","organ","event","array","count","proportion","running_prop","experiment")
embryos <- embryos[,c("sample","organ","event","array","count","proportion","running_prop","experiment")]

# merge it all together
total.data <- rbind(embryos,sevenB)

# ---------------------------------------------------------------------------------
# setup the correct names and colors for all the samples 
# ---------------------------------------------------------------------------------
total.data$organ[total.data$organ == "4.3HR"] = "4.3HR"
total.data$organ[total.data$organ == "9HR"] = "9HR"
total.data$organ[total.data$organ == "30HR"] = "30HR"
total.data$organ[total.data$organ == "72HR"] = "72HR"
total.data$organ[total.data$organ == "Blood"] = "Blood"
total.data$organ[total.data$organ == "Brain"] = "Brain" 
total.data$organ[total.data$organ == "Eye1"] = "Left eye"
total.data$organ[total.data$organ == "Eye2"] = "Right eye" 
total.data$organ[total.data$organ == "Gills"] = "Gills"
total.data$organ[total.data$organ == "Heart_chunk"] = "Heart"
total.data$organ[total.data$organ == "Heart_diss"] = "DHC" 
total.data$organ[total.data$organ == "Heart_GFP-"] = "NC" 
total.data$organ[total.data$organ == "Heart_GFP+"] = "Cardiomyocytes" 
total.data$organ[total.data$organ == "Intestine"] = "Post. intestine"
total.data$organ[total.data$organ == "Upper_GI"] = "Intestinal bulb"
total.data$organ[total.data$organ == "ALL"] = "All organs"

organColors = c("4.3HR" = "#CECECE",
                "9HR" = "#A2A2A2",
                "30HR" = "#767676",
                "72HR" = "#403F3F", 
                "Blood" = "#FF0000", 
                "Brain" = "#4F6128", 
                "Left eye" = "#77933C", 
                "Right eye" = "#C3D69B", 
                "Gills" = "#FFC000",
                "Heart" = "#632523", 
                "DHC" = "#943735", 
                "Cardiomyocytes" = "#E6B9B8", 
                "NC" = "#D99795", 
                "Intestinal bulb" = "#558ED5", 
                "Post. intestine" = "#8EB3E3",
                "All organs" = "black")

# set the text size
text.size = 16

# now actually plot the organ as a eCDF using the step function to be more clear that our data is discrete
organ_cdfs_focused = ggplot() + 
  geom_step(data=subset(total.data,experiment == "adult"),aes(array,running_prop,col=organ,group=organ),direction="vh",size=1.5) + 
  geom_step(data=subset(total.data,experiment == "embryo"),aes(array,running_prop,col=organ,group=sample),alpha=0.8, direction="vh",size=1.5,linetype=5) + 
  theme_tufte() + 
  scale_color_manual(values=organColors) + xlim(c(0,100)) + 
  xlab("Alleles") +
  ylab("Cumulative density of HMIDs") +
  theme(legend.text=element_text(family = "sans",  size=text.size)) + 
  theme(legend.title=element_text(family = "sans", size=text.size)) + 
  theme(axis.text=element_text(family = "sans",    size=text.size)) + 
  theme(axis.title=element_text(family = "sans",   size=text.size)) +
  theme(legend.position="none") 

ggsave(organ_cdfs_focused,file="~/Desktop/organ_focused_eCDF_fish7B.png",width=6,height=6)


#
#organ_cdfs_full = ggplot() + 
#  geom_step(data=subset(total.data,experiment == "adult"),aes(array,running_prop,col=organ,group=organ),direction="vh",size=1.5) + 
#  geom_step(data=subset(total.data,experiment == "embryo"),aes(array,running_prop,col=organ,group=sample),alpha=0.8, direction="vh",size=1.5) + 
#  theme_tufte() + 
#  scale_color_manual(values=organColors) +
#  xlab("Alleles") +
#  ylab("Cumulative density of HMIDs") +
#  theme(legend.text=element_text(family = "sans", size=text.size)) + 
#  theme(legend.title=element_text(family = "sans", size=text.size)) + 
#  theme(axis.text=element_text(family = "sans", size=text.size)) + 
#  theme(axis.title=element_text(family = "sans", size=text.size)) + 
#  theme(plot.background = element_rect(colour = "white")) +
#  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
#  theme(legend.position="none")

#vp <- viewport(width = 0.4, height = 0.5, x = .55,
#               y = .35)

#full <- function() {
#  print(organ_cdfs_focused)
#  theme_set(theme_bw(base_size = 8))
#  theme_classic()
#  print(organ_cdfs_full, vp = vp)
#  theme_set(theme_bw())
#}
#full()

