merged.data <- read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_02_10_Adult_Fish_7_9_12/test.txt",stringsAsFactors = F)


seven.b <- merged.data[merged.data$sample == "7B",]
seven.b <- seven.b[grep("to",as.character(seven.b$organ),invert=TRUE,value=FALSE),]

# heart reds
heart_colors = c("#FF8370","#FF2300","#D91D00","#A71700")

# eyes green
eye_colors = c("#00FF4D","#007B25")

# GI colors - blues
gi_colors = c("#08C9FF","#006987")

# brain color 
brain_color = c("brown")

# gills
gills_color = c("gray")

seven.b$color = "black"
seven.b[seven.b$organ=="Gills",]$color = gills_color[1]
seven.b[seven.b$organ=="Brain",]$color = brain_color[1]
seven.b[seven.b$organ=="Upper_GI",]$color = gi_colors[1]
seven.b[seven.b$organ=="Intestine",]$color = gi_colors[2]
seven.b[seven.b$organ=="Eye1",]$color = eye_colors[1]
seven.b[seven.b$organ=="Eye2",]$color = eye_colors[2]

seven.b[seven.b$organ=="Heart_chunk",]$color = heart_colors[1]
seven.b[seven.b$organ=="Heart_diss",]$color = heart_colors[2]
seven.b[seven.b$organ=="Heart_GFP-",]$color = heart_colors[3]
seven.b[seven.b$organ=="Heart_GFP+",]$color = heart_colors[4]
all_colors = c("Heart_chunk"="#A20FFF","Heart_diss"="#CB78FF","Heart_GFP-"="#590091",
               "Heart_GFP+"="#150070","Eye1"="#00FF4D","Eye2"="#007B25","Upper_GI"="#08C9FF","Intestine"="#006987",
               "Brain"="#D6FF00","Gills"="gray","Blood"="#FF2300")


# order by the total proportion that's going on our plot
df2 <- aggregate(proportion ~ event, data=seven.b, sum)
seven.b$event <-factor(seven.b$event, levels=df2[order(df2$proportion), "event"])


gg = ggplot(seven.b[seven.b$proportion >= 0.01,]) + 
  geom_bar(aes(x=event,y=proportion,fill=organ),col="black",stat="identity") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values=all_colors) + xlab("HMID") + ylab("Summed proportion over organs") + ggtitle("Common HMID sharing across organs")

ggsave(gg,file="histogram_of_common_alleles_fish_7B.pdf",width=20,height=15)
