library(ggplot2)
library(ggthemes)
library(reshape2)

# all --> samples = c("7B_1_to_100_blood", "7B_1_to_20_blood","7B_1_to_500_blood","7B_Blood","7B_Brain","7B_Eye1","7B_Eye2","7B_Gills","7B_Heart_chunk","7B_Heart_diss","7B_Heart_GFP-","7B_Heart_GFP+","7B_Intestine","7B_Upper_GI")
samples = c("7B_Brain","7B_Eye1","7B_Eye2","7B_Gills","7B_Intestine","7B_Upper_GI","7B_Heart_chunk","7B_Heart_diss","7B_Heart_GFP-","7B_Heart_GFP+","7B_Blood")
new_samples = c("Brain","Eye #1","Eye #2","Gills","Intestine","Upper GI","Intact Heart","DHC","NCHC","Cardomyocytes","Blood")

# run against my mount location -- change as needed
hmid_base_location = "/mount/www/2016_02_10_Adult_Fish_7_9_12/"
hmids_per_sample = NA
for (i in seq(1,length(samples))) {
  sample = samples[i]
  new_sample = new_samples[i]
  
  filename = paste(hmid_base_location,"/",sample,"/",sample,".allReadCounts",sep="")
  print(filename)
  tmp = read.delim(filename)
  
  # setup some new columns
  tmp$sample = new_sample
  tmp$cprop = cumsum(tmp$proportion)
  
  if (is.na(hmids_per_sample)) {
    hmids_per_sample = tmp
  } else {
    hmids_per_sample = as.data.frame(rbind(hmids_per_sample,tmp))
  }
}

# now figure out the proportion of the top 5 blood alleles in other organs -- first get the top 5
blood_alleles = as.data.frame(hmids_per_sample[hmids_per_sample$sample == "Blood" & hmids_per_sample$count >= 3521,])
blood_alleles$event = factor(blood_alleles$event,levels=blood_alleles$event)

# for each organ, get the proportion of these blood alleles
blood_proportions_in_organs = data.frame(row.names=c(as.character(blood_alleles[,"event"]),"Non-blood"))

for (sample in new_samples) {
  blood_props = sapply(as.character(blood_alleles[,"event"]),function(event) {
    if (is.element(event,hmids_per_sample[hmids_per_sample$sample == sample,"event"])) {
      return (hmids_per_sample[hmids_per_sample$sample == sample & hmids_per_sample$event == event,"proportion"])
    } else {
      return (0.0)
    }
  })
  blood_props = c(blood_props,1.0 - sum(blood_props))
  blood_proportions_in_organs[,sample] = blood_props
}

blood_proportions_in_organs = as.data.frame(t(blood_proportions_in_organs))
colnames(blood_proportions_in_organs) = c(as.character(blood_alleles[,"event"]),"Non-blood")
blood_proportions_in_organs$sample = rownames(blood_proportions_in_organs)

blood_proportions_in_organs_melt = melt(id.vars=c("sample"),blood_proportions_in_organs)
colnames(blood_proportions_in_organs_melt) <- c("Sample","HMID","Proportion")

# set the ordering of the events so that ggplot displays it correctly with the corresponding color
# gray plus this palette: http://paletton.com/palette.php?uid=1000u0koptWdxFXjxwDtOqXwllJ
colors = c("#DDDDDD","#EF3939","#FF9393","#FF6363","#D70F0F","#AD0000")
  
blood_proportions_in_organs_melt$HMID = factor(as.character(blood_proportions_in_organs_melt$HMID),levels=c("Non-blood",as.character(blood_alleles[,"event"])))
blood_proportions_in_organs_melt$Sample = factor(as.character(blood_proportions_in_organs_melt$Sample),levels=c(new_samples))

blood_per_organ = ggplot(blood_proportions_in_organs_melt) + 
  geom_bar(aes(y=Proportion,x=Sample,fill=HMID),stat="identity") + 
  scale_fill_manual(values=colors,name="Concentration",guide=FALSE) + 
  theme_tufte() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.4, hjust=1)) +
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) 


ggsave(blood_per_organ,file="Blood_proportions_in_other_organs.png",width=5,height=3)