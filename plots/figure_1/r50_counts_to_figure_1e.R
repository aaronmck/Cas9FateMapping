# make the deletion and insertion size plots for the t1high (V1) data
library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)

events <- read.delim("event_histogram_r50high.txt")

blues <- rev(RColorBrewer::brewer.pal(3, "Blues"))
reds <- rev(RColorBrewer::brewer.pal(8, "Reds"))

events <- events[events$event.count < 10,]

glog = ggplot() + 
  geom_bar(data=events[events$type == "D",], aes(x=bin, y=-1.0 * log(-1.0 * count), fill=interaction(type,event.count)),size=0.7,col="black", stat="identity", position="stack") + 
  geom_bar(data=events[events$type == "I",], aes(x=bin, y=log(count), fill=interaction(type,event.count)),size=0.7,col="black", stat="identity", position="stack") + 
  theme_classic() +
  scale_fill_manual(values = c(reds,blues)) + 
  ylab("Log occurance") + 
  xlab("Insertion or deletion size")
ggsave(glog,file="log_edit_sizes_r50high.pdf",width=15,height=6)

gnorm = ggplot() + 
  geom_bar(data=events[events$type == "D",], aes(x=bin, y=count, fill=interaction(type,event.count)),size=0.7, col="black", stat="identity", position="stack") + 
  geom_bar(data=events[events$type == "I",], aes(x=bin, y=count, fill=interaction(type,event.count)),size=0.7, col="black", stat="identity", position="stack") + 
  theme_classic() +
  scale_fill_manual(values = c(reds,blues)) + 
  ylab("Occurance") + 
  xlab("Insertion or deletion size")
ggsave(gnorm,file="nonlog_edit_sizes_r50high.pdf",width=15,height=6)

