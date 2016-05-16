library(ggplot2)
library(ggthemes)
library(reshape2)

# developmental time course colors 
# dev_colors = c("#F26E12","#0AA600","#0086FF","#515251")
# dev_colors = c("#515251","#0086FF","#0AA600","#F26E12")
dev_colors = c("#0086FF","#F26E12","#515251","#0AA600")
tbl2 <- read.delim("~/google_drive/UW/shendure/CRISPR/2016_03_22_embryo_eCDF_code/test_simulation_no_replace_version2.txt")

tbl2$hours = factor(as.character(tbl2$hours), levels=c("4.3H","9H","30H","72H"), ordered=T)

one_third = subset(tbl2, concentration == "1/3X")
one_x = subset(tbl2, concentration == "1X")

line_width = 3

one_third_plot = ggplot() +
  geom_step(data=subset(one_third,hours == "72H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  geom_step(data=subset(one_third,hours == "30H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  geom_step(data=subset(one_third,hours == "9H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  geom_step(data=subset(one_third,hours == "4.3H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  scale_color_manual(values=dev_colors) +
  theme_tufte() + xlab("Cells (sampled)") + ylab("Average unique HMIDs") +
  theme(legend.text=element_text(family = "sans", size=16)) +
  theme(legend.title=element_text(family = "sans", size=16)) +
  theme(axis.text=element_text(family = "sans", size=16)) +
  theme(axis.title=element_text(family = "sans", size=16)) +
  theme(strip.text=element_text(family = "sans", size=16)) +
  theme(legend.position="none") +
  xlim(c(0,35000)) +
  ylim(c(0,4500))

ggsave(one_third_plot,file="~/Desktop/one_third_simulated.pdf",width=8,height=8)
ggsave(one_third_plot,file="~/Desktop/one_third_simulated.png",width=8,height=8)
one_x_plot = ggplot() +
  geom_step(data=subset(one_x,hours == "72H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  geom_step(data=subset(one_x,hours == "30H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  geom_step(data=subset(one_x,hours == "9H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  geom_step(data=subset(one_x,hours == "4.3H"),aes(index,mean,col=hours,group=sample),direction="vh",size=line_width) + 
  scale_color_manual(values=dev_colors) +
  theme_tufte() + xlab("Cells (sampled)") + ylab("Average unique HMIDs") +
  theme(legend.text=element_text(family = "sans", size=16)) +
  theme(legend.title=element_text(family = "sans", size=16)) +
  theme(axis.text=element_text(family = "sans", size=16)) +
  theme(axis.title=element_text(family = "sans", size=16)) +
  theme(strip.text=element_text(family = "sans", size=16)) +
  theme(legend.position="none") +
  xlim(c(0,35000)) +
  ylim(c(0,4500))

ggsave(one_x_plot,file="~/Desktop/one_x_simulated.pdf",width=8,height=8)
ggsave(one_x_plot,file="~/Desktop/one_x_simulated.png",width=8,height=8)