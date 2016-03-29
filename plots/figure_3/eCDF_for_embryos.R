library(ggplot2)
library(ggthemes)
library(reshape2)

# developmental time course simulations
dev_colors = c("#515251","#0086FF","#0AA600","#F26E12")

# # our data
# tbl <- read.delim("~/google_drive/UW/shendure/CRISPR/2016_03_22_embryo_eCDF_code/test_simulation_100K.txt")
# 
# # plot by concentation and timepoint
# gg = ggplot(tbl[tbl$index < 30000,]) + 
#   geom_line(aes(index,mean,group=sample,col=hours),size=2) + 
#   facet_grid(. ~ concentration) + 
#   scale_color_colorblind() + 
#   theme_tufte() + xlab("Cells (sampled)") + ylab("average unique HMIDs")
# ggsave(gg,file="~/Desktop/average_unique_hmids_sampled.png",width=9,height=5)
# 
# raw_counts_rev = read.delim("/mount/vol10/CRISPR.lineage/nobackup/2016_03_05_Embryos/embryos_all_reads_Mar_23_2016_reverse.txt")
# 
# gg2 = ggplot(raw_counts_rev) + 
#   geom_line(aes(running_count,array,group=sample,col=stage),size=2,) + 
#   scale_color_colorblind() + 
#   facet_grid(. ~ rate) + 
#   theme_tufte() + ylab("unique HMIDs") + xlab("Cells")
# ggsave(gg2,file="~/Desktop/average_unique_hmids_total.png",width=9,height=5)

# our data
tbl2 <- read.delim("~/google_drive/UW/shendure/CRISPR/2016_03_22_embryo_eCDF_code/test_simulation_no_replace.txt")

tbl2$hours = factor(as.character(tbl2$hours), levels=c("72H","30H","9H","4.3H"), ordered=T)

one_third = subset(tbl2, concentration == "1/3X")
one_x = subset(tbl2, concentration == "1X")

one_third_full = ggplot(one_third) + 
  geom_line(aes(index,mean,group=sample,col=hours),size=3) + 
  geom_line(data = subset(one_third, hours == "4.3H"),aes(x = index, y = mean, group = sample), col = "#F26E12",size=1.5) + 
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


one_third_full_subset = ggplot(one_third) + 
  geom_line(aes(index,mean,group=sample,col=hours),size=1.5) + 
  geom_line(data = subset(one_third, hours == "4.3H"),aes(x = index, y = mean, group = sample), col = "#F26E12",size=1.5) + 
  scale_color_manual(values=dev_colors) + 
  theme_tufte() + xlab("") + ylab("") +
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) +
  theme(strip.text=element_text(family = "sans", size=12)) +
  xlim(c(0,1500)) + 
  ylim(c(0,750)) +
  theme(legend.position="none")

vp <- viewport(width = 0.3, height = 0.4, x = .75,
               y = .3)

full <- function() {
  print(one_third_full)
  theme_set(theme_bw(base_size = 8))
  theme_classic()
  print(one_third_full_subset, vp = vp)
  theme_set(theme_bw())
}
full()

one_full = ggplot(one_x) + 
  geom_line(aes(index,mean,group=sample,col=hours),size=3) + 
  geom_line(data = subset(one_x, hours == "4.3H"),aes(x = index, y = mean, group = sample), col = "#F26E12",size=1.5) + 
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

one_full_subset = ggplot(one_x) + 
  geom_line(aes(index,mean,group=sample,col=hours),size=1.5) + 
  geom_line(data = subset(one_x, hours == "4.3H"),aes(x = index, y = mean, group = sample), col = "#F26E12",size=1.5) + 
  scale_color_manual(values=dev_colors) + 
  theme_tufte() + xlab("") + ylab("") +
  theme(legend.text=element_text(family = "sans", size=12)) + 
  theme(legend.title=element_text(family = "sans", size=12)) + 
  theme(axis.text=element_text(family = "sans", size=12)) + 
  theme(axis.title=element_text(family = "sans", size=12)) +
  theme(strip.text=element_text(family = "sans", size=12)) +
  xlim(c(0,1500)) + 
  ylim(c(0,750)) +
  theme(legend.position="none")

vp <- viewport(width = 0.3, height = 0.4, x = .75,
               y = .3)

full <- function() {
  print(one_full)
  theme_set(theme_bw(base_size = 8))
  theme_classic()
  print(one_full_subset, vp = vp)
  theme_set(theme_bw())
}
full()


# ggsave(gg,file="~/Desktop/average_unique_hmids_no_resampling.png",width=6,height=3)
