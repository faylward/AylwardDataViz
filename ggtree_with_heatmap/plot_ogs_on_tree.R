
#library(ape)
#library(castor)
library(ggtree)
#library(tidyverse)
setwd("/home/aylward/PLV_pangenomics/ggtree/")

og <- t(read.table(file="revised_orthogroup_matrix.txt", sep="\t", header=T, row.names=1, check.names=F))

ogs <- c("OG0000001_C", "OG0000000_B", "OG0000025", "OG0000047", "OG0000049", "OG0000003", "OG0000015",  "OG0000045")
samp <- data.frame(og[,ogs])
samp[samp >0] <- "Yes"
samp[samp == 0] <- "No"

tree <- read.tree(file="rooted.nwk")

p <- ggtree(tree, size=1.3) + theme(plot.margin = unit(c(14,8,34,8), "mm")) +
  geom_treescale(offset=2)
 gheatmap(p, samp, offset=0, width=0.5, font.size=3, colnames_angle=40, hjust=1) + scale_x_ggtree() + 
  scale_fill_manual(breaks=c("Yes", "No"), values=c("black", "#FFFFFF"), name="OG occurrence") + scale_y_continuous(limits = c(-5, NA))

ggsave(file="tree_heatmap.png", width=12, height=10, units="in", dpi=600)
