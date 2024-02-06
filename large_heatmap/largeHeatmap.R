
setwd("/Users/faylward.AYLWARD-9BCZN53/Documents/research/vogclust/large_heatmap/")
library(pheatmap)
library(RColorBrewer)

fam <- read.table(file="Phage_families.tsv", sep="\t", row.names=1, header=T)
abund <- names(head(sort(table(fam$Family), decreasing=T), n=15))
annotrow <- data.frame(fam[which(fam$Family %in% abund), 3])
row.names(annotrow) <- row.names(fam)[which(fam$Family %in% abund)]
colnames(annotrow) <- "fam"

x <- read.table(file="cp2_inphared_05p_binary_profile_min5hit.tsv", header=T, row.names=1, sep="\t")
subset <- x[,colSums(x)>50]
pheatmap(subset, show_rownames = F, show_colnames = F, cluster_rows=T, color=c("white", "dodgerblue"), annotation_row=annotrow)

