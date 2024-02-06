setwd("/Users/faylward.AYLWARD-9BCZN53/Documents/research/chord_diagram/")
library(circlize)
library(pheatmap)
library(usedist)

x <- read.table(file="revised_orthogroup_matrix.txt", header=T, row.names=1)
x[x>1] <- 1
#x[x==0] <- NA
ogs <- read.table(file="ogs_to_include.txt", header=F)
ogn <- ogs$V1
subset <- x[ogn,]
ham <- function(vect1, vect2) {
  s <- sum(vect1 == vect2 & vect1==1)
  return(s)
}
d <- as.matrix(dist_make(subset, ham))
pols <- c("OG0000003", "OG0000015", "OG0000122", "OG0000045", "OG0000020")
mcps <- c("OG0000011", "OG0000012", "OG0000025", "OG0000047", "OG0000049", "OG0000073", "OG0000001_C", "OG0000084", "OG0000057")

d2 <- d[pols, mcps]
d2
chordDiagram(d2, transparency = 0.2)


