setwd("/Users/faylward.AYLWARD-9BCZN53/Documents/genome_plots/")
library(genoPlotR)

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# workflow for genome sets
# predict genes in genomes, merge GFF files, make final dna_seqs file (last part needs python script)
# decide order, all-vs-all blastp, python script to sort output into big input df
# so both input dna_seqs and comparison files are large dfs that are split in R

library(genoPlotR)
x <- read.table(file="ostreo.txt", header=T)
out <- split(x , f = x$set)
out <- out[c('NC_028092.1', 'NC_028091.1', 'JF974316.1', 'JN225873.1')]

for(i in 1:length(out)) {
  new <- dna_seg(out[[i]][2:12])
  out[[i]] <- new
}
plot_gene_map(out)
tree <- newick2phylog("(((NC_028094:4.2,NC_024697:3.9):3.1,HQ634144:7.3):1);")
c1 <- as.comparison(read.table(file="OlV1__OlV2.ortho", header=T))
c2 <- as.comparison(read.table(file="OlV2__OlV4.ortho", header=T))
c3 <- as.comparison(read.table(file="OlV4__OtV6.ortho", header=T))
compare <- list(c1, c2, c3)

for(i in 1:length(compare)) {
  compare[[i]]$col <- makeTransparent('dodgerblue', 70)
}

#pdf("ostreo.pdf", width=14, height=4)
plot_gene_map(out, comparison=compare, dna_seg_labels=c('OlV1', 'OlV2', 'OlV4', 'OtV6'), xlims=list(c(-Inf, Inf), c(-Inf, Inf), c(-Inf, 5000, 5000, Inf), c(-Inf, Inf)))
#dev.off()


