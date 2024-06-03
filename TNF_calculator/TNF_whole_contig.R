
setwd("/Users/faylward.AYLWARD-9BCZN53/Documents/research/tnf/TNF_calculator_sliding_window/")
library(Biostrings)
library(seqinr)
library(zoo)

ref <- read.fasta(file="ostreo_pseudocontigs.fasta", as.string=T, set.attributes=F)
string <- ref[[1]]
st <- strsplit(string, "")
item <- st[[1]]
tnf <- scale(oligonucleotideFrequency(DNAString(string),4))

record <- read.fasta(file="ostreo.fasta", as.string=T, set.attributes=F)
name_list <- list()
cor_list <- list()
for(i in 1:length(record)) {
  string <- record[[i]]
  name <- names(record)[i]
  st <- strsplit(string, "")
  item <- st[[1]]
  
  c <- compare(item, tnf)
  cor_list[i] <- c
  name_list[i] <- name
}
name_list <- as.character(name_list)
cor_list <- as.character(cor_list)
t <- cbind(name_list, cor_list)
write.table(t, file="TNF.correlations.txt", quote=F)

