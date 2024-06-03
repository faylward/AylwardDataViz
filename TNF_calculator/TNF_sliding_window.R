
setwd("/Users/faylward.AYLWARD-9BCZN53/Documents/research/tnf/TNF_calculator_sliding_window/")
library(Biostrings)
library(seqinr)
library(zoo)

# You need two files to run this script. In the example below, one file (ostreo_pseudocontigs.fasta) 
# is a reference of contigs/chromosomes that other contigs/chromosomes will be compared against. 
# The reference FASTA file should have only one record, but you can combine multiple records together 
# (i.e. pseudocontigs) if you want to include multiple contigs or chromosomes. 
# The second file (i.e. the target file, ostreo.fasta) can contain several records. For each record, a plot of the sliding 
# of the TNF correlation compared to the reference file is generated. 

# The default sliding window is 5000 bp - you may wish to alter this depending on your use case. 

# make function to calculate correlation - this will be used later
compare <- function(string, tnf) {
  dna <- DNAString(paste(string, collapse=""))
  newtnf <- scale(oligonucleotideFrequency(dna, 4), center=F)
  corr <- abs(cor(newtnf, tnf))
  return(corr)
}

# load in and process the reference file
ref <- read.fasta(file="ostreo_pseudocontigs.fasta", as.string=T, set.attributes=F)
string <- ref[[1]]
st <- strsplit(string, "")
item <- st[[1]]
tnf <- scale(oligonucleotideFrequency(DNAString(string),4))

# load in and process the target file
record <- read.fasta(file="ostreo.fasta", as.string=T, set.attributes=F)
for(i in 1:length(record)) {
  string <- record[[i]]
  name <- names(record)[i]
  st <- strsplit(string, "")
  item <- st[[1]]
  
  info <- rollapply(item, 5000, by=5000, compare, tnf=tnf) - 1
  xaxis <- seq(from=4999, length.out=length(info), by=5000)
 
  jpeg(paste(name, ".jpg", sep="_"), width=7, height=5, units="in", quality=100, res=200)
  plot(info, x=xaxis, type="l", col="#2c7bb6", lwd=3, ylim=c(-1, 1), main=name)
  abline(h=0)
  dev.off()
}
