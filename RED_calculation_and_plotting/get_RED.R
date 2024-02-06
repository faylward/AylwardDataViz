setwd("/home/aylward/RED_for_Jonatas/")
library(ape)
library(castor)
library(ggplot2)

#tree <- read.tree("tree.norogues.nwk")
#meta <- read.table(file="taxa_annotation.txt", header=T, sep="\t", comment.char = "")

tree <- read.tree("tree.norogues.nwk")
meta <- read.table(file="taxa_annotation.txt", header=T, sep="\t", comment.char = "")

##################################################################################################################
################################################### Family level ###################################################
##################################################################################################################

names(table(meta$Family))
clades <-  c("AF_01", "AF_02", "AG_01", "AG_02", "AG_03", "AG_04", "PV_02", "PV_03", "PV_04", "PV_05", "IM_01", "IM_02", "IM_06", "IM_07", "IM_08", "IM_09", "IM_12", "IM_13", "IM_14", "IM_16", "IM_18", "PM_01", "PM_02", "PM_03", "PM_04", "PM_05", "PM_06", "PM_07", "PM_08", "PM_09", "PX_01")
red_list <- list()
reds <- get_reds(tree)
for(i in 1:length(clades)){
  genomes <- as.character(meta$genome[which(meta$Family == clades[i])])
  lca <- getMRCA(tree, tip=which(tree$tip.label %in% genomes)) - 1382
  red <- reds[lca]
  red_list[i] <- red
  print(paste(clades[i], " ", lca, " ", red))
}

redvals <- as.numeric(red_list)
mat <- data.frame(cbind(clades, rep("family", length(redvals)), redvals))
mat[order(mat$redvals),]
mat$redvals <- as.numeric(as.character(mat$redvals))
color_dict <-  c("PX_01" = "PX", "AF_01" = "AF", "AF_02" = "AF", "AF_03" = "AF", "PM_01" = "PM", "PM_02"= "PM", "PM_03"= "PM", "PM_04"= "PM", "PM_05"= "PM", "PM_06"= "PM", "PM_07"= "PM", "PM_08"= "PM", "PM_09"= "PM", "AG2_02"= "AG2", "AG2_03"= "AG2", "AG2_04"= "AG2", "AG2_05"= "AG2", "PV_06"= "PV", "AG1_01"= "AG1", "AG1_02"= "AG1", "AG1_03"= "AG1", "AG1_04"= "AG1", "IM_01"= "IM", "IM_02"= "IM", "IM_06"= "IM", "IM_07"= "IM", "IM_08"= "IM", "IM_09"= "IM", "IM_12"= "IM", "IM_13"= "IM", "IM_14"= "IM", "IM_16"= "IM", "IM_18"= "IM")
mat2 <- cbind(mat, color_dict[as.character(mat$clades)]); colnames(mat2) <- c('clades', 'group', 'redvals', 'color')

clades <-  c("PX_01", "AF_01", "AF_02", "PM_01", "PM_02", "PM_03", "PM_04", "PM_05", "PM_06", "PM_07", "PM_08", "PM_09", "PV_02", "PV_03", "PV_04", "PV_05", "PV_06", "AG_01", "AG_02", "AG_03", "AG_04", "IM_01", "IM_02", "IM_06", "IM_07", "IM_08", "IM_09", "IM_12", "IM_13", "IM_14", "IM_16", "IM_18")
pdf("family_RED.pdf", height=7, width=6)
ggplot(mat2, aes(y=redvals, x=clades)) + geom_bar(stat="identity") + coord_flip() +  scale_y_continuous(limits=c(0, 0.8)) + scale_x_discrete(limits=rev(clades)) + theme_bw()
dev.off()

family <- mat

color <- rep("black", length(tree$tip.label))
color[which(tree$tip.label %in% genomes)] <- "red"
plot(tree, edge.color = color)

meta[which(meta$genome  %in% genomes), ]

##################################################################################################################
################################################### Orderlevel ###################################################
##################################################################################################################

clades <- c("Chitovirales", "Asfuvirales", "Pimascovirales", "Pandoravirales", "Algavirales",  "Imitervirales")

red_list <- list()
reds <- get_reds(tree)
for(i in 1:length(clades)){
  genomes <- as.character(meta$genome[which(meta$Order == clades[i])])
  lca <- getMRCA(tree, tip=which(tree$tip.label %in% genomes)) - 1382
  red <- reds[lca]
  red_list[i] <- red
  print(paste(clades[i], " ", lca, " ", red))
}

redvals <- as.numeric(red_list)
mat <- data.frame(cbind(clades, rep("order", length(redvals)), redvals))
mat$redvals <- as.numeric(as.character(mat$redvals))
color <- c('firebrick', 'darkgoldenrod', 'dodgerblue', 'green4', 'lightgreen', 'purple')
mat2 <- cbind(mat, color)

#mat[order(mat$redvals),]
pdf("order_RED.pdf", height=1.8, width=6)
ggplot(mat2, aes(y=redvals, x=clades)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + scale_y_continuous(limits=c(0, 0.8)) + scale_x_discrete(limits=rev(clades))
dev.off()

order <- mat

##################################################################################################################
################################################### Class level ###################################################
##################################################################################################################

### ClassLevel
clades <- c("Megaviricetes", "Pokkesviricetes")

red_list <- list()
reds <- get_reds(tree)
for(i in 1:length(clades)){
  genomes <- as.character(meta$genome[which(meta$Class == clades[i])])
  lca <- getMRCA(tree, tip=which(tree$tip.label %in% genomes)) - 1382
  red <- reds[lca]
  red_list[i] <- red
  print(paste(clades[i], " ", lca, " ", red))
}

redvals <- as.numeric(red_list)
mat <- data.frame(cbind(clades, rep("class", length(redvals)), redvals))
mat[order(mat$redvals),]
mat$redvals <- as.numeric(as.character(mat$redvals))

pdf("class_RED.pdf", height=1.1, width=6)
ggplot(mat, aes(y=redvals, x=clades)) + geom_bar(stat="identity") + coord_flip() + theme_bw() + scale_y_continuous(limits=c(0, 0.8))
dev.off()

class <- mat

##################################################################################################################
################################################### Genus level? ###################################################
##################################################################################################################
library(gplots)
library(WGCNA)

##################### produce genera level classifications ########################################
#tree <- read.tree("news_tree.nwk")
umtree <- chronos(tree)
write.tree(umtree, file="umtree.nwk")
umtree <- read.tree("umtree.nwk")

hclust <- as.hclust(umtree)
c <- cutree(hclust, h=0.13)

for(i in 1:length(c)) {
  c[i] <- paste("g", c[i], sep="")
}


color <- labels2colors(c)
hex <- col2hex(color)
t <- cbind(names(c), hex)
write.table(t, file="color_lookup.txt", sep="\t", quote=F, row.names=F)

meta <- read.table(file="taxa_annotation.txt", header=T, sep="\t", row.names=1)
genomes <- row.names(meta)
meta2 <- cbind(meta[genomes,], c[genomes])
colnames(meta2)[7] <- "genus"

write.table(meta2, file="genus_level_classification.txt", sep="\t", quote=F)

##############################################################################################

color <- labels2colors(meta$Genus)
hex <- col2hex(color)
t <- cbind(as.character(meta$genome), hex)
write.table(t, file="color_lookup.txt", sep="\t", quote=F, row.names=F)

clades <- names(table(c)[table(c) > 1])

clades <- names(table(meta$Genus)[table(meta$Genus) > 1])
red_list <- list()
reds <- get_reds(tree)
for(i in 1:length(clades)){
  genomes <- as.character(meta$genome[which(meta$Genus %in% clades[i])])
  lca <- getMRCA(tree, tip=which(tree$tip.label %in% genomes)) - 1382
  red <- reds[lca]
  red_list[i] <- red
  print(paste(clades[i], " ", lca, " ", red))
}

redvals <- as.numeric(red_list)
mat <- data.frame(cbind(clades, rep("genus", length(redvals)), redvals))
mat[order(mat$redvals),]
mat$redvals <- as.numeric(as.character(mat$redvals))

genus <- mat


complete <- rbind(genus, family, order, class)
colnames(complete) <- c('clade', 'level', 'RED')
complete$RED <- as.numeric(complete$RED)

pdf(file="RED.pdf", height=2, width=6)
ggplot(complete, aes(x=RED, y=level)) + geom_jitter(colour="dodgerblue4", alpha=0.6, height=0.2, width=0) + theme_bw()
dev.off()




meta <- read.table(file="taxa_annotation.txt", header=T, sep="\t")
color <- labels2colors(meta$subclade_new_notation)
hex <- col2hex(color)
t <- cbind(as.character(meta$genome), hex)
write.table(t, file="family_level_lookup.txt", sep="\t", quote=F, row.names=F)

##################################################################################################################
#################################################### Treemap #####################################################
##################################################################################################################
library(treemap)
library(gplots)
meta <- read.table(file="taxa_annotation.txt", header=T, sep="\t")
write.table(table(meta$subclade_new_notation), file="breakdown.txt", sep="\t", quote=F)
write.table(table(meta$genus), file="breakdown.txt", sep="\t", quote=F)


x <- read.table("order_family_treeplot.txt", sep="\t", header=T)
color <- c('firebrick', 'darkgoldenrod', 'dodgerblue', 'green4', 'lightgreen', 'purple')
colhex <- col2hex(color)  

pdf("treemap_family.pdf", height=4, width=6)
p <- treemap(x,
             index=c("Order","Family"),
             vSize="Value",
             type="index",
             #vColor = "color",
             #palette = "Set2",
             bg.labels=c("white"),
             align.labels=list(
               c("center", "center"), 
               c("right", "bottom"),
               #algorithm="pivotSize",
               palette=colhex
             )  
) 
dev.off()


gdp_json <- vt_export_json(
  
vt_input_from_df(x, hierachyVar0 = "h1", hierachyVar1 = "h2", colorVar = "color", weightVar="weight", labelVar = "codes")

##################################################################################################################
################################################### Alternative Root ###################################################
##################################################################################################################


tree <- read.tree("poxroot.nwk")
clades <-  c("AF_01",  "AF_03",  "EP_03",  "EP_04",  "EP_05",  "EP_06",  "IR_01",  "IR_02",  "IR_03",  "IR_04",  "LP_01",  "LP_02",  "LP_03",  "LP_04",  "MA_01",  "MA_02", "MC_01",  "MC_05",  "MC_07",  "MC_08",  "MC_09",  "MC_12",  "MC_13",  "MC_14",  "MC_16",  "MC_18",  "PT_01",  "PT_02",  "PT_03",  "PX_01")

red_list <- list()
reds <- get_reds(tree)
for(i in 1:length(clades)){
  genomes <- as.character(meta$genome[which(meta$subclade2 == clades[i])])
  lca <- getMRCA(tree, tip=which(tree$tip.label %in% genomes)) - 1380
  red <- reds[lca]
  red_list[i] <- red
  print(paste(clades[i], " ", lca, " ", red))
}

redvals_poxroot <- as.numeric(red_list)
mat2 <- data.frame(cbind(clades, redvals_poxroot))


#### merge the two
diff <- as.numeric(as.character(mat$redvals)) - as.numeric(as.character(mat2$redvals_poxroot))
comb <- cbind(mat, mat2, diff)
diff <- as.numeric(as.character(comb$redvals)) - as.numeric(as.character(comb$redvals_poxroot))
diff









hist(redvals, xlim=c(0, 1), breaks = 10, col="dodgerblue")


i <- "MC_16"
genomes <- as.character(meta$genome[which(meta$subclade == i)])
meta[meta$genome %in% genomes,]
