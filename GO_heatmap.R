### This script will filter DEGs annotated with a GO-term of interest and plot as a heatmap ###

library(dplyr)
library(tidyr)
library(gplots)
library(tidyverse)
library(pheatmap)

### filtering genes by GO term ###
iso2gene = read.csv("iso2gene.csv", sep=",", header = TRUE)
BP_Starved_GO=read.table("BP_Starved_GO.csv",sep="\t")
colnames(BP_Starved_GO)=c("name","term","lev","Iso","value")
filtered=filter(BP_Starved_GO,grepl('GO:0006955',term)) #change to GO term of interest
filtered=filtered[!duplicated(filtered$Iso), ] %>%
  dplyr::select("Iso") %>%
  left_join(iso2gene)

### get pvals and add gene names ###
rldpvals=read.csv("Starved_RLDandPVALS.csv")
head(rldpvals)
colnames(rldpvals)=c("Iso","10F","10S","1S","1F","2S","2F","4S","4F","6F","6S","7F","7S","8F","8S","9F","9S","pval.Starved","padj.Starved")
head(rldpvals)
rldpvals = rldpvals %>%
  left_join(filtered) %>%
  drop_na() %>%
  mutate(gene_unique = make.names(gene, unique = TRUE, )) %>%
  dplyr::select(-Iso, -gene) %>%
  column_to_rownames(var = "gene_unique")

topnum= 100 # number of DEGS
top100=head(rldpvals[order(rldpvals$padj.Starved), ],topnum)
head(top100)
length(top100[,1])
summary(top100)

### plotting expression ###
p.val=0.1 # FDR cutoff
conds=top100[top100$padj.Starved<=p.val & !is.na(top100$padj.Starved),]
length(conds[,1]) #31
exp=conds[,1:16]
means=apply(exp,1,mean) # means of rows
explc=exp-means # subtracting them
head(explc)
col0=colorRampPalette(rev(c("#ea6227","#f09167","white","#7caeff","#0666ff")))(100)
pheatmap(explc,cluster_cols=T,scale="row",color=col0, show_rownames = T)
