# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

library(DESeq2) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library("ggplot2")
library(dplyr)
library(RColorBrewer)
library(gplots)
library("pheatmap")
library(vegan)
library(ggrepel)
library(tidyverse)

###conduct array quality metrics to detect and remove outliers
#read in counts 
countData <- read.table("nve_counts.txt")
head(countData)
#22857
#read in counts 
countData <- read.table("nve_counts.txt")
head(countData)
length(countData[,1])
#22857
names(countData)=c("10AF","10OS","1AS","1OF","2AS","2OF","4AS","4OF","6AF","6OS","7AF","70S","8AF","8OS","9AF","9OS")
head(countData)
treat=c("F","S","S","F","S","F","S","F","F","S","F","S","F","S","F","S")
genet=c("10","10","1","1","2","2","4","4","6","6","7","7","8","8","9","9")
g=data.frame(treat, genet)
g$genet=as.factor(g$genet)
str(g)
g
colData= g

dds=DESeqDataSetFromMatrix(countData=countData,
                           colData = g,
                           #design = ~treat
                           design = ~treat+genet
                          )
vsd.ge=assay(vst(dds))
rl=vst(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,outdir=v,intgroup=c("treat"),force=T)
#dev.off() only for windows bug

###No Outliers detected###

countData <- read.table("nve_counts.txt")
head(countData)
length(countData[,1])
#22857

names(countData)=c("10AF","10OS","1AS","1OF","2AS","2OF","4AS","4OF","6AF","6OS","7AF","70S","8AF","8OS","9AF","9OS")
head(countData)
totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, ylab="raw counts")
#10AF    10OS     1AS     1OF     2AS     2OF     4AS     4OF     6AF     6OS 
#1410416 1601085 1224858 1828381 1152596 1133896 1494445 1982020 1194611 1243135 
#7AF     70S     8AF     8OS     9AF     9OS 
#1254392 1612949 1204938 1784150 1670875 1456236 
min(totalCounts) #1133896
max(totalCounts)  #1982020

treat=c("F","S","S","F","S","F","S","F","F","S","F","S","F","S","F","S")
genet=c("10","10","1","1","2","2","4","4","6","6","7","7","8","8","9","9")
g=data.frame(treat, genet)
g$genet=as.factor(g$genet)
str(g)
g
colData<- g
dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~treat+genet)

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
head(dds)
res<- results(dds)

###### rlogged dds for GO down the line #####
starved_rlogged = DESeq2::rlog(dds, blind = TRUE)
write.csv(assay(starved_rlogged), "starved_rlogged.csv")
#######

#Look at dispersions plot
plotDispEsts(dds, main="Dispersion plot Snails")

##### S vs F pairwise comparisons #####
colData$Starved<-factor(colData$treat, levels=c("S","F"))
resStarved <- results(dds, contrast=c("treat","S","F"))
table(resStarved$padj<0.1)
#FALSE  TRUE 
#7945   288
summary(resStarved) #uses padj < 0.1
#593 (2.6%) genes downregulated in S relative to F
#118 (0.52%) genes upregulated
nrow(resStarved[resStarved$padj<0.1 & !is.na(resStarved$padj),])  # Num significantly deferentially expressed genes excluding the no/low count genes   #711
#711 significant DEGs
#dev.off() only for windows bug
plotMA(resStarved, main="Fed vs Starved", ylim=c(-4,3))
results <- as.data.frame(resStarved)
head(results)
write.table(resStarved, file="Starved.txt", quote=F, sep="\t")

cd <- read.table("Starved.txt")
head(cd)

##### make the GO table for MWU #####
go_input_Starved = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)
head(go_input_Starved)
colnames(go_input_Starved) <- c("gene", "pval")
head(go_input_Starved)
write.csv(go_input_Starved, file="Starved_GO.csv", quote=F, row.names=FALSE)

##### get pvals #####
valStarved=cbind(resStarved$pvalue, resStarved$padj)
head(valStarved)
colnames(valStarved)=c("pval.Starved", "padj.Starved")
length(valStarved[,1]) #22857
table(complete.cases(valStarved)) 
#FALSE TRUE
#14624 8233

##### making rlogdata and pvals table #####
rlog=rlogTransformation(dds, blind=TRUE) #takes a bit
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])
#22857
rldpvals=cbind(rld,valStarved)
head(rldpvals)
dim(rldpvals)
# [1] 22857    18
table(complete.cases(rldpvals))
#FALSE  TRUE 
#14624  8233 
write.csv(rldpvals, "Starved_RLDandPVALS.csv", quote=F)

#### heat map of sample distances for pca #####
rldpvals <- read.csv(file="Starved_RLDandPVALS.csv", row.names=1)
head(rldpvals)
rld=rldpvals[,1:16]
head(rld)
sampleDists <- dist(t(rld))
sampleDistMatrix <- as.matrix( sampleDists )
names=c("10F","10S","1S","1F","2S","2F","4S","4F","6F","6S","7F","7S","8F","8S","9F","9S")
colnames(sampleDistMatrix)=paste(names) # treat
rownames(sampleDistMatrix)=paste(names) # treat
heat.colors = colorRampPalette(rev(c("blue","yellow","red")),bias=0.3)(100)
pheatmap(sampleDistMatrix,color = heat.colors,cex=0.9,border_color=NA,cluster_rows=T,cluster_cols=T)

rld_t=t(rld)
pca <- prcomp(rld_t,center = TRUE, scale. = TRUE)
head(pca)
li <- pca$sdev^2 / sum(pca$sdev^2)
pc1v <- round(li[1] * 100, 1)
pc2v <- round(li[2] * 100, 1)
pca_s <- as.data.frame(pca$x)
head(pca_s)
pca_s <- pca_s[,c(1,2)]
pca_s$Samples = c("10F","10S","1S","1F","2S","2F","4S","4F","6F","6S","7F","7S","8F","8S","9F","9S")
pca_s$treat=colData$treat
pca_s$genet=c("10","10","1","1","2","2","4","4","6","6","7","7","8","8","9","9")
head(pca_s)

cbPalette <- c("darkorchid3","gray0")
ggplot(pca_s, aes(PC1, PC2, color = treat, pch = genet)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(15,19,17,18,0,1,2,5))+
  geom_text_repel(aes(label=Samples)) + #adds labels
  stat_ellipse(type="t",aes(group = treat)) + #adds ellipse based of treatment
  scale_colour_manual(values=cbPalette)+
  theme_bw() +
  xlab(paste0("PC1: ",pc1v,"% variance")) +
  ylab(paste0("PC2: ",pc2v,"% variance")) 

head(pca)
adonis(pca$x ~ treat+genet, data = pca_s, method='eu', na.rm = TRUE)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#treat      1     29208   29208  1.6246 0.08519  0.005 ** 
#genet      7    187802   26829  1.4923 0.54776  0.001 ***
#Residuals  7    125845   17978         0.36705           
#Total     15    342855                 1.00000           
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

