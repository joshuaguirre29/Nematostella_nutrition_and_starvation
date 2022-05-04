#source("http://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
#biocLite("flashClust")

library(DESeq2)
library(WGCNA)
library(flashClust)
library(pheatmap)

### Get all genes into WGCNA but remove the genes with low basemean values ###
#navigate into the correct directory
countData <- read.table("nve_counts.txt")
head(countData)
length(countData[,1])
#22857
names(countData)
names(countData)=c( "F_10A", "S_10O", "S_1A", "F_1O", "S_2A", "F_2O", "S_4A", "F_4O",  "F_6A", "S_6O", "F_7A", "S_7O", "F_8A", "S_8O", "F_9A", "S_9O")
head(countData)

treat=c( "F", "S", "S", "F", "S", "F", "S", "F",  "F", "S", "F", "S", "F", "S", "F", "S")
g=data.frame(treat)
g
colData<- g

dds<-DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~ treat) 
dds<-DESeq(dds)

res<- results(dds)
#filter for contigs with average(baseMean) >3
res3<-res[res$baseMean>3, ]
dim(res) #22857 6

# get rlog data (better transformation when size factors vary across samples)
rld <- rlogTransformation(dds, blind=TRUE, fitType="local")
head(rld)
rld_wg=(assay(rld))
head(rld_wg)
nrow(rld_wg)
#22857
rldFiltered=(assay(rld))[(rownames((assay(rld))) %in% rownames(res3)),]
nrow(rldFiltered)
#15204
write.csv( rldFiltered,file="stella_wgcna_allgenes.csv",quote=F,row.names=T)

### Data input and cleaning ###
options(stringsAsFactors=FALSE)
allowWGCNAThreads()
dat=read.csv("stella_wgcna_allgenes.csv")
head(dat) 
rownames(dat)<-dat$X
head(dat)
dat$X=NULL
head(dat)
names(dat)
nrow(dat)
#15204
datExpr0 = as.data.frame(t(dat))

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #if TRUE, no outlier genes, if false run the script below
if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse=", ")))
  datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg=goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK 
dim(datExpr0) 
#12585  No outliers detected

### Outlier detection incorporated into trait measures. ###
traitData= read.csv("stella_wgcna_traits.csv", row.names=1)
dim(traitData)
head(traitData)
names(traitData)

### Form a data frame analogous to expression data that will hold the clinical traits. ###
dim(datExpr0)
rownames(datExpr0)
rownames(traitData)=rownames(datExpr0)
traitData$Sample= NULL
datTraits=traitData
table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
head(datTraits)
head(datExpr0)

#sample dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0),type="unsigned")
# this calculates the whole network connectivity
k=as.numeric(apply(A,2,sum))-1
# standardized connectivity
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
# Convert traits to a color representation where red indicates high values
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

save(datExpr0, datTraits, file="stella_Samples_Traits_ALL.RData")

##### Network construction and module detection - this section can take a lot of time you might consider running it on a cluster for a larger dataset ###
options(stringsAsFactors = FALSE)
#enableWGCNAThreads() use this in base R
allowWGCNAThreads() 
lnames = load(file="stella_Samples_Traits_ALL.RData")

# Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=2), seq(15,30, by=0.5)); #may need to adjust values for other analysis
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="unsigned", verbose = 2) #want smallest value, closest to 0.9 (but still under)

# Plot the results:
quartz()
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower=7 #smallest value to plateau at ~7
adjacency=adjacency(datExpr0, power=softPower,type="unsigned")
#translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "unsigned")
dissTOM= 1-TOM

geneTree= flashClust(as.dist(dissTOM), method="average")
quartz()
sizeGrWindow(10,6)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)

minModuleSize=90
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)

#dynamicMods
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17 
# 106 2426 1830 1492 1412 1317  912  637  591  464  442  429  422  417  385  350  350  346 
# 18   19   20 
# 307  290  279

dynamicColors= labels2colors(dynamicMods)
quartz()
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merg modules whose expression profiles are very similar
#calculate eigengenes
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 7)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_stella_nomerge.RData")

lnames = load(file = "Network_stella_nomerge.RData")
#plot
quartz()
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")

#here I am not collapsing any modules.
MEDissThres= 0
abline(h=MEDissThres, col="red")

merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

pdf(file="MergeNetwork.pdf", width=20, height=20)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
dev.off()

moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_unsigned_0.RData")

##### Relating modules to traits and finding important genes - This can be done locally #####
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "stella_Samples_Traits_ALL.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Network_unsigned_0.RData");
lnames = load(file = "Network_stella_nomerge.RData");
lnames

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)
#moduleColors
# moduleColors
# black         blue        brown         cyan        green  greenyellow         grey 
# 637         1830         1492          385         1317          429          106 
# grey60    lightcyan   lightgreen  lightyellow      magenta midnightblue         pink 
# 346          350          307          290          464          350          591 
# purple          red    royalblue       salmon          tan    turquoise       yellow 
# 442          912          279          417          422         2426         1412

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#represent module trait correlations as a heatmap
quartz()
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

### Gene relationship to trait and important modules:
weight = as.data.frame(datTraits$starved);
names(weight) = "starved"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")

### Gene-trait significance correlation plots
module = "green" #change to module of interest for others
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("ModMem in", module, "module"),
                   ylab = "Gene Sig for Starvation",
                   main = paste("MM vs. GS\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

### Making VSD files by module for GO plot functions
vs=t(datExpr0)
cands=names(datExpr0[moduleColors=="green"])

c.vsd=vs[rownames(vs) %in% cands,]
head(c.vsd)
nrow(c.vsd) #should correspond to module size, green: 1317
head(c.vsd)
write.csv(c.vsd,"rlog_MMgreen.csv",quote=F)

##### heatmap of module expression with bar plot of eigengene ####
sizeGrWindow(8,7);
which.module="green" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr0[,moduleColors==which.module ] #replace where says subgene below to plot all rather than just subset

quartz()
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")

### Gene relationship to trait and important modules: Gene Significance and Module membership
dat=read.csv("stella_wgcna_allgenes.csv")
rownames(dat)<-dat$X
dat$X=NULL
allkME =as.data.frame(signedKME(t(dat), MEs))
head(allkME)
vsd=read.csv(file="rlog_MMgreen.csv", row.names=1)
head(vsd)
gg=read.csv("iso2gene.csv", sep=",", header = TRUE)
head(gg)

### Plot
whichModule="green"
top=50

datME=MEs
vsd <- read.csv("stella_wgcna_allgenes.csv", row.names=1)
head(vsd)
datExpr=t(vsd)
modcol=paste("kME",whichModule,sep="")
head(vsd)
sorted=vsd[order(allkME[,modcol],decreasing=T),]
hubs=sorted[1:top,]

## attaching gene names, ignore to make heatmap with isogroup names instead of gene names ##
summary(hubs)

gnames=c();counts=0
for(i in 1:length(hubs[,1])) {
  if (row.names(hubs)[i] %in% gg$Iso) { 
    counts=counts+1
    gn=gg[gg$Iso==row.names(hubs)[i],2]
    if (gn %in% gnames) {
      gn=paste(gn,counts,sep=".")
    }
    gnames=append(gnames,gn)
  } else { 
    gnames=append(gnames,i)
  }
} 
row.names(hubs)=gnames
length(hubs)
##

contrasting = colorRampPalette(rev(c("#ea6227","#f09167","white","#7caeff","#0666ff")))(100)
pheatmap(hubs,scale="row",col=contrasting,border_color=NA, main=paste(whichModule,"top",top,"kME",sep=""))

###fisher for GO
vsd <- read.csv("stella_wgcna_allgenes.csv", row.names=1)
head(vsd)
options(stringsAsFactors=FALSE)
data=t(vsd)
allkME =as.data.frame(signedKME(data, MEs))

whichModule="green" # change module of interest

length(moduleColors) #15204
inModule=data.frame("module"=rep(0,nrow(vsd)))
row.names(inModule)=row.names(vsd)
head(inModule)
nrow(inModule)
genes=row.names(vsd)[moduleColors == whichModule]
table(moduleColors == whichModule)
head(genes)
length(genes)
inModule[genes,1]=1
sum(inModule[,1])
head(inModule)
write.csv(inModule,file=paste(whichModule,"_fisher.csv",sep=""),quote=F)

modColName=paste("kME",whichModule,sep="")
modkME=as.data.frame(allkME[,modColName])
row.names(modkME)=row.names(allkME)
names(modkME)=modColName
write.csv(modkME,file=paste(whichModule,"_kME.csv",sep=""),quote=F)

###-end-###
