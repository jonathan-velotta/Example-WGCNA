library(edgeR)
library(WGCNA)
library(Hmisc)
library(reshape)

counts <- read.csv("epas1_cat.counts.csv", stringsAsFactors=FALSE, row.names=1)#read in raw read count table
adrenal <- counts[, c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 34, 36, 37)] #dataset for adrenal tissue only
##################################EdgeR is used to create log transformed and normalized counts for WGCNA
#adrenal1 <- adrenal[, -c(5,8,17)]
adrenal1 <- adrenal
id <- read.csv("sample_id.csv")#dataset containing individual level information
adrenal.id <- subset(id, tissue == "adrenal")
#adrenal.id1 <- adrenal.id[-c(5,8,17),]
adrenal.id1 <- adrenal.id
match(adrenal.id1$mouse, colnames(adrenal1)) #verify that colnames and sample id names are in the same order
genotype <- adrenal.id1$geno
type <- adrenal.id1$geno2
table <- data.frame(genotype, type)
group <- factor(paste(table$genotype, table$type, sep="_"))
cbind(table, group=group)
table$genotype = as.factor(table$genotype)
table$type = as.factor(table$type)
table$genotype <- relevel(table$genotype, ref="LL")
table$type <- relevel(table$type, ref="L")
design <- model.matrix(~genotype, data=table)
design2 <- model.matrix(~type, data=table)
#rowMeans takes means of each row
adrenal1$mean = rowMeans(adrenal1)
head(adrenal1)
#filter by means
dim(adrenal1)
keep_adrenal1 = subset(adrenal1, mean >= 5)
dim(keep_adrenal1)
#clean up dataset
keep_adrenal1$mean = NULL
head(keep_adrenal1)
dim(keep_adrenal1)
y <- DGEList(counts=keep_adrenal1, group=group) # make a DGE list
dim(y)
y <- calcNormFactors(y) # normalize
#exploration plots
plotMDS(y)
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

################################################################# Begin WGCNA code
adrenal.norm <- cpm(y, log=TRUE, prior.count=2, normalized.lib.sizes=TRUE) #cpm normalized and log transformed expression dataset
#write.csv(adrenal.norm, file="Epas1.adrenal.norm.csv")
head(adrenal.norm)
adrenal.norm1 <- adrenal.norm
#transpose expression data for further analysis
Expr0 = as.data.frame(t(adrenal.norm1));
rownames(Expr0) = colnames(adrenal.norm1)
#check for genes and samples with too many missing values
gsg = goodSamplesGenes(Expr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(Expr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(Expr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  Expr0 = Expr0[gsg$goodSamples, gsg$goodGenes]
}
Expr = Expr0
#cluster the samples to check for outliers
sampleTree = hclust(dist(Expr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# load trait data
traitData <- read.csv("traitData_adrenal.csv")[, -2]

# Form a data frame analogous to expression data that will hold the clinical traits.
Samples = rownames(Expr);
traitRows = match(Samples, traitData$mouse)# make sure IDs are in the same order as Expr dataset
Traits0 = traitData[traitRows, -1];
Traits0 <- na.omit(Traits0)
rownames(Traits0) = na.omit(traitData[traitRows, 1]);
Traits = Traits0[, -c(2:5)]
Traits$geno <- NULL
collectGarbage();
match(rownames(Traits), rownames(Expr))# should be in order if datasets are matched

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(Expr, powerVector = powers, verbose = 5)
# Plot the results:
#pdf('wgcna/rv_beta_plot.pdf', h=4, w=7)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.93,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#dev.off()

# construct gene network
Net = blockwiseModules(Expr, power = 9, maxBlockSize = 15000,
                       TOMType = "signed", networkType = "signed",
                       minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ExprTOM",
                       verbose = 3)
table(Net$colors)
moduleLabels = Net$colors
moduleColors = labels2colors(Net$colors)
MEs = Net$MEs;
geneTree = Net$dendrograms[[1]];
table(moduleColors)
dim(table(moduleColors))

#######################################################################################################################
# Trait Module Associations
# Define numbers of genes and samples
nGenes = ncol(Expr);
nSamples = nrow(Expr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(Expr, moduleColors)$eigengenes
MEs = orderMEs(MEs0) #the rownames of this dataset are equal to Expr
MEs1 <- MEs[-c(1, 10, 13),]
match(rownames(Traits), rownames(Expr))
moduleTraitCor = cor(MEs1, Traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
corr.table <- as.data.frame(cbind(moduleTraitCor, moduleTraitPvalue))
colnames(corr.table) <- c("cor", "p")
corr.table$p.adjust <- p.adjust(corr.table$p, method="fdr")
corr.table <- corr.table[order(corr.table$p),]

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(Expr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
# Create a dataset containing all gene-specific information
genes=names(Expr)
geneInfo0 = data.frame(Gene = genes,
                       moduleColor = moduleColors)
# Add module membership information in the chosen order
oldNames = names(geneInfo0)
geneInfo0 = data.frame(geneInfo0, geneModuleMembership,
                       MMPvalue);
names(geneInfo0) = c(oldNames, paste("MM.", modNames, sep=""),
                     paste("p.MM.", modNames, sep=""))
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor);
geneInfo = geneInfo0[geneOrder, ]

## ANOVA on ME
## code from Plachetzki et al 2014
head(adrenal.id1)
match(rownames(Expr), adrenal.id1$mouse)
setSamples = rownames(Expr);

ME = MEs
modNames = substring(names(ME), 3)
anovas<-list()
anovasP<-list()
anovasF<-list()
for (i in 1:ncol(ME)){
  rME<-rank(ME[,i])
  anova<-summary(aov(rME ~ adrenal.id1$geno))
  anovasP[[i]]<-c(anova[[1]][["Pr(>F)"]][1]); names(anovasP)[i]<-modNames[i]
  anovasF[[i]]<-c(anova[[1]][["F value"]][1]); names(anovasF)[i]<-modNames[i]
  anovas[i]<-anova; names(anovas)[i]<-modNames[i]
}
P_vals<-data.frame(matrix(round(unlist(anovasP),4),dim(table(moduleColors))[1] ,1,byrow=T, dimnames=list(names(MEs), "P-Value")))
P <- p.adjust(P_vals$P.Value, method = "fdr")
P_vals$p.cor <- data.frame(P)
P_vals <- P_vals[order(P_vals$P.Value),]