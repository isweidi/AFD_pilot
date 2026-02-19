##WGCNA pipeline for somascan 20 samples##
setwd('/Users/WGCNA')

##############################################PCA##############################################
library(tidyverse)
expr_matrix <- read.csv("genename_20sample_pp.csv",header=T)
gene_ids <- expr_matrix[, 1] 
expr_matrix2 <- expr_matrix[, -1]
rownames(expr_matrix2) <- gene_ids
expr_matrix3 <- t(expr_matrix2)
colnames(expr_matrix3) <- gene_ids
constant_cols <- apply(expr_matrix3, 2, function(col) var(col) == 0)
expr_matrix3_filtered <- expr_matrix3[, !constant_cols]
dataExpr1 <- expr_matrix3_filtered


# Perform PCA
pca_result3 <- prcomp(expr_matrix3_filtered, scale. = TRUE)

# Create a data frame for PCA results
pca_data <- data.frame(pca_result3$x)
pca_data$Sample <- rownames(pca_data)

pdf("./image/pca_plot.pdf", width = 8, height = 6)

# Plot PCA with sample labels
ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point() +
  geom_text(vjust = 0.15, hjust = 0.15) +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2") +
  theme_minimal()
dev.off()

##############################################WGCNA analysis##########################################
library(WGCNA)
head(expr_matrix3_filtered)
expr2<-expr_matrix3_filtered

dataExpr1<-expr2
gsg=goodSamplesGenes(dataExpr1,verbose=3)
gsg$allOK
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr1 = dataExpr1[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(dataExpr1)
nGenes

nSamples = nrow(dataExpr1)
nSamples

traitData=read.table("TraitData.txt",header=T,sep="\t",row.names=1)
traitData

fpkmSamples<-rownames(dataExpr1)

traitSamples<-rownames(traitData)

traitRows<-match(fpkmSamples,traitSamples)

dataTraits<-traitData[traitRows,]

sampleTree2<-hclust(dist(dataExpr1),method="average")
traitColors<-numbers2colors(dataTraits,signed=FALSE)

png("./image/sample-subtype-cluster.png",width = 1800,height = 1600)
plotDendroAndColors(sampleTree2,traitColors,groupLabels=names(dataTraits),main="Sample dendrogram and trait heatmap",
          cex.colorLabels=1.5,cex.dendroLabels=1,cex.rowText=1)

dev.off()


powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr1, powerVector = powers, verbose = 5)

png("./image/step2-beta-value.png",width = 800,height = 600)
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
dev.off()

cor<-WGCNA::cor

net = blockwiseModules(dataExpr1, power = 12, maxBlockSize = nGenes,
                       TOMType ='unsigned', minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, saveTOMFileBase = "drought",
                       verbose = 3)
#time-consuming

table(net$colors)
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18 
# 132 3272 1434  915  674  602  573  522  510  479  474  383  232  215  141   63   58   57   40 


cor<-stats::cor
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
png("./image/step4-genes-modules.png",width = 800,height = 300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

nGenes = ncol(dataExpr1)
nSamples = nrow(dataExpr1)
design=read.table("TraitData3.txt", sep = '\t', header = T, row.names = 1)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr1, moduleColors)$eigengenes

MEs = orderMEs(MEs0);
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
png("./image/step5-Module-trait-relationships.png",width = 1800,height = 1000,res = 120)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(30), #greenWhiteRed(50)
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 cex.lab.x = 0.5,    
                 cex.lab.y = 0.5,
                 zlim = c(-1,1),
                 xLabelsAngle = 30,    
                 main = paste("Module-trait relationships"))
dev.off()



modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(dataExpr1, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(dataExpr1, design, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", colnames(design), sep = "")
names(GSPvalue) = paste("p.GS.", colnames(design), sep = "")

for (module in modNames) {
  column = match(module, modNames) 
  moduleGenes = moduleColors == module  

  png(paste0("./image/step6-Module_membership-gene_significance", module, ".png"),
      width = 800, height = 600)

  par(mfrow = c(1,1)) 
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for Basal",
                     main = paste("Module membership vs. gene significance\n", module),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

  dev.off()
}

################################trait-module plot################################
modNames = substring(names(MEs), 3) 
geneModuleMembership = as.data.frame(cor(dataExpr1, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

for (trait in colnames(design)) {
  traitData = as.data.frame(design[, trait])
  names(traitData) = trait

  geneTraitSignificance = as.data.frame(cor(dataExpr1, traitData, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", trait, sep = "")
  names(GSPvalue) = paste("p.GS.", trait, sep = "")

  for (module in modNames) {
    column = match(module, modNames)
    moduleGenes = moduleColors == module

    if (sum(moduleGenes) < 5) next

    file_name = paste0("./image/MM_vs_GS_", module, "_vs_", trait, ".png")
    png(file_name, width = 800, height = 600)

    par(mfrow = c(1,1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", trait),
                       main = paste("MM vs. GS\nModule:", module, "| Trait:", trait),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
                       col = module)

    dev.off()
  }
}


geneTree = net$dendrograms[[1]];  
dissTOM = 1-TOMsimilarityFromExpr(dataExpr1, power = 12); 

nSelect = 4000
set.seed(1)
selectGenes = sample(ncol(dataExpr1), size = nSelect)
selectTOM = dissTOM[selectGenes, selectGenes]
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[selectGenes]


png(filename = "./image/TOMplot_heatmap_selected_genes.png", width = 1200, height = 1000)
TOMplot(selectTOM^7, selectTree, selectColors,
        main = paste("Network heatmap plot (top", nSelect, "genes)"))
dev.off()

RA = as.data.frame(design[,31])
names(RA) = "Stearic acid"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, RA))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
png(filename = paste0("./imagestep7-Eigengene-dendrogram",names(RA),".png"),width = 800,height = 600)
plotEigengeneNetworks(MET, "", 
                        marDendro =c(0,5,1,5),
                        marHeatmap = c(5,6,1,2), cex.lab = 0.8,
                        xLabelsAngle = 90) 
dev.off()


## ---- Signed ME–trait correlation heatmap (shows negatives) ----
RA <- as.data.frame(design[, 15, drop = FALSE])
colnames(RA) <- "Indole-3-acetonitrile"

MET <- orderMEs(cbind(MEs, RA))

# Spearman（Pearson：method="pearson"）
MEcor <- cor(MET, use = "pairwise.complete.obs", method = "spearman")
MEp   <- corPvalueStudent(MEcor, nrow(MET))  # 样本量 = 行数

textMatrix <- paste(signif(MEcor, 2), "\n(",
                    signif(MEp, 1), ")", sep = "")
dim(textMatrix) <- dim(MEcor)

png("./image/step7-ME_trait_signedCorr_Indole-3-acetonitrile.png",
    width = 1200, height = 900, res = 120)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix      = MEcor,              
  xLabels     = colnames(MET),
  yLabels     = rownames(MEcor),
  ySymbols    = rownames(MEcor),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),      # （greenWhiteRed）
  textMatrix  = textMatrix,
  setStdMargins = FALSE,
  cex.text    = 0.4,
  cex.lab.x   = 0.7,
  cex.lab.y   = 0.9,
  xLabelsAngle = 45,                   
  zlim        = c(-1, 1),            
  main        = "Signed correlations: Module eigengenes vs Stearic acid"
)
dev.off()


################################every trait plot################################
for (trait in colnames(design)) {
  
  trait_df <- as.data.frame(design[, trait])
  names(trait_df) <- trait

  MET <- orderMEs(cbind(MEs, trait_df))

  file_name <- paste0("./image/step7-Eigengene-dendrogram", trait, ".png")

  # plot
  png(file_name, width = 800, height = 600)
  par(cex = 0.9)
  plotEigengeneNetworks(MET, "",
                        marDendro = c(0, 5, 1, 5),
                        marHeatmap = c(5, 6, 1, 2),
                        cex.lab = 0.8,
                        xLabelsAngle = 90)
  dev.off()
}

TOM = TOMsimilarityFromExpr(dataExpr1, power = 12);  ##time consuming
module = "brown";##"purple",blue/red
# Select module probes
probes = colnames(dataExpr1)
inModule = (moduleColors==module);
modProbes = probes[inModule];

modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

cyt = exportNetworkToCytoscape(  
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""), 
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""), 
    weighted = TRUE,
    threshold = 0.02,
    nodeNames = modProbes,
    nodeAttr = moduleColors[inModule]
  )



########################exportModuleNetworks() function############################
exportModuleNetworks <- function(dataExpr, moduleColors, power = 12,
                                 outDir = "WGCNA_Module_Export", threshold = 0.02) {
  # WGCNA package
  if (!require("WGCNA")) stop("Please install/load the WGCNA package first.")

  dir.create(outDir, showWarnings = FALSE)
  dir.create(file.path(outDir, "Cytoscape"), showWarnings = FALSE)
  dir.create(file.path(outDir, "GeneLists"), showWarnings = FALSE)

  if (length(moduleColors) == nrow(dataExpr)) {
    dataExpr_t <- dataExpr  
    } else if (length(moduleColors) == ncol(dataExpr)) {
    dataExpr_t <- t(dataExpr) 
  } else {
   stop("moduleColors not match!")
  }

  probes <- rownames(dataExpr_t)

  if (length(moduleColors) != nrow(dataExpr_t)) {
    stop("Length of moduleColors does not match number of genes in expression data (rows).")
  }

  TOM = TOMsimilarityFromExpr(dataExpr, power = power)

  allModules <- unique(moduleColors)
  allModules <- allModules[allModules != "grey"]

  for (module in allModules) {
    inModule <- (moduleColors == module)
    modProbes <- probes[inModule]

    if (length(modProbes) < 2) {
      next
    }

    modTOM <- TOM[inModule, inModule]

    if (length(modProbes) != nrow(modTOM)) {
      next
    }

    dimnames(modTOM) <- list(modProbes, modProbes)

    edgeFile <- file.path(outDir, "Cytoscape", paste0("CytoscapeInput-edges-", module, ".txt"))
    nodeFile <- file.path(outDir, "Cytoscape", paste0("CytoscapeInput-nodes-", module, ".txt"))
    geneListFile <- file.path(outDir, "GeneLists", paste0(module, "_genes.txt"))

    # export Cytoscape file
    tryCatch({
      exportNetworkToCytoscape(modTOM,
                               edgeFile = edgeFile,
                               nodeFile = nodeFile,
                               weighted = TRUE,
                               threshold = threshold,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])
    }, error = function(e) {
      cat("Errors", module, "failure", e$message, "\n")
    })

    # 导出基因名
    write.table(modProbes, file = geneListFile,
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  cat("Successed! The file is in：", normalizePath(outDir), "\n")
}
exportModuleNetworks(dataExpr = dataExpr1,
                     moduleColors = moduleColors,
                     power = 12,
                     outDir = "WGCNA_Export",
                     threshold = 0.02)


moduleColors <- net$colors
geneNames <- colnames(dataExpr1) 
str(geneNames)
moduleGenesList <- split(geneNames, moduleColors)

for (module in names(moduleGenesList)) {
  cat("Module:", module, "\n")
  cat("Genes:\n")
  cat(moduleGenesList[[module]], "\n\n")
}

for (module in names(moduleGenesList)) {
  write.table(moduleGenesList[[module]], 
              file = paste0("Module_", module, "_genes.txt"), 
              row.names = TRUE, col.names = TRUE, quote = FALSE)
}


MEs <- moduleEigengenes(dataExpr1, moduleColors)$eigengenes

for (module in names(moduleGenesList)) {
  moduleGenes <- moduleGenesList[[module]]
  moduleExpr <- dataExpr1[, moduleGenes, drop = FALSE]
  
  write.table(moduleExpr, 
              file = paste0("Module_", module, "_gene_expression.txt"), 
              row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  moduleEigengene <- MEs[, paste0("ME", module)]
  moduleEigengeneData <- data.frame(Sample = rownames(dataExpr1), Eigengene = moduleEigengene)
  write.table(moduleEigengeneData, 
              file = paste0("Module_", module, "_eigengene_expression.txt"), 
              row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")
}


moduleGenesList <- split(geneNames, moduleColors)
moduleIndices <- as.data.frame(table(moduleColors))

cat("Module Colors and Indices:\n")
print(moduleIndices)
