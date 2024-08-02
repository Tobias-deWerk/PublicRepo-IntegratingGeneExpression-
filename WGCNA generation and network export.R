#Manual Network Construction

#Setting up session
workingDir = 'C:/WGCNA/output for visualization2'
setwd(workingDir)
library(WGCNA)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Load the data saved in the first part
lnames = load(file = "Drought ExprData TraitData.RData")
lnames

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(ExprData, powerVector = powers, verbose = 5, networkType = 'signed')

# Plot the results:
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
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Best power is 20
softPower = 20;
adjacency = adjacency(ExprData, power = softPower, type = 'signed')

# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency, TOMType = 'signed');
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 60;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(ExprData, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.01
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(ExprData, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "ManualData DroughtSeries NoMerge - signed soft20.RData")

#Export to cytoscape - PER MODULE
# Export the gene list of new modules 
for (i in 1:length(MEs)){
  modules = c(substring(names(MEs)[i], 3));
  genes = colnames(ExprData)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = 0.3, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}


#Cytoscape stuff
library(RCy3)
geneInfo = read.csv(file = "geneInfocsv.csv", sep = ",")

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

edge <- read.delim("output_for_cytoscape/merge_CytoscapeInput-edges-greenyellow.txt")
node <- read.delim("output_for_cytoscape/merge_CytoscapeInput-nodes-greenyellow.txt")

colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName") 
colnames(node) <- c("id","altName","node_attributes")
 
nodeWithInfo <- merge(node, geneInfo, by.x = "id", by.y = "gene", all.x = TRUE)

str(edge)
str(node)
str(nodeWithInfo)

createNetworkFromDataFrames(nodes = nodeWithInfo, edges = edge, title = "greenyellow_With_Info", collection = "DataFrame Example")

#Export to cytoscape - TFs of interest
library(RCy3)
threshold <- 0.4

TFlist = c("Solyc08g007830.1", "Solyc08g078180.1", "Solyc02g078285.1", "Solyc12g007070.2", 
           "Solyc02g083220.2", "Solyc02g089340.4", "Solyc07g063420.3", "Solyc07g054220.1", 
           "Solyc08g081650.3", "Solyc04g009450.1", "Solyc02g085560.3", "Solyc06g083505.2", 
           "Solyc03g044300.3", "Solyc01g109880.3", "Solyc03g097650.3", "Solyc01g106040.3", 
           "Solyc02g089790.4", "Solyc06g034340.3", "Solyc05g010516.1", "Solyc08g079120.4", 
           "Solyc03g093540.1", "Solyc11g044560.3", "Solyc06g060940.2", "Solyc03g097320.3", 
           "Solyc01g081490.3", "Solyc02g037530.3", "Solyc02g090310.1", "Solyc03g119580.1", 
           "Solyc10g078610.1", "Solyc04g015360.3", "Solyc01g104650.3", "Solyc08g075940.4", 
           "Solyc09g066010.3", "Solyc12g042070.3", "Solyc03g121940.3", "Solyc06g062900.4", 
           "Solyc03g093550.1", "Solyc02g088070.3", "Solyc05g052410.3", "Solyc07g007120.3", 
           "Solyc01g087240.3", "Solyc07g053590.5", "Solyc11g068620.2", "Solyc03g093560.1", 
           "Solyc02g021680.3", "Solyc12g008830.3", "Solyc12g015710.2", "Solyc01g087990.3", 
           "Solyc03g115850.3", "Solyc03g097120.3", "Solyc04g078640.3", "Solyc09g009760.1", 
           "Solyc02g088180.3", "Solyc06g059970.4", "Solyc04g016000.3", "Solyc12g015640.2", 
           "Solyc03g114440.1", "Solyc06g071820.3", "Solyc07g047960.3", "Solyc02g092090.3", 
           "Solyc12g009240.1", "Solyc02g084350.3", "Solyc07g042230.1", "Solyc02g089210.4", 
           "Solyc03g026280.3", "Solyc04g071770.3", "Solyc06g050520.3", "Solyc06g066540.1", 
           "Solyc03g120840.3", "Solyc10g006880.3", "Solyc07g063410.3", "Solyc01g100460.3", 
           "Solyc04g005610.3", "Solyc06g074320.3", "Solyc01g095460.3", "Solyc06g072710.3", 
           "Solyc06g075370.3", "Solyc03g124110.2", "Solyc02g072560.1", "Solyc06g069710.3", 
           "Solyc12g013620.2", "Solyc04g016290.3", "Solyc12g044390.3", "Solyc06g053960.3", 
           "Solyc08g081960.3", "Solyc04g078840.3", "Solyc10g081840.3", "Solyc07g040680.3", 
           "Solyc10g083450.3", "Solyc02g087920.3", "Solyc09g005610.4", "Solyc05g056477.1", 
           "Solyc05g050220.3", "Solyc10g079150.3")

for (tf in TFlist) {
  # Find the index of the TF in the gene list
  tfIndex <- which(genes == tf)
  
    # Extract connections for this TF with a TOM score above a threshold, e.g., 0.2
  connectedIndices <- which(TOM[tfIndex, ] > threshold)
  
  # Get the gene names of the connected genes
  connectedGenes <- genes[connectedIndices]
  
  # Prepare the edge list for the subnetwork
  edges <- data.frame(
    source = rep(tf, length(connectedGenes)),
    target = connectedGenes,
    weight = TOM[tfIndex, connectedIndices]
  )
  
  # Prepare the node list for the subnetwork, including the TF itself
  nodes <- data.frame(
    id = c(tf, connectedGenes)
    # Add additional node attributes here, if available
  )
  
  # Export the node and edge lists to files
  write.table(nodes, file = paste0("nodes_", tf, ".txt"), row.names = FALSE, sep = "\t")
  write.table(edges, file = paste0("edges_", tf, ".txt"), row.names = FALSE, sep = "\t")
}

geneInfo = read.csv(file = "geneInfocsv.csv", sep = ",")

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

Day4_list = c("Solyc06g059970.4", "Solyc05g010516.1", "Solyc02g090310.1")

for (tf in 1:length(Day1_list)){

 TF_ID = Day1_list[tf]
 edge <- read.delim(paste0("edges_", TF_ID, ".txt"))
 node <- read.delim(paste0("nodes_", TF_ID, ".txt"))

 createNetworkFromDataFrames(nodes = node, edges = edge, title = TF_ID, collection = "DataFrame Example")

}
tf = 1
