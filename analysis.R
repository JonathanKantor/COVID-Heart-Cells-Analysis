# Import Libraries
library(Seurat)
library(ggplot2)
library(sctransform)
library(patchwork)
library(cowplot)
library(readxl)

# Function to process a single directory
process_directory <- function(directory) {
  m <- strsplit(directory, '_')
  cond <- m[[1]][3]
  patient <- paste(m[[1]][1], m[[1]][2], sep = '_')
  subdir <- paste("C:/Users/Jonathan/Downloads/IPSC_data/", directory, sep = "")

  if (cond == "Mock") {
    # Process Mock Data
    Mock.data <- Read10X(data.dir = paste(subdir, "/raw_feature_bc_matrix", sep = ""))
    Mock <- CreateSeuratObject(counts = Mock.data, min.cells = 3, min.features = 1, project = patient)
    Mock[["percent.mt"]] <- PercentageFeatureSet(Mock, pattern = "^MT-")
    Mock <- subset(Mock, subset = nFeature_RNA > 500)
    Mock$stim <- "Mock"
    Mock <- SCTransform(Mock, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
    return(list(Mock))
  } 
else if{
 # Process CoV-ILs Data
    CoV-ILs.data <- Read10X(data.dir = paste(subdir, "/raw_feature_bc_matrix", sep = ""))
    CoV-ILs <- CreateSeuratObject(counts = CoV-ILs.data, min.cells = 3, min.features = 1, project = patient)
    CoV-ILs[["percent.mt"]] <- PercentageFeatureSet(CoV-ILs, pattern = "^MT-")

    # Remove hCov-19 Gene
    indexhCov19 <- grep("hCov-19", rownames(CoV))
    indexKeep <- 1:length(rownames(CoV))
    indexKeep <- indexKeep[!indexKeep %in% indexhCov19]
    CoV-ILs <- CoV-ILs[indexKeep,]
    CoV-ILs <- subset(CoV-ILs, subset = nFeature_RNA > 500)
    CoV-ILs$stim <- "CoV-ILs"
    CoV-ILs <- SCTransform(CoV-ILs, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
    return(list(CoV-ILs))
}
else{
    # Process CoV Data
    CoV.data <- Read10X(data.dir = paste(subdir, "/raw_feature_bc_matrix", sep = ""))
    CoV <- CreateSeuratObject(counts = CoV.data, min.cells = 3, min.features = 1, project = patient)
    CoV[["percent.mt"]] <- PercentageFeatureSet(CoV, pattern = "^MT-")

    # Remove hCov-19 Gene
    indexhCov19 <- grep("hCov-19", rownames(CoV))
    indexKeep <- 1:length(rownames(CoV))
    indexKeep <- indexKeep[!indexKeep %in% indexhCov19]
    CoV <- CoV[indexKeep,]

    CoV <- subset(CoV, subset = nFeature_RNA > 500)
    CoV$stim <- "CoV"
    CoV <- SCTransform(CoV, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
    return(list(CoV))
  }
}

dirs <- list.files(path = "C:/Users/Jonathan/Downloads/IPSC_data/", recursive = F, full.names = F)
ipsc.list <- lapply(dirs, process_directory)

# Select Integration Features
features <- SelectIntegrationFeatures(object.list = ipsc.list, nfeatures = 3000)
ipsc.list <- PrepSCTIntegration(object.list = ipsc.list, anchor.features = features)

# Integrate Data
ipsc.anchors <- FindIntegrationAnchors(object.list = ipsc.list, normalization.method = "SCT", anchor.features = features)
ipsc.combined.sct <- IntegrateData(anchorset = ipsc.anchors, normalization.method = "SCT")

# Perform Dimensionality Reduction
ipsc.combined.sct <- RunPCA(ipsc.combined.sct, verbose = FALSE)
ipsc.combined.sct <- RunUMAP(ipsc.combined.sct, reduction = "pca", dims = 1:30)

# Find Neighbors and Cluster
ipsc.combined.sct <- FindNeighbors(ipsc.combined.sct, reduction = "pca", dims = 1:20)
ipsc.combined.sct <- FindClusters(ipsc.combined.sct, resolution = 0.5)

# Seperate Between Treatment Conditions
DimPlot(ipsc.combined.sct, reduction = "umap", split.by = "stim")

# Create Table of Cluster and Number of Cells
pt <- table(Idents(ipsc.combined.sct), ipsc.combined.sct$orig.ident, ipsc.combined.sct$stim)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)

# Calculate Percent of Condition
pt <- pt %>%
  group_by(Var3) %>%
  mutate(Sum = sum(Freq),
         PercentTotal = 100 * Freq / Sum) %>%
  ungroup()

#Split by Treatment
intCoV <- subset(x = ipsc.combined.sct, subset = stim == "CoV")
intCoV-ILs <- subset(x = ipsc.combined.sct, subset = stim == "CoV-ILs")
intMock <- subset(x = ipsc.combined.sct, subset = stim == "Mock")

#Percent hCov-19 per Cluster and per Treatment
VlnPlot(intCoV, features = "percent.cov") 
VlnPlot(intCoV-ILs, features = "percent.cov")

#Percent ACE2
VlnPlot(intCoV, features = "ACE2")
VlnPlot(intCoV-ILs, features = "ACE2")

#Plot Cluster Dendrogram
clustree(intMock)

#Calculating Mock Expression Markers
DefaultAssay(intMock) <- "RNA"
intMock<- NormalizeData(intMock)
all.markers <- FindAllMarkers(intMock, assay = "RNA", slot = “data”, min.pct = 0.25, thresh.use = 0.25)

#Filter Top 300 Lowest adj_p_value for each cluster
sortedMock <- head(all.markers, 300)
IPSCclusters = unique(all.markers$cluster)
for(i in IPSCclusters){
  indx <- match(i, all.markers$cluster)
  sortedMock <- rbind(sortedMock, all.markers[indx:(indx+299),])
  }
overlapComparison <- function(sortedMock, pathway, colnames) {
  num_sample_clusters <- length(unique(sortedMock$seurat_clusters))
  num_pathway_clusters <- length(colnames)

  pDf <- data.frame(matrix(ncol = 8, nrow = num_sample_clusters * num_pathway_clusters))
  colnames(pDf) <- c("Sample", "Pathway", "NumClusterMarkerGenes", "NumPathwayMarkerGenes",
                     "Overlap", "BackgroundGeneCount", "OverlapP-value", "-log(p-value)")

  count <- 1

  for (i in 0:num_sample_clusters-1) {
    ipscGenes <- sortedMock$gene[sortedMock$seurat_clusters == i]
    ipscNum <- length(ipscGenes)

    for (j in 0:num_pathway_clusters-1) {
      pathwayGenes <- pathway$Gene[pathway$Cluster == j]
      pathwayNum <- length(pathwayGenes)

      overlapGenes <- intersect(ipscGenes, pathwayGenes)
      overlapNum <- length(overlapGenes)
      backgroundNum <- length(setdiff(ipscGenes, pathwayGenes))

      # Create Contigency Table
      dat <- matrix(c(overlapNum, pathwayNum - overlapNum, ipscNum - overlapNum,
                      num_features - pathwayNum - ipscNum + overlapNum),
                    nrow = 2, byrow = TRUE)
      colnames(dat) <- c("Sample DE", "Sample not DE")
      rownames(dat) <- c("Pathway DE", "Pathway not DE")

      # Perform Fisher's Exact Test
      pTest <- fisher.test(dat)

      # Store results in pDf
      pDf$Sample[count] <- i
      pDf$Pathway[count] <- colnames[j + 1]
      pDf$NumClusterMarkerGenes[count] <- ipscNum
      pDf$NumPathwayMarkerGenes[count] <- pathwayNum
      pDf$Overlap[count] <- overlapNum
      pDf$BackgroundGeneCount[count] <- backgroundNum
      pDf$`OverlapP-value`[count] <- pTest$p.value
      pDf$`-log(p-value)`[count] <- -log10(pTest$p.value)

      count <- count + 1
    }
  }

  pMat <- matrix(pDf$`OverlapP-value`, nrow = num_sample_clusters, ncol = num_pathway_clusters)
  colnames(pMat) <- colnames

  return(list(pDf = pDf, pMat = pMat))
}

#Compare Against Patrick Ellinor Healthy Heart Data Set
HF <- read_excel("C:/Users/Jonathan/Downloads/mature.xlsx")
HF = subset(PE, Marker>0)
colnames <- c("Fibroblast I","Fibroblast II", "Atrial Cardiomyocyte", "Ventricular Cardiomyocyte I", "Ventricular Cardiomyocyte II", "Pericyte","Macrophage", "Endothelium I", "Endothelium II", "Adipocyte", "Vascular Smooth Muscle", "Fibroblast III", "Ventricular Cardiomyocyte III", "Neuronal", "Lymphocyte")
overlapComparison(sortedMock, HF, colnames)

# Cluster 5 Selected as Fibroblast-like Cardiomyocyte qnd Cluster 7 Selected as Ventricular Cardiomyocyte
covidMarkers <- function(selected_cluster, clusterMarkerData){
  combinedCluster <- subset(ipsc.combined.sct, seurat_cluster == selected_cluster)
  combinedCluster$nominatorStatus <- "null"
  for(i in 1:length(combinedCluster5$stim)){
    if(combinedCluster$percent.cov[i] > 1 & combinedCluster$percent.cov[i] < 62){
      combinedCluster$nominatorStatus[i] <- "nominator"
    }
    else if(combinedCluster$percent.cov<=1 || combinedCluster5$percent.cov == "null"){
      combinedCluster$nominatorStatus[i] <- "denominator"
    }
  }
  clusterMarkerData <- FindMarkers(combinedCluster, ident.1 = "nominator", ident.2 = “denominator”, group.by = "nominatorStatus", min.pct = 0.25)
  clusterMarkerData <- as.data.frame(row.names(clusterMarkerData))
  print(df[1], row.names = FALSE)
  }

  DefaultAssay(intCoV) <- "RNA"
  intCoV<- NormalizeData(intCoV)
  all.markersCoV <- FindAllMarkers(intCoV, assay = "RNA", slot = “data”, min.pct = 0.25, thresh.use = 0.25)
  cluster_5_markers_CoV <- all.markersCoV[all.markersCoVs$cluster == 5, ]
  cluster_7_markers_CoV <- all.markersCoV[all.markersCoV$cluster == 7, ]


#Finding DEG for Covid Infected Groups
  covidMarkers(5, cluster_5_markers_CoV)
  covidMarkers(7, cluster_7_markers_CoV)
  
EnhancedVolcano(cluster_5_markers_CoV , 
                lab = cluster_5_markers_CoV$gene,
                title = 'COVID Fibroblast-like Cluster 5',
                x ="avg_log2FC", 
                y ="p_val_adj")

EnhancedVolcano(cluster_7_markers_CoV , 
                lab = cluster_7_markers_CoV$gene,
                x ="avg_log2FC",
                title = 'COVID Cardiomyocyte Cluster 7',
                y ="p_val_adj")

DefaultAssay(intCoV-ILs) <- "RNA"
  intCoV<- NormalizeData(intCoV-ILs)
  all.markersCoV-ILs <- FindAllMarkers(intCoV, assay = "RNA", slot = “data”, min.pct = 0.25, thresh.use = 0.25)
  cluster_5_markers_CoV-ILs <- all.markersCoV-ILs[all.markersCoV-ILs$cluster == 5, ]
  cluster_7_markers_CoV-ILs <- all.markersCoV-ILs[all.markersCoV-ILs$cluster == 7, ]


#Finding DEG for Covid Infected Groups
  covidMarkers(5, cluster_5_markers_CoV-ILs)
  covidMarkers(7, cluster_7_markers_CoV-ILs)
  
EnhancedVolcano(cluster_5_markers_CoV-ILs , 
                lab = cluster_5_markers_CoV-ILs$gene,
                title = 'COVID + ILs Fibroblast-like Cluster 5',
                x ="avg_log2FC", 
                y ="p_val_adj")

EnhancedVolcano(cluster_7_markers_CoV-ILs , 
                lab = cluster_7_markers_CoV-ILs$gene,
                x ="avg_log2FC",
                title = 'COVID Cardiomyocyte Cluster 7',
                y ="p_val_adj")





