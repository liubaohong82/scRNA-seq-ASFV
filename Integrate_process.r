library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(harmony)

dir <- c("Day0_248","Day0_249","Day0_258","Day3_246","Day3_250","Day3_256","Day5_244","Day5_251","Day5_252","Day7_254","Day7_255","Day7_257")
dir <- paste0("/liubh/scRNA_seq/",dir)
sample_name <- c("Day0_248","Day0_249","Day0_258","Day3_246","Day3_250","Day3_256","Day5_244","Day5_251","Day5_252","Day7_254","Day7_255","Day7_257")

################QC##########################
scRNAlist <- list()
for(i in 1:length(dir)){
counts <- Read10X(data.dir = dir[i])
scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name[i], min.cells=3, min.features = 200) 
scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = percent.mt < 5) 
}   
#saveRDS(scRNAlist, "scRNAlist.rds")
###QC
sceList = lapply(scRNAlist, function(x) {
  subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 400 & nCount_RNA < 20000)
})

##saveRDS(sceList, "sceList.rds")

### Normalization and variable feature selection
scRNA_harmony <- merge(sceList[[1]], y=c(sceList[[2]], sceList[[3]], sceList[[4]],sceList[[5]], sceList[[6]], sceList[[7]],sceList[[8]], sceList[[9]], sceList[[10]],sceList[[11]], sceList[[12]]))
scRNA_harmony <- NormalizeData(scRNA_harmony) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

##Harmony
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")

#Dimensional reduction and Cell cluster identification
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30) %>% FindClusters()

##Visulization 
Idents(scRNA_harmony)="seurat_clusters"
plot1 = DimPlot(scRNA_harmony, reduction = "umap", label=T) 
plot2 = DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident') 
plotc <- plot1+plot2
ggsave("scRNA_harmony.pdf", plot = plotc, width = 10, height = 5)
saveRDS(scRNA_harmony, 'scRNA_harmony.rds')