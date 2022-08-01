library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(harmony)

#cell type annotation
markers= c("MS4A1","MZB1","CD3E","CD3G","NKG7",'KLRB1','NCR1','CD68','C1QA','C1QB',"C1QC","LYZ","S100A9", "S100A8","FLT3","XCR1","IRF8","CXCL8","SELL")
pdf("DotPlot_cluster.pdf",width=10,height=5)
DotPlot(scRNA_harmony, group.by="cluster_dot",features = markers, cols = c('yellow', 'red')) + RotatedAxis() + coord_flip()
dev.off()

pdf("ViolinPlot_cluster.pdf",width=6,height=6)
VlnPlot(scRNA_harmony,group.by="RNA_snn_res.0.8",features = markers,pt.size = 0,stack = T,flip=T) &
  scale_fill_manual(values = c(brewer.pal(12,"Set3"),brewer.pal(5,"Dark2"),brewer.pal(5,"Dark2")))&
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(angle = 45,size = 10,hjust = 0,vjust = 0),
    legend.position = "none"
  )
dev.off()
###########
markers= c("MS4A1","MZB1","CD3E","NKG7",'CD68',"LYZ","FLT3","IRF8","CXCL8")
pdf("FeaturePlot_scRNA_harmony.pdf",width=7,height=7)
FeaturePlot(scRNA_harmony, features = markers,cols = c("lightgrey", "#ff0000"),combine=T)
dev.off()

#################celltype 
####vlnplot
library(ggplot2)
a= c("MS4A1","MZB1","CD3E","CD3G","NKG7",'KLRB1','NCR1','CD68','C1QA','C1QB',"C1QC","LYZ","S100A9", "S100A8","XCR1","IRF8","CXCL8","SELL")
pdf("ViolinPlot_celltype_ordered.pdf",width=6,height=6)
VlnPlot(scRNA_harmony,group.by="celltype2",features = a,pt.size = 0,stack = T,flip=F) &
  scale_fill_manual(values = c("#A6CEE3","#A6CEE3","#1F78B4","#1F78B4","#B2DF8A","#B2DF8A",
  "#B2DF8A","#33A02C","#33A02C","#33A02C","#33A02C","#FB9A99","#FB9A99",
  "#FB9A99","#E31A1C","#E31A1C","#FDBF6F","#FDBF6F"))&
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(angle = 90,size = 8,hjust = 0,vjust = 0),
    legend.position = "none"
  )
dev.off()

############Umap
pdf("Umap_combine_Paired_nolabel.pdf",width=7.5,height=6)
DimPlot(object = scRNA_harmony,label=F,label.size = 3,pt.size=0.6,cols=brewer.pal(9, "Paired"),group.by="celltype")
dev.off()

#########################Cell Abundance
cellcount_sample=table(scRNA_harmony$orig.ident,scRNA_harmony@meta.data$celltype)
sum=apply(cellcount_sample,1,sum)
cbind(cellcount_sample,sum)->cellcount_sample
cellcount_sample[,1:8]/cellcount_sample[,9]->Percent_samples
write.table(Percent_samples,"cell_Percent_samples.txt",sep="\t",quote=F)

cellcount_group=table(scRNA_harmony$group,scRNA_harmony@meta.data$celltype)
sum=apply(cellcount_group,1,sum)
cbind(cellcount_group,sum)->cellcount_group
cellcount_group[,1:8]/cellcount_group[,9]->Percent_group
write.table(Percent_group,"cell_Percent_group.txt",sep="\t",quote=F)