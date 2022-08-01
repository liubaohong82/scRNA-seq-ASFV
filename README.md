library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(harmony)

###Part1 
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

####Part2
library(Seurat)
library(ggplot2)
library(RColorBrewer)

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

##Part3
####################### MAGIC
 library(Rmagic)
 library(ggplot2)
scRNA_harmony_count = as.data.frame(scRNA_harmony[["RNA"]]@counts)
bmmsc = scRNA_harmony_count
# keep genes expressed in at least 10 cells
keep_cols <- colSums(bmmsc > 0) > 10
bmmsc <- bmmsc[,keep_cols]
keep_rows <- rowSums(bmmsc) > 200
bmmsc <- bmmsc[keep_rows,]
#Normalizeing
bmmsc <- library.size.normalize(bmmsc)
bmmsc <- sqrt(bmmsc)
bmmsc_MAGIC <- magic(bmmsc)
scRNA_MAGIC = bmmsc_MAGIC$result

#######################Scatter plot for CD14 and CD16 in virus infected cells 
 CD14_magic = which(colnames(scRNA_MAGIC)=="CD14") #a
 CD16_magic= which(colnames(scRNA_MAGIC)=="FCGR3A") #b

data_extracted = cbind(scRNA_meta_count_virus[,c(1:9)],scRNA_MAGIC[,c(CD14_magic,CD16_magic)])
data_extracted_mono =data_extracted[data_extracted[,9]=="Monocyte",]
data_extracted_mono_Day0 = data_extracted_mono[data_extracted_mono[,8]=="Day0",]
data_extracted_mono_Day3 = data_extracted_mono[data_extracted_mono[,8]=="Day3",]
data_extracted_mono_Day5 = data_extracted_mono[data_extracted_mono[,8]=="Day5",]
data_extracted_mono_Day7 = data_extracted_mono[data_extracted_mono[,8]=="Day7",]

library(ggplot2)
 pdf("Mono_CD14_CD16.pdf",width=16,height=4)
 p1 = ggplot(data_extracted_mono_Day0,aes(x=CD14,y=FCGR3A,color=Infect_status)) + geom_point(shape=19)+theme(legend.position="top") +  xlab("CD14") + ylab("CD16")+ ylim(0,2.25) + scale_color_brewer(palette = "Paired") 
 p2 = ggplot(data_extracted_mono_Day3,aes(x=CD14,y=FCGR3A,color=Infect_status)) + geom_point(shape=19)+theme(legend.position="top") +  xlab("CD14") + ylab("CD16")+ ylim(0,2.25) + scale_color_brewer(palette = "Paired")
 p3 = ggplot(data_extracted_mono_Day5,aes(x=CD14,y=FCGR3A,color=Infect_status)) + geom_point(shape=19)+theme(legend.position="top") +  xlab("CD14") + ylab("CD16")+ ylim(0,2.25) + scale_color_brewer(palette = "Paired")
 p4 = ggplot(data_extracted_mono_Day7,aes(x=CD14,y=FCGR3A,color=Infect_status)) + geom_point(shape=19)+theme(legend.position="top")+  xlab("CD14") + ylab("CD16")+ ylim(0,2.25) + scale_color_brewer(palette = "Paired")  
library(plyr)
library(gridExtra)
 grid.arrange(p1,p2,p3,p4,ncol=4,nrow=1)
 dev.off()

###Part4
##########################ISG scores
gene = read.table("ISG_score_genes.txt")
isg_features <- list(c(gene))
isg_features=gene
 scRNA.isg <- AddModuleScore( 
    object = scRNA_harmony,
   features = isg_features,
   ctrl = 41,
   name = 'ISG_Features'
 )

#scRNA_ISG_score = read.table("scRNA_ISG_score.txt",header=T,sep="\t",row.names=1)
#ISG_data = scRNA_ISG_score
pdf("ISG_celltype_boxplot.pdf",height=3.5,width=10)
ggplot(data = scRNA.isg@meta.data) + geom_boxplot(aes(x = celltype1, y = ISG_Features1,fill=factor(group)))+
scale_fill_manual(name = "", values = c(brewer.pal(7,"Paired")[c(1,2,4,5)]))
dev.off()

pdf("ISG_baseline_infected_bystanders_mono_boxplot.pdf",height=3.5,width=10)
ggplot(data = ISG_boxplot) + geom_boxplot(aes(x = Infect_status, y = ISG_Features1,fill=factor(Infect_status)))+
scale_fill_manual(name = "", values = c(brewer.pal(7,"Paired")[c(1,2,4)]))
dev.off()

aa=scRNA.isg@meta.data[scRNA.isg@meta.data[,8]=="Monocyte",]

pdf("ISG_baseline_infected_bystanders_mono_boxplot_0411.pdf",height=3.5,width=4)
 #scRNA_harmony= scRNA.isg
ggplot(data = aa) + geom_boxplot(aes(x = Infect_status, y = ISG_Features1,fill=factor(Infect_status))) +
scale_fill_manual(name = "", values = c(brewer.pal(9,"Paired")[c(1,2,3,4,5)]))
dev.off()
pdf("Magic_baseline_infected_bystanders_mono_scatterPlot_0411.pdf",width=4.5,height=3.5)
ggplot(aa,aes(x=CD14,y=FCGR3A,color=Infect_status)) + geom_point(shape=19)+theme(legend.position="right") +  xlab("CD14") + ylab("CD16")+ ylim(0,2.25) + 
scale_color_brewer(palette = "Paired") # scale_y_continuous(limits=c(0, 2.25), breaks=0.25)
dev.off()
 grid.arrange(p2,p1,ncol=2,nrow=1)
dev.off()

pdf("correlation_ISG_Virus_scatterPlot_0411.pdf",width=4,height=4)
plot(log10(ISG_boxplot[,19]),ISG_boxplot[,20],ylab="ISG Score",xlab="log10 sum(Virus UMI count)",pch=16,col="#A6CEE3")
text(3.5,1,"p < 2.2e-16")
text(3.5,1.15,"r = -0.32")


###Part4
###DE analysis for Day3 vs. Day0, Dy5 vs. Day0 and Day7 vs. Day0 for each celltype. 
scRNA_harmony@meta.data$celltype_condition=paste(scRNA_harmony@meta.data$celltype,scRNA_harmony@meta.data$group,sep = "_")
## The following code is for Day3 vs. Day0. 
marker_condition=data.frame()
Idents(scRNA_harmony)="celltype_condition"
for ( ci in sort(as.character(unique(scRNA_harmony@meta.data$celltype))) ) {
  tmp.marker <- FindMarkers(
    scRNA_harmony, logfc.threshold = 0.25, min.pct = 0.1,
    only.pos = F, test.use = "MAST",
    ident.1=paste0(ci,"_Day3"),ident.2=paste0(ci,"_Day0")
  )
  
  tmp.marker$gene=rownames(tmp.marker)
  tmp.marker$condition=ifelse(tmp.marker$avg_log2FC > 0,paste0(ci,"_Day3"),paste0(ci,"_Day0"))
  tmp.marker$cluster=ci
  
  tmp.marker=tmp.marker%>%filter(p_val_adj < 0.01)
  tmp.marker=as.data.frame(tmp.marker)
  tmp.marker=tmp.marker%>%arrange(desc(avg_log2FC))
  
  marker_condition=marker_condition%>%rbind(tmp.marker)
}
#table(marker_condition$condition)
#write.table(marker_condition,file = "markers.BasedOncondition_log2fc0.25_padj0.01.txt",sep="\t",row.names = F,col.names = T)

celltype.Day3=marker_condition[str_detect(marker_condition$condition,"_Day3"),]
celltype.Day0=marker_condition[str_detect(marker_condition$condition,"_Day0"),]
dim(marker_condition)
dim(celltype.Day3)
dim(celltype.Day0)


### 
celltype.order.name=names(sort(table(celltype.Day0$cluster))) 
celltype.order=1:length(unique(scRNA_harmony@meta.data$celltype))
names(celltype.order)=celltype.order.name

gene.shared=sort(table(celltype.Day0$gene))[sort(table(celltype.Day0$gene)) >= 2]
gene.shared=names(gene.shared)
gene.unique=sort(table(celltype.Day0$gene))[sort(table(celltype.Day0$gene)) < 2]
gene.unique=names(gene.unique)

Day0.gene.shared=celltype.Day0[celltype.Day0$gene %in% gene.shared,c("gene","condition","cluster")]

gene.shared.freq=as.data.frame(sort(table(Day0.gene.shared$gene)))
colnames(gene.shared.freq)=c("gene","freq")
gene.shared.freq$gene=as.character(gene.shared.freq$gene)
gene.order=c()
for (freqi in rev(unique(gene.shared.freq$freq))) {
  one.freq.df=as.data.frame(gene.shared.freq %>% filter(freq == freqi))
  if (dim(one.freq.df)[1] > 1) {
    one.gene.left=c()
    for (genei in one.freq.df[,"gene"]) {
      one.gene.distribution=Day0.gene.shared[Day0.gene.shared$gene %in% genei,]
      one.gene.left=append(one.gene.left,min(celltype.order[one.gene.distribution$cluster]))
    }
    one.gene.leftcelltype.df=as.data.frame(one.gene.left)
    one.gene.leftcelltype.df$gene=one.freq.df[,"gene"]
    one.gene.leftcelltype.df=one.gene.leftcelltype.df%>%arrange(one.gene.left)
    gene.order=append(gene.order,one.gene.leftcelltype.df$gene)
  }else{
    gene.order=append(gene.order,one.freq.df[1,"gene"])
  }
}

Day0.gene.unique=celltype.Day0[celltype.Day0$gene %in% gene.unique,c("gene","condition","cluster")]
Day0.gene.unique$cluster=factor(Day0.gene.unique$cluster,levels = celltype.order.name)
Day0.gene.unique=Day0.gene.unique%>%arrange(cluster)
gene.order=append(gene.order,Day0.gene.unique$gene)

Day0.plot.data=celltype.Day0[,c("gene","condition","cluster")]
Day0.plot.data$gene=factor(Day0.plot.data$gene,levels = rev(gene.order))
Day0.plot.data$cluster=factor(Day0.plot.data$cluster,levels = celltype.order.name)


p1=Day0.plot.data%>%ggplot(aes(x=cluster,y=gene))+
  geom_stripped_cols(odd = "#d9d9d9",even ="white",alpha=0.5)+
  geom_tile(aes(width=0.9),color="#A6CEE3",fill="#A6CEE3")+
  geom_hline(yintercept = length(gene.unique)+0.5,linetype=5)+
  scale_y_discrete(paste0("upper: ",length(gene.shared)," genes; lower: ",length(gene.unique)," genes"))+
  scale_x_discrete(expand = c(0,0))+
  labs(title = "Day3 downregulated")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.y.left = element_text(size = 16),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    plot.title = element_text(size = 18,hjust = 0.5)
  )

lower.text.df=as.data.frame(table(Day0.gene.unique$cluster))
Day0.gene.shared$cluster=factor(Day0.gene.shared$cluster,levels = celltype.order.name)
upper.text.df=as.data.frame(table(Day0.gene.shared$cluster))

p1+geom_text(data = lower.text.df,mapping = aes(x=Var1,y=10,label=Freq),size=5)+
  geom_text(data = upper.text.df,mapping = aes(x=Var1,y=10+length(gene.unique),label=Freq),size=5)

ggsave("Day3_downregulated.pdf",width = 10,height = 20,units = "cm")

celltype.order.name.copy=celltype.order.name 
celltype.order.name=celltype.order.name.copy 
celltype.order=1:length(unique(scRNA_harmony@meta.data$celltype))
names(celltype.order)=celltype.order.name

gene.shared=sort(table(celltype.Day3$gene))[sort(table(celltype.Day3$gene)) >= 2]
gene.shared=names(gene.shared)
gene.unique=sort(table(celltype.Day3$gene))[sort(table(celltype.Day3$gene)) < 2]
gene.unique=names(gene.unique)

Day3.gene.shared=celltype.Day3[celltype.Day3$gene %in% gene.shared,c("gene","condition","cluster")]

gene.shared.freq=as.data.frame(sort(table(Day3.gene.shared$gene)))
colnames(gene.shared.freq)=c("gene","freq")
gene.shared.freq$gene=as.character(gene.shared.freq$gene)
gene.order=c()
for (freqi in rev(unique(gene.shared.freq$freq))) {
  one.freq.df=as.data.frame(gene.shared.freq %>% filter(freq == freqi))
  if (dim(one.freq.df)[1] > 1) {
    one.gene.left=c()
    for (genei in one.freq.df[,"gene"]) {
      one.gene.distribution=Day3.gene.shared[Day3.gene.shared$gene %in% genei,]
      one.gene.left=append(one.gene.left,min(celltype.order[one.gene.distribution$cluster]))
    }
    one.gene.leftcelltype.df=as.data.frame(one.gene.left)
    one.gene.leftcelltype.df$gene=one.freq.df[,"gene"]
    one.gene.leftcelltype.df=one.gene.leftcelltype.df%>%arrange(one.gene.left)
    gene.order=append(gene.order,one.gene.leftcelltype.df$gene)
  }else{
    gene.order=append(gene.order,one.freq.df[1,"gene"])
  }
}

Day3.gene.unique=celltype.Day3[celltype.Day3$gene %in% gene.unique,c("gene","condition","cluster")]
Day3.gene.unique$cluster=factor(Day3.gene.unique$cluster,levels = celltype.order.name)
Day3.gene.unique=Day3.gene.unique%>%arrange(cluster)
gene.order=append(gene.order,Day3.gene.unique$gene)

Day3.plot.data=celltype.Day3[,c("gene","condition","cluster")]
Day3.plot.data$gene=factor(Day3.plot.data$gene,levels = rev(gene.order))
Day3.plot.data$cluster=factor(Day3.plot.data$cluster,levels = celltype.order.name)

p1=Day3.plot.data%>%ggplot(aes(x=cluster,y=gene))+
  geom_stripped_cols(odd = "#d9d9d9",even ="white",alpha=0.5)+
  geom_tile(aes(width=0.9),color="#FB9A99",fill="#FB9A99")+
  geom_hline(yintercept = length(gene.unique)+0.5,linetype=5)+
  scale_y_discrete(paste0("upper: ",length(gene.shared)," genes; lower: ",length(gene.unique)," genes"))+
  scale_x_discrete(expand = c(0,0))+
  labs(title = "Day3 upregulated")+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    axis.title.y.left = element_text(size = 16),
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,size = 14,color = "black"),
    axis.title.x.bottom = element_blank(),
    plot.title = element_text(size = 18,hjust = 0.5)
  )

lower.text.df=as.data.frame(table(Day3.gene.unique$cluster))
Day3.gene.shared$cluster=factor(Day3.gene.shared$cluster,levels = celltype.order.name)
upper.text.df=as.data.frame(table(Day3.gene.shared$cluster))

p1+geom_text(data = lower.text.df,mapping = aes(x=Var1,y=12,label=Freq),size=5)+
  geom_text(data = upper.text.df,mapping = aes(x=Var1,y=12+length(gene.unique),label=Freq),size=5)
ggsave("Day3_upregulated.pdf",width = 10,height = 20,units = "cm")

####DE analysis for infected and bystanders for Monocyte and Macrophage
Monocyte.sub <- subset(scRNA_harmony, celltype=="Monocyte")
Idents(Monocyte.sub)="Infect_status"
diff_mono_Infected_vs_Bystanders  <- FindMarkers(Monocyte.sub,
                     ident.1 = "Infected", 
                     ident.2 = "Bystanders",assay = "RNA",test.use="MAST")

Monocyte.sub.day5 <- subset(Monocyte.sub, group=="Day5")
diff_mono_Infected_vs_Bystanders_day5  <- FindMarkers(Monocyte.sub.day5,
                     ident.1 = "Infected", 
                     ident.2 = "Bystanders",assay = "RNA",test.use="MAST")

Monocyte.sub.day7 <- subset(Monocyte.sub, group=="Day7")
diff_mono_Infected_vs_Bystanders_day7  <- FindMarkers(Monocyte.sub.day7,
                     ident.1 = "Infected", 
                     ident.2 = "Bystanders",assay = "RNA",test.use="MAST")

write.table(diff_mono_Infected_vs_Bystanders,"diff_mono_Infected_vs_Bystanders.txt",sep="\t",quote=F)
write.table(diff_mono_Infected_vs_Bystanders_day7,"diff_mono_Infected_vs_Bystanders_day7.txt",sep="\t",quote=F)
write.table(diff_mono_Infected_vs_Bystanders_day5,"diff_mono_Infected_vs_Bystanders_day5.txt",sep="\t",quote=F)

Macrophage.sub <- subset(scRNA_harmony, celltype=="Macrophage")
Idents(Macrophage.sub)="Infect_status"
diff_macro_Infected_vs_Bystanders  <- FindMarkers(Macrophage.sub,
                     ident.1 = "Infected", 
                     ident.2 = "Bystanders",assay = "RNA",test.use="MAST")
write.table(diff_macro_Infected_vs_Bystanders,"diff_macro_Infected_vs_Bystanders.txt",sep="\t",quote=F)

















 