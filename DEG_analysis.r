###DE analysis for Day3 vs. Day0, Dy5 vs. Day0 and Day7 vs. Day0 for each celltype. 
library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(harmony)
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