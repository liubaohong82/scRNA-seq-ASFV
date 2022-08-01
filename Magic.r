library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(harmony)
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