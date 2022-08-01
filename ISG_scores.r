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
