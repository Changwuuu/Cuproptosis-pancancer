RNC_cu_pancan=read.table(file = "RNC_cu_pancan.txt",sep = "\t",header = T,check.names = F)

RNC_cu_pancan=as.data.frame(t(RNC_cu_pancan))
colnames(RNC_cu_pancan)=RNC_cu_pancan[1,]
RNC_cu_pancan=RNC_cu_pancan[-1,] 

#####cluster, it will need long time

library("ConsensusClusterPlus")

data=as.matrix(t(RNC_cu_pancan))

results = ConsensusClusterPlus(data,maxK=10,reps=1000,pItem=0.8,pFeature=1,
                               clusterAlg="pam",distance="pearson",corUse="everything",seed=1262118388.71279,plot="pdf",writeTable = T)


icl <- calcICL(results, 
               plot = "pdf")

#########
RNC_cu_pancan=tibble::rownames_to_column(RNC_cu_pancan)
colnames(RNC_cu_pancan)[1]="sample"


library(pheatmap)
cu_cluster=merge(k2,RNC_cu_pancan,by="sample")
cu_cluster=cu_cluster[order(cu_cluster[,2]),]

gsva_matrix1=cu_cluster[,-2]
rownames(gsva_matrix1)=gsva_matrix1[,1]
gsva_matrix1=gsva_matrix1[,-1]
write.table(gsva_matrix1,"gsva_matrix1.txt", sep="\t", quote=FALSE, row.names = TRUE)

gsva_matrix1=read.table(file = "gsva_matrix1.txt",sep = "\t",header = T,check.names = F,row.names = 1)

gsva_matrix1<- t(scale(gsva_matrix1))

gsva_matrix1[gsva_matrix1< -2] <- -2
gsva_matrix1[gsva_matrix1>2] <- 2

annotation_col = data.frame(Group=cu_cluster$Cluster)
rownames(annotation_col)<-cu_cluster$sample

pheatmap(gsva_matrix1,
         show_colnames = F,
         cluster_row = T,cluster_cols = F,colorRampPalette(colors = c("blue","white","red"))(100),
         annotation_col = annotation_col,
         cellwidth=0.015,cellheight=15,
         fontsize=8,treeheight_row=20,gaps_col=c(7125),
         width =7)
write.table(cu_cluster,"cu_cluster.txt", sep="\t", quote=FALSE, row.names = F, col.names = T)


#####combine cluster, cuscore and clinical information. Using Sangerbox to plot box plots and KM curve. Please cite:Shen, et al. 2022. Sangerbox: A comprehensive, interaction-friendly clinical bioinformatics analysis platform. iMeta 1(3): e36. https://doi.org/10.1002/imt2.36
cu_cluster_score=merge(cu_cluster,RNA_cu_clinical_cuscore,by="sample")
write.table(cu_cluster_score,"cu_cluster_score.txt", sep="\t", quote=FALSE, row.names = F, col.names = T)
