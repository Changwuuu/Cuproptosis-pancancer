library(data.table)
RNASeq2=fread("RNASeq2.txt", data.table=FALSE)
rownames(RNASeq2)=RNASeq2$sample

RNASeq2=RNASeq2[,-1]


#######GSVA
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(GSVA)    
library(pheatmap)

##########input data ,construct gmt
##
library(openxlsx)
hallmark_cu_genesets<- read.xlsx("hallmark_cu_genesets.xlsx",sheet = "Sheet1")


sel_gmt=hallmark_cu_genesets
head(sel_gmt)
dim(sel_gmt)
sets=as.list(sel_gmt)
sets=lapply(sets, function(x) x[!is.na(x)])
#sets[1]

#########calculation GSVA score

exprMatrix=RNASeq2

head(exprMatrix)
dim(exprMatrix)
gsva_matrix<- gsva(as.matrix(exprMatrix), sets,method='gsva',kcdf='Gaussian',abs.ranking=F)
gsva_matrix1<- t(scale(t(gsva_matrix)))
#head(gsva_matrix1)  normalization
normalization<-function(x){
  return((x-min(x))/(max(x)-min(x)))}
nor_gsva_matrix1 <- normalization(gsva_matrix1) 
dim(nor_gsva_matrix1)

score.gsva.metab=as.data.frame(t(nor_gsva_matrix1))
head(score.gsva.metab)