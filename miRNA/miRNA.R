library(data.table)
miRNA_expr=fread("miRNA_EXP.txt", data.table=FALSE)

library(openxlsx)

miRNA_mRNA <- read.xlsx("miRNA_mRNA.xlsx",sheet = "Sheet1")
RNA_cu_clinical<- read.xlsx("RNA_cu_clinical.xlsx",sheet = "Sheet1")

FDX1=miRNA_mRNA[,1]
colnames(FDX1)[1]="sample"
FDX1=merge(FDX1,miRNA_expr,by="sample")
FDX1=as.data.frame(t(FDX1))
FDX1=tibble::rownames_to_column(FDX1)
colnames(FDX1)=FDX1[1,]
FDX1=FDX1[-1,]
FDX1=merge(FDX1,RNA_cu_clinical,by="sample")
write.table(FDX1,file="FDX1.txt",col.names = T,row.names = F,sep = "\t", quote = F)

LIAS=miRNA_mRNA[,2]
colnames(LIAS)[1]="sample"
LIAS=merge(LIAS,miRNA_expr,by="sample")
LIAS=as.data.frame(t(LIAS))
LIAS=tibble::rownames_to_column(LIAS)
colnames(LIAS)=LIAS[1,]
LIAS=LIAS[-1,]
LIAS=merge(LIAS,RNA_cu_clinical,by="sample")
write.table(LIAS,file="LIAS.txt",col.names = T,row.names = F,sep = "\t", quote = F)

LIPT1=miRNA_mRNA[,3]
colnames(LIPT1)[1]="sample"
LIPT1=merge(LIPT1,miRNA_expr,by="sample")
LIPT1=as.data.frame(t(LIPT1))
LIPT1=tibble::rownames_to_column(LIPT1)
colnames(LIPT1)=LIPT1[1,]
LIPT1=LIPT1[-1,]
LIPT1=merge(LIPT1,RNA_cu_clinical,by="sample")
write.table(LIPT1,file="LIPT1.txt",col.names = T,row.names = F,sep = "\t", quote = F)

DLD=miRNA_mRNA[,4]
colnames(DLD)[1]="sample"
DLD=merge(DLD,miRNA_expr,by="sample")
DLD=as.data.frame(t(DLD))
DLD=tibble::rownames_to_column(DLD)
colnames(DLD)=DLD[1,]
DLD=DLD[-1,]
DLD=merge(DLD,RNA_cu_clinical,by="sample")
write.table(DLD,file="DLD.txt",col.names = T,row.names = F,sep = "\t", quote = F)

DLAT=miRNA_mRNA[,5]
colnames(DLAT)[1]="sample"
DLAT=merge(DLAT,miRNA_expr,by="sample")
DLAT=as.data.frame(t(DLAT))
DLAT=tibble::rownames_to_column(DLAT)
colnames(DLAT)=DLAT[1,]
DLAT=DLAT[-1,]
DLAT=merge(DLAT,RNA_cu_clinical,by="sample")
write.table(DLAT,file="DLAT.txt",col.names = T,row.names = F,sep = "\t", quote = F)

PDHA1=miRNA_mRNA[,6]
colnames(PDHA1)[1]="sample"
PDHA1=merge(PDHA1,miRNA_expr,by="sample")
PDHA1=as.data.frame(t(PDHA1))
PDHA1=tibble::rownames_to_column(PDHA1)
colnames(PDHA1)=PDHA1[1,]
PDHA1=PDHA1[-1,]
PDHA1=merge(PDHA1,RNA_cu_clinical,by="sample")
write.table(PDHA1,file="PDHA1.txt",col.names = T,row.names = F,sep = "\t", quote = F)

PDHB=miRNA_mRNA[,7]
colnames(PDHB)[1]="sample"
PDHB=merge(PDHB,miRNA_expr,by="sample")
PDHB=as.data.frame(t(PDHB))
PDHB=tibble::rownames_to_column(PDHB)
colnames(PDHB)=PDHB[1,]
PDHB=PDHB[-1,]
PDHB=merge(PDHB,RNA_cu_clinical,by="sample")
write.table(PDHB,file="PDHB.txt",col.names = T,row.names = F,sep = "\t", quote = F)

MTF1=miRNA_mRNA[,8]
colnames(MTF1)[1]="sample"
MTF1=merge(MTF1,miRNA_expr,by="sample")
MTF1=as.data.frame(t(MTF1))
MTF1=tibble::rownames_to_column(MTF1)
colnames(MTF1)=MTF1[1,]
MTF1=MTF1[-1,]
MTF1=merge(MTF1,RNA_cu_clinical,by="sample")
write.table(MTF1,file="MTF1.txt",col.names = T,row.names = F,sep = "\t", quote = F)

GLS=miRNA_mRNA[,9]
colnames(GLS)[1]="sample"
GLS=merge(GLS,miRNA_expr,by="sample")
GLS=as.data.frame(t(GLS))
GLS=tibble::rownames_to_column(GLS)
colnames(GLS)=GLS[1,]
GLS=GLS[-1,]
GLS=merge(GLS,RNA_cu_clinical,by="sample")
write.table(GLS,file="GLS.txt",col.names = T,row.names = F,sep = "\t", quote = F)

CDKN2A=miRNA_mRNA[,10]
colnames(CDKN2A)[1]="sample"
CDKN2A=merge(CDKN2A,miRNA_expr,by="sample")
CDKN2A=as.data.frame(t(CDKN2A))
CDKN2A=tibble::rownames_to_column(CDKN2A)
colnames(CDKN2A)=CDKN2A[1,]
CDKN2A=CDKN2A[-1,]
CDKN2A=merge(CDKN2A,RNA_cu_clinical,by="sample")
write.table(CDKN2A,file="CDKN2A.txt",col.names = T,row.names = F,sep = "\t", quote = F)

LIAS_original=LIAS
########FDX1

LIAS=FDX1
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])

colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "LIPT1_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)


##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,278]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(277:310)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

FDX1_caner=LIAS_caner

write.table(FDX1_caner,file = "FDX1_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)
###
#######
###LIAS 
##ACC
LIAS=LIAS_original
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "LIAS_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,79]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(78:111)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

write.table(LIAS_caner,file = "LIAS_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)
LIAS_caner_original=LIAS_caner

#######
###LIPT1 
LIAS_original=LIAS
LIAS=LIPT1
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "LIPT1_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)


##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,47]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(46:79)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

LIPT1_caner=LIAS_caner

write.table(LIPT1_caner,file = "LIPT1_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#######
###DLD 

LIAS=DLD
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "DLD_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,334]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(333:366)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

DLD_caner=LIAS_caner

write.table(DLD_caner,file = "DLD_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#######
###DLD 

LIAS=DLAT
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "DLAT_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,158]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(157:190)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

DLAT_caner=LIAS_caner

write.table(DLAT_caner,file = "DLAT_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)



#######
###PDHA1 

LIAS=PDHA1
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "PDHA1_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,208]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(207:240)]

y <- as.numeric(exprSet[,"PDHA1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

PDHA1_caner=LIAS_caner

write.table(PDHA1_caner,file = "PDHA1_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#######
###MTF1

LIAS=MTF1
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "MTF1_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,133]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(132:165)]

y <- as.numeric(exprSet[,"MTF1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

MTF1_caner=LIAS_caner

write.table(MTF1_caner,file = "MTF1_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#######
###PDHB

LIAS=PDHB
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "PDHB_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,276]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(275:308)]

y <- as.numeric(exprSet[,"PDHB"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

PDHB_caner=LIAS_caner

write.table(PDHB_caner,file = "PDHB_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#######
###GLS

LIAS=GLS
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "GLS_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,273]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(272:305)]

y <- as.numeric(exprSet[,"GLS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

GLS_caner=LIAS_caner

write.table(GLS_caner,file = "GLS_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#######

LIAS=CDKN2A
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ACC_LIAS=cor_data_df
ACC_LIAS$type="ACC"


##BLCA
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BLCA_LIAS=cor_data_df
BLCA_LIAS$type="BLCA"

##BRCA
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

BRCA_LIAS=cor_data_df
BRCA_LIAS$type="BRCA"

##CESC
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CESC_LIAS=cor_data_df
CESC_LIAS$type="CESC"

##CHOL
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

CHOL_LIAS=cor_data_df
CHOL_LIAS$type="CHOL"

##COAD
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

COAD_LIAS=cor_data_df
COAD_LIAS$type="COAD"

##DLBC
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

DLBC_LIAS=cor_data_df
DLBC_LIAS$type="DLBC"

##ESCA
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

ESCA_LIAS=cor_data_df
ESCA_LIAS$type="ESCA"

##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_LIAS=cor_data_df
HNSC_LIAS$type="HNSC"

##KICH
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KICH_LIAS=cor_data_df
KICH_LIAS$type="KICH"

##KIRC
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_LIAS=cor_data_df
KIRC_LIAS$type="KIRC"

##KIRP
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRP_LIAS=cor_data_df
KIRP_LIAS$type="KIRP"

##LAML
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LAML_LIAS=cor_data_df
LAML_LIAS$type="LAML"

##LGG
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIAS=cor_data_df
LGG_LIAS$type="LGG"

##LIHC
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_LIAS=cor_data_df
LIHC_LIAS$type="LIHC"

##LUAD
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_LIAS=cor_data_df
LUAD_LIAS$type="LUAD"

##LUSC
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUSC_LIAS=cor_data_df
LUSC_LIAS$type="LUSC"

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

MESO_LIAS=cor_data_df
MESO_LIAS$type="MESO"

##OV
C_LIAS=LIAS[grep(pattern="OV",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

OV_LIAS=cor_data_df
OV_LIAS$type="OV"

##PAAD
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PAAD_LIAS=cor_data_df
PAAD_LIAS$type="PAAD"

##PCPG
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PCPG_LIAS=cor_data_df
PCPG_LIAS$type="PCPG"

##PRAD
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_LIAS=cor_data_df
PRAD_LIAS$type="PRAD"

##READ
C_LIAS=LIAS[grep(pattern="READ",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

READ_LIAS=cor_data_df
READ_LIAS$type="READ"

##SARC
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SARC_LIAS=cor_data_df
SARC_LIAS$type="SARC"

##SKCM
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

SKCM_LIAS=cor_data_df
SKCM_LIAS$type="SKCM"
write.table(SKCM_LIAS,file = "CDKN2A_SKCM.txt", sep = "\t", quote = F, col.names = T, row.names = F)

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

STAD_LIAS=cor_data_df
STAD_LIAS$type="STAD"

##TGCT
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

TGCT_LIAS=cor_data_df
TGCT_LIAS$type="TGCT"

##THCA
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df
THCA_LIAS$type="THCA"

##THYM
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THYM_LIAS=cor_data_df
THYM_LIAS$type="THYM"

##UCEC
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_LIAS=cor_data_df
UCEC_LIAS$type="UCEC"

##UCS
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCS_LIAS=cor_data_df
UCS_LIAS$type="UCS"

##UVM
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,83]),] 
exprSet=C_LIAS
exprSet=exprSet[,-1]
exprSet=exprSet[,-c(82:115)]

y <- as.numeric(exprSet[,"CDKN2A"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UVM_LIAS=cor_data_df
UVM_LIAS$type="UVM"

##
LIAS_caner=rbind(ACC_LIAS,BLCA_LIAS)
LIAS_caner=rbind(LIAS_caner,BRCA_LIAS)
LIAS_caner=rbind(LIAS_caner,CESC_LIAS)
LIAS_caner=rbind(LIAS_caner,CHOL_LIAS)
LIAS_caner=rbind(LIAS_caner,COAD_LIAS)
LIAS_caner=rbind(LIAS_caner,DLBC_LIAS)
LIAS_caner=rbind(LIAS_caner,ESCA_LIAS)
LIAS_caner=rbind(LIAS_caner,HNSC_LIAS)
LIAS_caner=rbind(LIAS_caner,KICH_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRC_LIAS)
LIAS_caner=rbind(LIAS_caner,KIRP_LIAS)
LIAS_caner=rbind(LIAS_caner,LAML_LIAS)
LIAS_caner=rbind(LIAS_caner,LGG_LIAS)
LIAS_caner=rbind(LIAS_caner,LIHC_LIAS)
LIAS_caner=rbind(LIAS_caner,LUAD_LIAS)
LIAS_caner=rbind(LIAS_caner,LUSC_LIAS)
LIAS_caner=rbind(LIAS_caner,MESO_LIAS)
LIAS_caner=rbind(LIAS_caner,OV_LIAS)
LIAS_caner=rbind(LIAS_caner,PAAD_LIAS)
LIAS_caner=rbind(LIAS_caner,PCPG_LIAS)
LIAS_caner=rbind(LIAS_caner,PRAD_LIAS)
LIAS_caner=rbind(LIAS_caner,READ_LIAS)
LIAS_caner=rbind(LIAS_caner,SARC_LIAS)
LIAS_caner=rbind(LIAS_caner,SKCM_LIAS)
LIAS_caner=rbind(LIAS_caner,STAD_LIAS)
LIAS_caner=rbind(LIAS_caner,TGCT_LIAS)
LIAS_caner=rbind(LIAS_caner,THCA_LIAS)
LIAS_caner=rbind(LIAS_caner,THYM_LIAS)
LIAS_caner=rbind(LIAS_caner,UCEC_LIAS)
LIAS_caner=rbind(LIAS_caner,UCS_LIAS)
LIAS_caner=rbind(LIAS_caner,UVM_LIAS)

CDKN2A_caner=LIAS_caner

write.table(CDKN2A_caner,file = "CDKN2A_cancer.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#####

