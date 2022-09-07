######
library(openxlsx)
######combine cuproptosis and immune cell inflitration
RNASeq2_gsva_cu<- read.xlsx("RNASeq2_gsva_cu.xlsx",sheet = "Sheet1")
TCGA_cibersort<- read.xlsx("TCGA_cibersort.xlsx",sheet = "Sheet1")

TCGA_cibersort=TCGA_cibersort[!duplicated(TCGA_cibersort$sample),]
RNA_Cuscore=merge(RNASeq2_gsva_cu,TCGA_cibersort,by="sample")
write.table(RNA_Cuscore,file = "RNA_Cuscore.txt", sep = "\t", quote = F, col.names = T, row.names = F)

######cuproptosis, immune cell inflitration, CORRELATION
LIAS=RNA_Cuscore
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
##GBM
C_LIAS=LIAS[grep(pattern="GBM",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

GBM_LIAS=cor_data_df
GBM_LIAS$type="GBM"


##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##MESO
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="OV",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="READ",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,2]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(5:27)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
LIAS_caner=rbind(LIAS_caner,GBM_LIAS)
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

CIBSORT_cu_Cor=LIAS_caner

write.table(CIBSORT_cu_Cor,file = "CIBSORT_cu_Cor.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#######HRD cuproptosis correlation
TCGA.HRD=read.table(file = "TCGA.HRD.txt",sep = "\t",header = T,check.names = F)

TCGA.HRD=as.data.frame(t(TCGA.HRD))
colnames(TCGA.HRD)=TCGA.HRD[1,]
TCGA.HRD=TCGA.HRD[-1,]
colnames(TCGA.HRD)[1]="sample"
colnames(mutation.load_updated)[3]="sample"

HRD_cuscore=merge(RNASeq2_gsva_cu,TCGA.HRD,by="sample")
HRD_cuscore=merge(HRD_cuscore,mutation.load_updated,by="sample")
write.table(HRD_cuscore,file = "HRD_cuscore.txt", sep = "\t", quote = F, col.names = T, row.names = F)
####
LIAS=HRD_cuscore
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
##GBM
C_LIAS=LIAS[grep(pattern="GBM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

GBM_LIAS=cor_data_df
GBM_LIAS$type="GBM"


##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##MESO
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="OV",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="READ",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]
y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(4,5,6,7,8,12)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
LIAS_caner=rbind(LIAS_caner,GBM_LIAS)
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

HRD_cu_Cor=LIAS_caner

write.table(HRD_cu_Cor,file = "HRD_cu_Cor.txt", sep = "\t", quote = F, col.names = T, row.names = F)

########immune genes
library(data.table)
RNASeq=fread("RNASeqV2.txt", data.table=FALSE)
RNASeq=RNASeq[!duplicated(RNASeq$sample),]
rownames(RNASeq)=RNASeq$sample
immunegenes=read.table(file = "immunegenes.txt",sep = "\t",header = T,check.names = F)

immune_expr=merge(immunegenes,RNASeq,by="sample")
immune_expr=as.data.frame(t(immune_expr))
immune_expr=tibble::rownames_to_column(immune_expr)
colnames(immune_expr)=immune_expr[1,]
immune_expr=immune_expr[-1,]
immune_expr_cuscore=merge(immune_expr,RNASeq2_gsva_cu,by="sample")
immune_expr_cuscore=immune_expr_cuscore[!duplicated(immune_expr_cuscore$sample),]
immune_expr_cuscore=merge(immune_expr_cuscore,XENA_clinical,by="sample")

write.table(immune_expr_cuscore,file = "immune_expr_cuscore.txt", sep = "\t", quote = F, col.names = T, row.names = F)

#####


######CUSCORE immunegenes CORRELATION
LIAS=immune_expr_cuscore
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
##GBM
C_LIAS=LIAS[grep(pattern="GBM",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

GBM_LIAS=cor_data_df
GBM_LIAS$type="GBM"


##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##MESO
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="OV",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="READ",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,27]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2:26)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
LIAS_caner=rbind(LIAS_caner,GBM_LIAS)
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

immunegenes_cu_Cor=LIAS_caner

write.table(immunegenes_cu_Cor,file = "immunegenes_cu_Cor.txt", sep = "\t", quote = F, col.names = T, row.names = F)


#####estimatescore 

library(IOBR)
library(estimate) 
library(tidyverse)
library(tidyHeatmap)

exprset=fread("RNASeq2.txt", data.table=FALSE)
rownames(exprset)=exprset[,1]
exprset=exprset[,-1]
estimate_TCGA<-deconvo_tme(eset = exprset, method = "estimate", arrays = FALSE, perm = 200 )
write.table(estimate_TCGA,"estimate_TCGA.txt",quote = F,sep = "\t",col.names = T,row.names = F)

RNASeq.estimatescore=read.xlsx("RNASeq.estimatescore.xlsx",sheet = "Sheet1")

estimatescore=merge(RNASeq.estimatescore,RNASeq2_gsva_cu,by="sample")
estimatescore=merge(estimatescore,XENA_clinical,by="sample")


######CUSCORE estimatescore CORRELATION
LIAS=estimatescore
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
##GBM
C_LIAS=LIAS[grep(pattern="GBM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

GBM_LIAS=cor_data_df
GBM_LIAS$type="GBM"


##HNSC
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##MESO
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="OV",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="READ",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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

##STAD
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]


y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,9]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(2,3,4,7)]

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
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
LIAS_caner=rbind(LIAS_caner,GBM_LIAS)
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

estimatescore_cu_Cor=LIAS_caner

write.table(estimatescore_cu_Cor,file = "estimatescore_cu_Cor.txt", sep = "\t", quote = F, col.names = T, row.names = F)

