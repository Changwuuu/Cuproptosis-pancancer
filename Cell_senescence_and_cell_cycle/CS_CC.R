library(data.table)
RNASeq=fread("RNASeq2.txt", data.table=FALSE)

library(openxlsx)
cs <- read.xlsx("cs_genes.xlsx",sheet = "Sheet1")
cc <- read.xlsx("cc_genes.xlsx",sheet = "Sheet1")

RNASeq2_gsva_cu=read.xlsx("RNASeq2_gsva_cu.xlsx",sheet = "Sheet1")
XENA_clinical=read.xlsx("XENA_clinical.xlsx",sheet = "Sheet1")
######
colnames(RNASeq)[1]="Symbol"
cs_TCGA=merge(RNASeq,cs,by="Symbol")
cc_TCGA=merge(RNASeq,cc,by="Symbol")

cs_TCGA=as.data.frame(t(cs_TCGA))
cs_TCGA=tibble::rownames_to_column(cs_TCGA)
colnames(cs_TCGA)=cs_TCGA[1,]
cs_TCGA=cs_TCGA[-1,]
colnames(cs_TCGA)[1]="sample"
write.table(cs_TCGA,"cs_TCGA.txt",quote = F,sep = "\t",col.names = T,row.names = F)
cs_TCGA=read.xlsx("cs_TCGA.xlsx",sheet = "Sheet1")


cc_TCGA=as.data.frame(t(cc_TCGA))
cc_TCGA=tibble::rownames_to_column(cc_TCGA)
colnames(cc_TCGA)=cc_TCGA[1,]
cc_TCGA=cc_TCGA[-1,]
colnames(cc_TCGA)[1]="sample"
write.table(cc_TCGA,"cc_TCGA.txt",quote = F,sep = "\t",col.names = T,row.names = F)
cc_TCGA=read.xlsx("cc_TCGA.xlsx",sheet = "Sheet1")


######CUSCORE CS_genes CORRELATION
cs_cu=merge(RNASeq2_gsva_cu,cs_TCGA,by="sample")
cs_cu=merge(XENA_clinical,cs_cu,by="sample")


LIAS=cs_cu
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]

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
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="GBM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]

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
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="OV",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="READ",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]

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
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:139)]


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

cs_cu_Cor=LIAS_caner

write.table(cs_cu_Cor,file = "cs_cu_Cor.txt", sep = "\t", quote = F, col.names = T, row.names = F)



######CUSCORE CC_genes CORRELATION
cc_cu=merge(RNASeq2_gsva_cu,cc_TCGA,by="sample")
cc_cu=merge(XENA_clinical,cc_cu,by="sample")


LIAS=cc_cu
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]

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
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="GBM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]

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
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="OV",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="READ",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]

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
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,3]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(37:155)]


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

cc_cu_Cor=LIAS_caner

write.table(cc_cu_Cor,file = "cc_cu_Cor.txt", sep = "\t", quote = F, col.names = T, row.names = F)
