library(data.table)
RNASeq=fread("RNASeqV2.txt", data.table=FALSE)
RNASeq=RNASeq[!duplicated(RNASeq$sample),]

CU=read.table(file = "CU.TXT",header = T, sep = "\t")

##
BLCA=read.table(file = "BLCA.TXT",header = T, sep = "\t")
BLCA=merge(BLCA,CU,by="Target")
BLCA=BLCA[,1:2]
colnames(BLCA)[2]="sample"
BLCA=merge(BLCA,RNASeq,by="sample")
BLCA_target=BLCA[,1:2]
BLCA_target$type="BLCA"

BRCA=read.table(file = "BRCA.TXT",header = T, sep = "\t")
BRCA=merge(BRCA,CU,by="Target")
BRCA=BRCA[,1:2]
colnames(BRCA)[2]="sample"
BRCA=merge(BRCA,RNASeq,by="sample")
BRCA_target=BRCA[,1:2]
BRCA_target$type="BRCA"

CESC=read.table(file = "CESC.TXT",header = T, sep = "\t")
CESC=merge(CESC,CU,by="Target")
CESC=CESC[,1:2]
colnames(CESC)[2]="sample"
CESC=merge(CESC,RNASeq,by="sample")
CESC_target=CESC[,1:2]
CESC_target$type="CESC"

HNSC=read.table(file = "HNSC.TXT",header = T, sep = "\t")
HNSC=merge(HNSC,CU,by="Target")
HNSC=HNSC[,1:2]
colnames(HNSC)[2]="sample"
HNSC=merge(HNSC,RNASeq,by="sample")
HNSC_target=HNSC[,1:2]
HNSC_target$type="HNSC"


KIRC=read.table(file = "KIRC.TXT",header = T, sep = "\t")
KIRC=merge(KIRC,CU,by="Target")
KIRC=KIRC[,1:2]
colnames(KIRC)[2]="sample"
KIRC=merge(KIRC,RNASeq,by="sample")
KIRC_target=KIRC[,1:2]
KIRC_target$type="KIRC"

KIRP=read.table(file = "KIRP.TXT",header = T, sep = "\t")
KIRP=merge(KIRP,CU,by="Target")
KIRP=KIRP[,1:2]
colnames(KIRP)[2]="sample"
KIRP=merge(KIRP,RNASeq,by="sample")
KIRP_target=KIRP[,1:2]
KIRP_target$type="KIRP"


LGG=read.table(file = "LGG.TXT",header = T, sep = "\t")
LGG=merge(LGG,CU,by="Target")
LGG=LGG[,1:2]
colnames(LGG)[2]="sample"
LGG=merge(LGG,RNASeq,by="sample")
LGG_target=LGG[,1:2]
LGG_target$type="LGG"


LIHC=read.table(file = "LIHC.TXT",header = T, sep = "\t")
LIHC=merge(LIHC,CU,by="Target")
LIHC=LIHC[,1:2]
colnames(LIHC)[2]="sample"
LIHC=merge(LIHC,RNASeq,by="sample")
LIHC_target=LIHC[,1:2]
LIHC_target$type="LIHC"


LUAD=read.table(file = "LUAD.TXT",header = T, sep = "\t")
LUAD=merge(LUAD,CU,by="Target")
LUAD=LUAD[,1:2]
colnames(LUAD)[2]="sample"
LUAD=merge(LUAD,RNASeq,by="sample")
LUAD_target=LUAD[,1:2]
LUAD_target$type="LUAD"

OV=read.table(file = "OV.TXT",header = T, sep = "\t")
OV=merge(OV,CU,by="Target")
OV=OV[,1:2]
colnames(OV)[2]="sample"
OV=merge(OV,RNASeq,by="sample")
OV_target=OV[,1:2]
OV_target$type="OV"


PRAD=read.table(file = "PRAD.TXT",header = T, sep = "\t")
PRAD=merge(PRAD,CU,by="Target")
PRAD=PRAD[,1:2]
colnames(PRAD)[2]="sample"
PRAD=merge(PRAD,RNASeq,by="sample")
PRAD_target=PRAD[,1:2]
PRAD_target$type="PRAD"


SKCM=read.table(file = "SKCM.TXT",header = T, sep = "\t")
SKCM=merge(SKCM,CU,by="Target")
SKCM=SKCM[,1:2]
colnames(SKCM)[2]="sample"
SKCM=merge(SKCM,RNASeq,by="sample")
SKCM_target=SKCM[,1:2]
SKCM_target$type="SKCM"

THCA=read.table(file = "THCA.TXT",header = T, sep = "\t")
THCA=merge(THCA,CU,by="Target")
THCA=THCA[,1:2]
colnames(THCA)[2]="sample"
THCA=merge(THCA,RNASeq,by="sample")
THCA_target=THCA[,1:2]
THCA_target$type="THCA"

UCEC=read.table(file = "UCEC.TXT",header = T, sep = "\t")
UCEC=merge(UCEC,CU,by="Target")
UCEC=UCEC[,1:2]
colnames(UCEC)[2]="sample"
UCEC=merge(UCEC,RNASeq,by="sample")
UCEC_target=UCEC[,1:2]
UCEC_target$type="UCEC"


LNC_cu=rbind(BLCA_target,BRCA_target)
LNC_cu=rbind(LNC_cu,CESC_target)
LNC_cu=rbind(LNC_cu,HNSC_target)
LNC_cu=rbind(LNC_cu,KIRC_target)
LNC_cu=rbind(LNC_cu,KIRP_target)
LNC_cu=rbind(LNC_cu,LGG_target)
LNC_cu=rbind(LNC_cu,LIHC_target)
LNC_cu=rbind(LNC_cu,LUAD_target)
LNC_cu=rbind(LNC_cu,OV_target)
LNC_cu=rbind(LNC_cu,PRAD_target)
LNC_cu=rbind(LNC_cu,SKCM_target)
LNC_cu=rbind(LNC_cu,THCA_target)
LNC_cu=rbind(LNC_cu,UCEC_target)

write.table(LNC_cu,file="LNC_cu.txt",col.names = T,row.names = F,sep = "\t", quote = F)

####
LNC_CU=read.table(file = "LNC_CU.TXT",header = T, sep = "\t")
library(openxlsx)

RNA_cu_clinical<- read.xlsx("RNA_cu_clinical.xlsx",sheet = "Sheet1")

###
FDX1=as.data.frame(LNC_CU[,4])
colnames(FDX1)[1]="sample"
FDX1=merge(FDX1,RNASeq,by="sample")
FDX1=as.data.frame(t(FDX1))
FDX1=tibble::rownames_to_column(FDX1)
colnames(FDX1)=FDX1[1,]
FDX1=FDX1[-1,]
FDX1=merge(FDX1,RNA_cu_clinical,by="sample")

LIAS=as.data.frame(LNC_CU[,6])
colnames(LIAS)[1]="sample"
LIAS=merge(LIAS,RNASeq,by="sample")
LIAS=as.data.frame(t(LIAS))
LIAS=tibble::rownames_to_column(LIAS)
colnames(LIAS)=LIAS[1,]
LIAS=LIAS[-1,]
LIAS=merge(LIAS,RNA_cu_clinical,by="sample")



LIPT1=as.data.frame(LNC_CU[,7])
colnames(LIPT1)[1]="sample"
LIPT1=merge(LIPT1,RNASeq,by="sample")
LIPT1=as.data.frame(t(LIPT1))
LIPT1=tibble::rownames_to_column(LIPT1)
colnames(LIPT1)=LIPT1[1,]
LIPT1=LIPT1[-1,]
LIPT1=merge(LIPT1,RNA_cu_clinical,by="sample")

DLD=as.data.frame(LNC_CU[,3])
colnames(DLD)[1]="sample"
DLD=merge(DLD,RNASeq,by="sample")
DLD=as.data.frame(t(DLD))
DLD=tibble::rownames_to_column(DLD)
colnames(DLD)=DLD[1,]
DLD=DLD[-1,]
DLD=merge(DLD,RNA_cu_clinical,by="sample")

DLAT=as.data.frame(LNC_CU[,2])
colnames(DLAT)[1]="sample"
DLAT=merge(DLAT,RNASeq,by="sample")
DLAT=as.data.frame(t(DLAT))
DLAT=tibble::rownames_to_column(DLAT)
colnames(DLAT)=DLAT[1,]
DLAT=DLAT[-1,]
DLAT=merge(DLAT,RNA_cu_clinical,by="sample")

PDHA1=as.data.frame(LNC_CU[,9])
colnames(PDHA1)[1]="sample"
PDHA1=merge(PDHA1,RNASeq,by="sample")
PDHA1=as.data.frame(t(PDHA1))
PDHA1=tibble::rownames_to_column(PDHA1)
colnames(PDHA1)=PDHA1[1,]
PDHA1=PDHA1[-1,]
PDHA1=merge(PDHA1,RNA_cu_clinical,by="sample")

PDHB=as.data.frame(LNC_CU[,10])
colnames(PDHB)[1]="sample"
PDHB=merge(PDHB,RNASeq,by="sample")
PDHB=as.data.frame(t(PDHB))
PDHB=tibble::rownames_to_column(PDHB)
colnames(PDHB)=PDHB[1,]
PDHB=PDHB[-1,]
PDHB=merge(PDHB,RNA_cu_clinical,by="sample")

MTF1=as.data.frame(LNC_CU[,8])
colnames(MTF1)[1]="sample"
MTF1=merge(MTF1,RNASeq,by="sample")
MTF1=as.data.frame(t(MTF1))
MTF1=tibble::rownames_to_column(MTF1)
colnames(MTF1)=MTF1[1,]
MTF1=MTF1[-1,]
MTF1=merge(MTF1,RNA_cu_clinical,by="sample")

GLS=as.data.frame(LNC_CU[,5])
colnames(GLS)[1]="sample"
GLS=merge(GLS,RNASeq,by="sample")
GLS=as.data.frame(t(GLS))
GLS=tibble::rownames_to_column(GLS)
colnames(GLS)=GLS[1,]
GLS=GLS[-1,]
GLS=merge(GLS,RNA_cu_clinical,by="sample")

CDKN2A=as.data.frame(LNC_CU[,1])
colnames(CDKN2A)[1]="sample"
CDKN2A=merge(CDKN2A,RNASeq,by="sample")
CDKN2A=as.data.frame(t(CDKN2A))
CDKN2A=tibble::rownames_to_column(CDKN2A)
colnames(CDKN2A)=CDKN2A[1,]
CDKN2A=CDKN2A[-1,]
CDKN2A=merge(CDKN2A,RNA_cu_clinical,by="sample")

####
C_FDX1=FDX1[grep(pattern="KIRC",FDX1[,20]),] 
exprSet=C_FDX1

exprSet=exprSet[,c(3:10)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

KIRC_FDX1=cor_data_df
KIRC_FDX1$type="KIRC"

#####
C_FDX1=FDX1[grep(pattern="LIHC",FDX1[,20]),] 
exprSet=C_FDX1

exprSet=exprSet[,c(3:10)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LIHC_FDX1=cor_data_df
LIHC_FDX1$type="LIHC"

#####
C_FDX1=FDX1[grep(pattern="LUAD",FDX1[,20]),] 
exprSet=C_FDX1

exprSet=exprSet[,c(3:10)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LUAD_FDX1=cor_data_df
LUAD_FDX1$type="LIHC"
#####
C_FDX1=FDX1[grep(pattern="PRAD",FDX1[,20]),] 
exprSet=C_FDX1

exprSet=exprSet[,c(3:10)]

y <- as.numeric(exprSet[,"FDX1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

PRAD_FDX1=cor_data_df
PRAD_FDX1$type="PRAD"

FDX1_LNC=rbind(KIRC_FDX1,LIHC_FDX1,LUAD_FDX1,PRAD_FDX1)

write.table(FDX1_LNC,file="FDX1_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)

#####Remember Change the cancer type everytime
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,28]),] #####Remember Change the cancer type everytime
exprSet=C_LIAS

exprSet=exprSet[,c(2:19)]

y <- as.numeric(exprSet[,"LIAS"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

THCA_LIAS=cor_data_df  #####Remember Change the cancer type everytime
THCA_LIAS$type="THCA"  #####Remember Change the cancer type everytime

LIAS_LNC=rbind(BLCA_LIAS,KIRC_LIAS,KIRP_LIAS,LGG_LIAS,LUAD_LIAS,PRAD_LIAS,SKCM_LIAS,THCA_LIAS,UCEC_LIAS)
write.table(LIAS_LNC,file="LIAS_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)


#########Remember Change the cancer type everytime

C_LIPT1=LIPT1[grep(pattern="LGG",LIPT1[,27]),] #####Remember Change the cancer type everytime
exprSet=C_LIPT1

exprSet=exprSet[,c(2:19)]

y <- as.numeric(exprSet[,"LIPT1"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

LGG_LIPT1=cor_data_df#####Remember Change the cancer type everytime
LGG_LIPT1$type="LGG"#####Remember Change the cancer type everytime

LIPT1_LNC=rbind(BLCA_LIPT1,BRCA_LIPT1,LGG_LIPT1,LUAD_LIPT1,PRAD_LIPT1,SKCM_LIPT1,UCEC_LIPT1)
write.table(LIPT1_LNC,file="LIPT1_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)


########Remember Change the cancer type everytime
C_DLD=DLD[grep(pattern="UCEC",DLD[,28]),] #####Remember Change the cancer type everytime
exprSet=C_DLD

exprSet=exprSet[,c(2:21)]

y <- as.numeric(exprSet[,"DLD"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_DLD=cor_data_df#####Remember Change the cancer type everytime
UCEC_DLD$type="UCEC"#####Remember Change the cancer type everytime

DLD_LNC=rbind(BRCA_DLD,KIRC_DLD,KIRP_DLD,PRAD_DLD,THCA_DLD,UCEC_DLD)
write.table(DLD_LNC,file="DLD_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)


##########Remember Change the cancer type everytime
C_DLAT=DLAT[grep(pattern="UCEC",DLAT[,37]),] #####Remember Change the cancer type everytime
exprSet=C_DLAT

exprSet=exprSet[,c(2:31)]

y <- as.numeric(exprSet[,"DLAT"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_DLAT=cor_data_df#####Remember Change the cancer type everytime
KIRC_DLAT$type="KIRC"#####Remember Change the cancer type everytime

DLAT_LNC=rbind(BLCA_DLAT,BRCA_DLAT,HNSC_DLAT,KIRC_DLAT,LGG_DLAT,PRAD_DLAT,THCA_DLAT,UCEC_DLAT)
write.table(DLAT_LNC,file="DLAT_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)


##########Remember Change the cancer type everytime
C_PDHA1=PDHA1[grep(pattern="HNSC",PDHA1[,42]),] #####Remember Change the cancer type everytime
exprSet=C_PDHA1

exprSet=exprSet[,c(2:37)]

colnames <- colnames(exprSet)
y <- as.numeric(exprSet[,"PDHA1"])
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

HNSC_PDHA1=cor_data_df#####Remember Change the cancer type everytime
HNSC_PDHA1$type="HNSC"#####Remember Change the cancer type everytime

PDHA1_LNC=rbind(BLCA_PDHA1,BRCA_PDHA1,KIRC_PDHA1,KIRP_PDHA1,LGG_PDHA1,LIHC_PDHA1,LUAD_PDHA1,OV_PDHA1,PRAD_PDHA1,SKCM_PDHA1,THCA_PDHA1)
write.table(PDHA1_LNC,file="PDHA1_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)

##########Remember Change the cancer type everytime
C_PDHB=PDHB[grep(pattern="UCEC",PDHB[,29]),] #####Remember Change the cancer type everytime
exprSet=C_PDHB

exprSet=exprSet[,c(2:25)]

colnames <- colnames(exprSet)
y <- as.numeric(exprSet[,"PDHB"])
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_PDHB=cor_data_df#####Remember Change the cancer type everytime
UCEC_PDHB$type="UCEC"#####Remember Change the cancer type everytime

PDHB_LNC=rbind(BLCA_PDHB,BRCA_PDHB,HNSC_PDHB,KIRC_PDHB,LUAD_PDHB,OV_PDHB,PRAD_PDHB,THCA_PDHB,UCEC_PDHB)
write.table(PDHB_LNC,file="PDHB_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)


##########Remember Change the cancer type everytime
C_MTF1=MTF1[grep(pattern="UCEC",MTF1[,187]),] #####Remember Change the cancer type everytime
exprSet=C_MTF1

exprSet=exprSet[,c(2:184)]

colnames <- colnames(exprSet)

y <- as.numeric(exprSet[,"MTF1"])
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_MTF1=cor_data_df#####Remember Change the cancer type everytime
UCEC_MTF1$type="UCEC"#####Remember Change the cancer type everytime

MTF1_LNC=rbind(BLCA_MTF1,BRCA_MTF1,CESC_MTF1,HNSC_MTF1,KIRC_MTF1,KIRP_MTF1,LGG_MTF1,LIHC_MTF1,LUAD_MTF1,OV_MTF1,PRAD_MTF1,SKCM_MTF1,THCA_MTF1,UCEC_MTF1)
write.table(MTF1_LNC,file="MTF1_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)


##########Remember Change the cancer type everytime
C_GLS=GLS[grep(pattern="UCEC",GLS[,113]),] #####Remember Change the cancer type everytime
exprSet=C_GLS

exprSet=exprSet[,c(2:111)]

colnames <- colnames(exprSet)

y <- as.numeric(exprSet[,"GLS"])
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_GLS=cor_data_df#####Remember Change the cancer type everytime
UCEC_GLS$type="UCEC"#####Remember Change the cancer type everytime

GLS_LNC=rbind(BLCA_GLS,BRCA_GLS,CESC_GLS,HNSC_GLS,KIRC_GLS,KIRP_GLS,LGG_GLS,LIHC_GLS,LUAD_GLS,OV_GLS,PRAD_GLS,SKCM_GLS,THCA_GLS,UCEC_GLS)
write.table(GLS_LNC,file="GLS_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)

##########Remember Change the cancer type everytime
C_CDKN2A=CDKN2A[grep(pattern="UCEC",CDKN2A[,50]),] #####Remember Change the cancer type everytime
exprSet=C_CDKN2A

exprSet=exprSet[,c(2:49)]

colnames <- colnames(exprSet)

y <- as.numeric(exprSet[,"CDKN2A"])
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")

UCEC_CDKN2A=cor_data_df#####Remember Change the cancer type everytime
UCEC_CDKN2A$type="UCEC"#####Remember Change the cancer type everytime

CDKN2A_LNC=rbind(BLCA_CDKN2A,BRCA_CDKN2A,CESC_CDKN2A,HNSC_CDKN2A,KIRC_CDKN2A,LGG_CDKN2A,LIHC_CDKN2A,PRAD_CDKN2A,THCA_CDKN2A,UCEC_CDKN2A)
write.table(CDKN2A_LNC,file="CDKN2A_LNC.txt",col.names = T,row.names = F,sep = "\t", quote = F)


