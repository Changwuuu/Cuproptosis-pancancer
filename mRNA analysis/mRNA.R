#####RNASeqV2 can be download in Xena(http://xena.ucsc.edu/).
library(data.table)
RNASeq=fread("RNASeqV2.txt", data.table=FALSE)
RNASeq=RNASeq[!duplicated(RNASeq$sample),]
rownames(RNASeq)=RNASeq$sample
cu_RNASeq=RNASeq[c("FDX1","LIAS","LIPT1","DLD","DLAT","PDHA1","PDHB","MTF1",
          "GLS","CDKN2A"),]
write.table(cu_RNASeq,"cu_RNASeq.txt",quote = F,sep = "\t",col.names = T,row.names = F)

na.omit(RNASeq)[1:4,1:4]
RNASeq2=na.omit(RNASeq)
write.table(RNASeq,"RNASeq.txt",quote = F,sep = "\t",col.names = T,row.names = F)
write.table(RNASeq2,"RNASeq2.txt",quote = F,sep = "\t",col.names = T,row.names = F)

RNA_cu=t(cu_RNASeq[,-1])
RNA_cu=as.data.frame(RNA_cu)

RNA_cu$sample=rownames(RNA_cu)                              


library(openxlsx)
######combine CNV and mRNA expression
RNA_cu <- read.xlsx("RNA_cu.xlsx",sheet = "Sheet1")
CNV_cu=read.table(file = "CNV_cu.txt",sep = "\t",header = T,check.names = F)
XENA_clinical=read.table(file = "XENA_clinical.txt",sep = "\t",header = T,check.names = F)

colnames(CNV_cu)[1]="sample"

RNA_CNV=merge(RNA_cu,CNV_cu,by="sample")
RNA_CNV_clinical=merge(RNA_CNV,XENA_clinical,by="sample")
write.table(RNA_CNV_clinical,"RNA_CNV_clinical.txt",quote = F,sep = "\t",col.names = T,row.names = F)#####Using "RNA_CNV_clinical" to plot box plots of CNV/mRNA in sangerbox. Please cite:Shen, et al. 2022. Sangerbox: A comprehensive, interaction-friendly clinical bioinformatics analysis platform. iMeta 1(3): e36. https://doi.org/10.1002/imt2.36


#####combine cuprotosis score, hallmark score and clinical information
RNA_cu_clinical<- read.xlsx("RNA_cu_clinical.xlsx",sheet = "Sheet1")

RNA_cu_clinical=merge(RNA_cu_clinical,RNASeq2_gsva_cu,by="sample")
RNA_cu_clinical=merge(RNA_cu_clinical,RNASeq2_gsva_hallmark,by="sample")

write.table(RNA_cu_clinical,"RNA_cu_clinical_cuscore.txt",quote = F,sep = "\t",col.names = T,row.names = F)

###### tumor vs normal, mRNA expression fold change, using sangerbox to caculated significance
T_N<- read.xlsx("T_N_cu.xlsx",sheet = "Sheet1")
indata=t(T_N[,c(1,5,6,7,8,9,10,11,12,13,14)])
indata=as.data.frame(indata)
colnames(indata)=indata[1,]
indata=indata[-1,]

##BLCA
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="BLCA");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="BLCA");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

BLCA=cbind(rownames(indata),result)
colnames(BLCA)[1]="SYMBOL"
colnames(BLCA)[2]="BLCA"

##BRCA
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="BRCA");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="BRCA");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

BRCA=cbind(rownames(indata),result)
colnames(BRCA)[1]="SYMBOL"
colnames(BRCA)[2]="BRCA"

##CHOL
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="CHOL");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="CHOL");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

CHOL=cbind(rownames(indata),result)
colnames(CHOL)[1]="SYMBOL"
colnames(CHOL)[2]="CHOL"

##COAD
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="COAD");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="COAD");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

COAD=cbind(rownames(indata),result)
colnames(COAD)[1]="SYMBOL"
colnames(COAD)[2]="COAD"

##ESCA
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="ESCA");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="ESCA");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

ESCA=cbind(rownames(indata),result)
colnames(ESCA)[1]="SYMBOL"
colnames(ESCA)[2]="ESCA"

##HNSC
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="HNSC");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="HNSC");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

HNSC=cbind(rownames(indata),result)
colnames(HNSC)[1]="SYMBOL"
colnames(HNSC)[2]="HNSC"

##KICH
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="KICH");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="KICH");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

KICH=cbind(rownames(indata),result)
colnames(KICH)[1]="SYMBOL"
colnames(KICH)[2]="KICH"

##KIRC
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="KIRC");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="KIRC");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

KIRC=cbind(rownames(indata),result)
colnames(KIRC)[1]="SYMBOL"
colnames(KIRC)[2]="KIRC"

##KIRP
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="KIRP");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="KIRP");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

KIRP=cbind(rownames(indata),result)
colnames(KIRP)[1]="SYMBOL"
colnames(KIRP)[2]="KIRP"

##LIHC
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="LIHC");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="LIHC");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

LIHC=cbind(rownames(indata),result)
colnames(LIHC)[1]="SYMBOL"
colnames(LIHC)[2]="LIHC"

##LUAD
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="LUAD");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="LUAD");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

LUAD=cbind(rownames(indata),result)
colnames(LUAD)[1]="SYMBOL"
colnames(LUAD)[2]="LUAD"

##LUSC
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="LUSC");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="LUSC");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

UCEC=cbind(rownames(indata),result)
colnames(LUSC)[1]="SYMBOL"
colnames(LUSC)[2]="LUSC"

##PRAD
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="PRAD");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="PRAD");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

PRAD=cbind(rownames(indata),result)
colnames(PRAD)[1]="SYMBOL"
colnames(PRAD)[2]="PRAD"

##READ
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="READ");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="READ");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

READ=cbind(rownames(indata),result)
colnames(READ)[1]="SYMBOL"
colnames(READ)[2]="READ"

##STAD
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="STAD");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="STAD");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

STAD=cbind(rownames(indata),result)
colnames(STAD)[1]="SYMBOL"
colnames(STAD)[2]="STAD"

##THCA
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="THCA");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="THCA");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

THCA=cbind(rownames(indata),result)
colnames(THCA)[1]="SYMBOL"
colnames(THCA)[2]="THCA"

##UCEC
sampleIdsCase<-which(T_N$TN=="tumor"&T_N$cancer=="UCEC");
sampleIdsControl<-which(T_N$TN=="normal"&T_N$cancer=="UCEC");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

UCEC=cbind(rownames(indata),result)
colnames(UCEC)[1]="SYMBOL"
colnames(UCEC)[2]="UCEC"

TN_CU=cbind(BLCA,BRCA,CHOL,COAD,ESCA,HNSC,KICH,KIRC,KIRP,LIHC,LUAD,LUSC,PRAD,READ,STAD,THCA,UCEC)

write.table(TN_CU, file = "TN_CU_FC.txt", sep = "\t", row.names = F,col.names = T, quote = F)


#####correlation between cuprotosis score and hallmark 

LIAS=RNA_cu_clinical
##ACC
C_LIAS=LIAS[grep(pattern="ACC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="BLCA",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="BRCA",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="CESC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="CHOL",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="COAD",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="DLBC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="ESCA",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="GBM",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="HNSC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="KICH",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="KIRC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="KIRP",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="LAML",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="LGG",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="LIHC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="LUAD",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="LUSC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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

##
C_LIAS=LIAS[grep(pattern="MESO",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="OV",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="PAAD",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="PCPG",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="PRAD",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="READ",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="SARC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="SKCM",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="STAD",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="TGCT",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="THCA",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="THYM",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="UCEC",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="UCS",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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
C_LIAS=LIAS[grep(pattern="UVM",LIAS[,12]),] 
exprSet=C_LIAS
exprSet=exprSet[,c(45:95)]


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

hallmark_cu_Cor=LIAS_caner

write.table(hallmark_cu_Cor,file = "hallmark_cu_Cor.txt", sep = "\t", quote = F, col.names = T, row.names = F)

###correlation among cuproptosis scores and regulators
cu_cuscore=RNA_cu_clinical[,c(2,3,4,5,6,7,8,9,10,11,43,44,45)]
library("Hmisc")
res2 <- rcorr(as.matrix(cu_cuscore))
res2
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
res2<-rcorr(as.matrix(MAPK))
flattenCorrMatrix(res2$r, res2$P)

library(corrplot)


corrplot(res2$r, type="lower", tl.col="black", tl.srt = 90, tl.pos = "lt",
         p.mat = res2$P, sig.level = 0.05, insig = "blank",tl.cex = 0.92,cl.cex = 1.1, cl.align.text = "l")
corrplot(res2$r, add = TRUE, type = "upper", method = "number",  tl.col="black" , tl.pos = "n",
         diag = FALSE,p.mat = res2$P, sig.level = 0.05, insig = "blank",tl.cex = 1.0,cl.pos = "n",number.cex=0.8)



