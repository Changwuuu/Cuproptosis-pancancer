library(openxlsx)
T_N<- read.xlsx("T_N.xlsx",sheet = "Sheet1")
Methy_cu<- read.xlsx("Methy_cu.xlsx",sheet = "Sheet3")
Methy_cu_TN=merge(T_N,Methy_cu,by="sample")
write.table(Methy_cu_TN,"Methy_cu_TN.txt",quote = F,sep = "\t",col.names = T,row.names = F)
Methy_cu_TN<-read.xlsx("Methy_cu_TN.xlsx",sheet = "Sheet1")

##### caculate tumor vs normal methylation fold change
indata=t(Methy_cu_TN[,c(1,5,6,7,8,9,10,11,12,13,14)])
indata=as.data.frame(indata)
colnames(indata)=indata[1,]
indata=indata[-1,]

####BLCA
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="BLCA");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="BLCA");
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

####BRCA
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="BRCA");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="BRCA");
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

####CHOL
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="CHOL");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="CHOL");
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

####COAD
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="COAD");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="COAD");
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

####ESCA
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="ESCA");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="ESCA");
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

####HNSC
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="HNSC");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="HNSC");
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

####KIRC
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="KIRC");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="KIRC");
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

####KIRP
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="KIRP");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="KIRP");
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

####LIHC
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="LIHC");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="LIHC");
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

####LUAD
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="LUAD");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="LUAD");
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

####LUSC
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="LUSC");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="LUSC");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

LUSC=cbind(rownames(indata),result)
colnames(LUSC)[1]="SYMBOL"
colnames(LUSC)[2]="LUSC" 

####PAAD
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="PAAD");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="PAAD");
probeFC<-rep(0,nrow(indata))
for(i in 1:nrow(indata))
{
  probeFC[i]<-mean(as.numeric(as.vector(indata[i,sampleIdsCase])))/mean(as.numeric(as.vector(indata[i,sampleIdsControl])));
}
probeFC<-log(probeFC,base=2);
result<-probeFC;

PAAD=cbind(rownames(indata),result)
colnames(PAAD)[1]="SYMBOL"
colnames(PAAD)[2]="PAAD" 

####PRAD
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="PRAD");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="PRAD");
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

####READ
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="READ");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="READ");
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

####THCA
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="THCA");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="THCA");
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


####UCEC
sampleIdsCase<-which(Methy_cu_TN$TN=="tumor"&Methy_cu_TN$cancer=="UCEC");
sampleIdsControl<-which(Methy_cu_TN$TN=="normal"&Methy_cu_TN$cancer=="UCEC");
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

Methy_TN_CU=cbind(BLCA,BRCA,CHOL,COAD,ESCA,HNSC,KIRC,KIRP,LIHC,LUAD,LUSC,PAAD,PRAD,READ,THCA,UCEC)

write.table(Methy_TN_CU, file = "Methy_TN_CU_FC.txt", sep = "\t", row.names = F,col.names = T, quote = F)
#########combine methylation and mRNA expression

RNA_cu=read.table(file = "RNA_cu.txt", header = T, check.names = F,row.names = 1)

RNA_cu=tibble::rownames_to_column(RNA_cu)
colnames(RNA_cu)[1]="sample"
RNA_cu=RNA_cu[,-12]
cu_RNASeq_methy=merge(RNA_cu,Methy_cu,by="sample")
cu_RNASeq_methy=merge(cu_RNASeq_methy,XENA_clinical,by="sample")

write.table(cu_RNASeq_methy, file = "cu_RNASeq_methy.txt", sep = "\t", row.names = F,col.names = T, quote = F)
####caculating the correlation of methylation and mRNA expression in tumor by Sangerbox. Please cite :Shen, et al. 2022. Sangerbox: A comprehensive, interaction-friendly clinical bioinformatics analysis platform. iMeta 1(3): e36. https://doi.org/10.1002/imt2.36
#####combine methylation and clinical information. Using Sangerbox to caculate the survival value.Please cite :Shen, et al. 2022. Sangerbox: A comprehensive, interaction-friendly clinical bioinformatics analysis platform. iMeta 1(3): e36. https://doi.org/10.1002/imt2.36

methy_cu_clinical=merge(Methy_cu,XENA_clinical,by="sample")

write.table(methy_cu_clinical, file = "methy_cu_clinical.txt", sep = "\t", row.names = F,col.names = T, quote = F)


######combine mRNA and clinical information
RNA_cu_clinical=merge(RNA_cu,XENA_clinical,by="sample")
write.table(RNA_cu_clinical, file = "RNA_cu_clinical.txt", sep = "\t", row.names = F,col.names = T, quote = F)


