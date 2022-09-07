
rm(list = ls()) 
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='./DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

testExpr<- read.table(file = "ACC_expr.txt",sep = "\t",header = T,check.names = F, row.names = 1)
testExpr=t(testExpr)
testExpr=as.matrix(testExpr)
testExpr[1:4,1:4]  
##colnames(testExpr)=paste0('test',colnames(testExpr))
dim(testExpr)  


calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

######
library(openxlsx)

ACC_IC50 <- read.xlsx("ACC_IC50.xlsx",sheet = "Sheet1")

####
exprSet=ACC_IC50[,-1]
y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
colnames <- colnames(exprSet)
cor_data_df <- data.frame(colnames)
for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
write.table(cor_data_df,"cor_data_df.txt",quote = F,sep = "\t",col.names = T,row.names = F)


####
ACC_IC50_cor <- read.xlsx("ACC_IC50_cor.xlsx",sheet = "Sheet1")
pathwaydurg <- read.xlsx("pathwaydurg.xlsx",sheet = "Sheet1")

pathway=merge(ACC_IC50_cor,pathwaydurg,by="Drug")
write.table(pathway,"pathway.txt",quote = F,sep = "\t",col.names = T,row.names = F)




