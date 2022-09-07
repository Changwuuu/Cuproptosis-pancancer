library(openxlsx)

GDSC2_cell_cuscore <- read.xlsx("GDSC2_cell_cuscore.xlsx",sheet = "Sheet1")
GDSC2_RES <- read.xlsx("GDSC2_RES.xlsx",sheet = "Sheet1")

GDSC2_cell_cuscore=as.data.frame(t(GDSC2_cell_cuscore))
GDSC2_cell_cuscore=tibble::rownames_to_column(GDSC2_cell_cuscore)
colnames(GDSC2_cell_cuscore)=GDSC2_cell_cuscore[1,]
GDSC2_cell_cuscore=GDSC2_cell_cuscore[-1,]
GDSC2_cuscore_ic50=merge(GDSC2_cell_cuscore,GDSC2_RES,by="sample")

exprSet=GDSC2_cuscore_ic50
exprSet=exprSet[,c(4:202)]

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
pathwaydurg <- read.xlsx("pathwaydurg.xlsx",sheet = "Sheet1")
library(readxl)
cu_ic50_39 <- read_excel("cu_ic50_cor.xlsx", 
                          sheet = "Sheet2")

pathway39=merge(cu_ic50_39,pathwaydurg,by="Drug")
write.table(pathway39,"pathway39.txt",quote = F,sep = "\t",col.names = T,row.names = F)
