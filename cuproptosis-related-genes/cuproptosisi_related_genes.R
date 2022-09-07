library(data.table)
RNASeq=fread("RNASeq.txt", data.table=FALSE)
RNASeq=RNASeq[!duplicated(RNASeq$sample),]
library(openxlsx)
RNASeq2_gsva_cu<- read.xlsx("RNASeq2_gsva_cu.xlsx",sheet = "Sheet1")
cu=RNASeq2_gsva_cu[,c(1,4)]
cu=read.xlsx("cu.xlsx",sheet = "Sheet1")

RNASeq=as.data.frame(t(RNASeq))
RNASeq=tibble::rownames_to_column(RNASeq)
colnames(RNASeq)=RNASeq[1,]
RNASeq=RNASeq[-1,]

RNASeq_cu=merge(cu,RNASeq,by="sample")


######
exprSet=RNASeq_cu[,-1]

colnames <- colnames(exprSet)

y <- as.numeric(exprSet[,"Cuproptosis_activity_score"])
cor_data_df <- data.frame(colnames)

for (i in 1:length(colnames)){
  test <- cor.test(as.numeric(exprSet[,i]),y,type="spearman")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
write.table(cor_data_df, "cor_data_df.txt",quote = F,sep = "\t",col.names = T,row.names = F)
#####
#####kegg go

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(GO.db)
library(clusterProfiler)
library(pathview)

# 导入数据
data <- read.table("gene.txt",header=TRUE)
data$GeneName <- as.character(data$GeneName)

# 转换基因名
transID = bitr(data$GeneName,
               fromType="SYMBOL",
               toType=c("ENSEMBL", "ENTREZID"),
               OrgDb="org.Hs.eg.db"
)

dir.create("GO")
dir.create("KEGG")

# GO_CC 注释
CC <- enrichGO(transID$ENTREZID,
               "org.Hs.eg.db",
               keyType="ENTREZID",
               ont="CC",
               pvalueCutoff=0.05,
               pAdjustMethod="BH",
               qvalueCutoff=0.1
)
CC <- setReadable(CC, OrgDb=org.Hs.eg.db)

pdf(file="./GO/GO_CC.pdf", bg="transparent")
barplot(CC, showCategory=10, title="GO_CC", font.size=12)

dev.off()

write.table(as.data.frame(CC@result), file="./GO/GO_CC.xls", sep="\t", row.names=F)

# GO_MF注释
MF <- enrichGO(transID$ENTREZID, "org.Hs.eg.db", keyType="ENTREZID", ont="MF", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
MF <- setReadable(MF, OrgDb=org.Hs.eg.db)

pdf(file="./GO/GO_MF.pdf", bg="transparent")
barplot(MF, showCategory=10, title="GO_MF", font.size=12)

dev.off()

write.table(as.data.frame(MF@result), file="./GO/GO_MF.xls", sep="\t", row.names=F)

# GO_BP注释
BP <- enrichGO(transID$ENTREZID, "org.Hs.eg.db", keyType="ENTREZID", ont="BP", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
BP <- setReadable(BP, OrgDb=org.Hs.eg.db)

pdf(file="./GO/GO_BP.pdf", bg="transparent")
barplot(BP, showCategory=10, title="GO_BP", font.size=12)

dev.off()

write.table(as.data.frame(BP@result), file="./GO/GO_BP.xls", sep="\t", row.names=F)

# KEGG 注释
kegg <- enrichKEGG(transID$ENTREZID, organism="hsa", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.1)
kegg <- setReadable(kegg, OrgDb=org.Hs.eg.db, keytype="ENTREZID")

pdf(file="./KEGG/KEGG.pdf", bg="transparent",width = 6,height=11)
dotplot(kegg, showCategory=30, font.size=12, title="KEGG") # + theme(axis.text.y = element_text(angle = 45))
barplot(kegg, showCategory=30, title="KEGG", font.size=12)
dev.off()

write.table(as.data.frame(kegg@result), file="./KEGG/kegg.xls", sep="\t", row.names=F)


