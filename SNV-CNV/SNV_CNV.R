library(openxlsx)
FDX1 <- read.xlsx("SNV.xlsx",sheet = "Sheet2")
LIAS <- read.xlsx("SNV.xlsx",sheet = "Sheet3")
LIPT1 <- read.xlsx("SNV.xlsx",sheet = "Sheet4")
DLD <- read.xlsx("SNV.xlsx",sheet = "Sheet5")
DLAT <- read.xlsx("SNV.xlsx",sheet = "Sheet6")
PDHA1 <- read.xlsx("SNV.xlsx",sheet = "Sheet7")
PDHB <- read.xlsx("SNV.xlsx",sheet = "Sheet8")
MTF1 <- read.xlsx("SNV.xlsx",sheet = "Sheet9")
GLS <- read.xlsx("SNV.xlsx",sheet = "Sheet10")
CDKN2A <- read.xlsx("SNV.xlsx",sheet = "Sheet11")

tumortype<- read.xlsx("tumortype.xlsx",sheet = "Sheet1")


library(dplyr)
SNV=dplyr::left_join(tumortype,FDX1,by="sample")
SNV=dplyr::left_join(SNV,LIAS,by="sample")
SNV=dplyr::left_join(SNV,LIPT1,by="sample")
SNV=dplyr::left_join(SNV,DLD,by="sample")
SNV=dplyr::left_join(SNV,DLAT,by="sample")
SNV=dplyr::left_join(SNV,PDHA1,by="sample")
SNV=dplyr::left_join(SNV,PDHB,by="sample")
SNV=dplyr::left_join(SNV,MTF1,by="sample")
SNV=dplyr::left_join(SNV,GLS,by="sample")
SNV=dplyr::left_join(SNV,CDKN2A,by="sample")
write.table(SNV,"SNV2.txt", sep="\t", quote=FALSE, col.names = TRUE,row.names = F)


library(ComplexHeatmap)
library(circlize)
#SNV
mut <- SNV[,-2]
cli <- SNV

rownames(mut) <- mut$sample
mat <- mut[,-1]
mat[is.na(mat)]<-""
mat[1:6,1:6]
mat=t(mat)
mat[1:6,1:6]

######
col <- c( "Missense_Mutation" = "#27AE60","Nonsense_Mutation" = "blue","Nonstop_Mutation" = "#B452CD",
          "Splice_Site" = "red", "In_Frame_Del" = "#EEC900","Frame_Shift_Del" = "#FF6A6A","Frame_Shift_Ins" = "#9400D3",
          "Mutiple" = "black")

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7, 
              gp = gpar(fill = "#e8e7e3", col = NA))
  },
  
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Missense_Mutation"], col = NA))
  },
  
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Nonsense_Mutation"], col = NA))
  },
  Nonstop_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Nonstop_Mutation"], col = NA))
  },
  
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Splice_Site"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["In_Frame_Del"], col = NA))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Frame_Shift_Del"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Frame_Shift_Ins"], col = NA))
  },
  Mutiple= function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.3,
              gp = gpar(fill = col["Mutiple"], col = NA))
  }
)
heatmap_legend <- list(title = "SNV",
                       at = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                              "Splice_Site", "In_Frame_Del","Frame_Shift_Del","Frame_Shift_Ins",
                              "Mutiple"),
                       labels = c("Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation",
                                  "Splice_Site", "In_Frame_Del","Frame_Shift_Del","Frame_Shift_Ins",
                                  "Mutiple"))
oncoPrint(mat ,  #数据
          alter_fun = alter_fun, col = col,
          #column_title = column_title,
          remove_empty_columns = T, 
          remove_empty_rows = T, 
          row_names_side = "left", #我们可以把基因信息放左边
          pct_side = "right",
          alter_fun_is_vectorized = FALSE,
          heatmap_legend_param = heatmap_legend,  #设置图例
)
pdata <- cli
head(pdata)
ha<-HeatmapAnnotation(Cancer_type=pdata$Cancer.type,
                      show_annotation_name = TRUE,
                      annotation_name_gp = gpar(fontsize = 7))


oncoPrint(mat ,
          alter_fun = alter_fun, col = col,
          remove_empty_columns = T, 
          remove_empty_rows = T, 
          row_names_side = "left", 
          pct_side = "right",
          alter_fun_is_vectorized = FALSE,
          heatmap_legend_param = heatmap_legend,
          bottom_annotation = ha 
)
#自定义样本顺序
s <- pdata[order(pdata$Cancer.type,pdata$CDKN2A,pdata$DLD,pdata$MTF1),]
sample_order <- as.character(s$sample)

#自定义颜色
#连续性变量设置颜色（外）
ha<-HeatmapAnnotation(Cancer_type=pdata$Cancer.type,
                      col = list(Cancer_type = c("ACC" =  "#FFB6C1", "BLCA" =  "#4E72B8", "BRCA" =  "#FAEBD7"
                                                 ,"CESC" =  "#00FFFF", "CHOL" =  "#7FFFAA", "COAD" =  "#FFE4C4"
                                                 ,"DLBC" =  "#0000FF", "ESCA" =  "#8A2BE2", "GBM" =  "#A52A2A"
                                                 ,"HNSC" =  "#7FFF00", "KICH" =  "#5F9EA0", "KIRC" =  "#DEB887"
                                                 ,"KIRP" =  "#FF7F50", "LAML" =  "#6495ED", "LGG" =  "#FFF8DC"
                                                 ,"LIHC" =  "#00008B", "LUAD" =  "#CD5C5C", "LUSC" =  "#A9A9A9"
                                                 ,"MESO" =  "#8B008B", "OV" =  "#FF8C00", "PAAD" =  "#E9967A"
                                                 ,"PCPG" =  "#FF1493", "PRAD" =  "#00BFFF", "READ" =  "#FFFAF0"
                                                 ,"SARC" =  "#FFD700", "SKCM" =  "#FF6347", "STAD" =  "#FF69B4"
                                                 ,"TGCT" =  "#66CDAA", "THCA" =  "#BA55D3", "THYM" =  "#FF4500"
                                                 ,"UCEC" =  "#DDA0DD", "UCS" =  "#008B8B", "UVM" =  "#FFFF00"
                      )),
                      show_annotation_name = TRUE,
                      annotation_name_gp = gpar(fontsize = 7))


oncoPrint(mat ,
          alter_fun = alter_fun, col = col,
          remove_empty_columns = T, 
          remove_empty_rows = T, 
          column_order = sample_order,
          row_names_side = "left", 
          pct_side = "right",
          alter_fun_is_vectorized = FALSE,
          heatmap_legend_param = heatmap_legend,
          bottom_annotation = ha 
)    
##8.5*4.5

#CNV
mut <- read.xlsx("CNV.xlsx",sheet = "Sheet1")
cli <- read.xlsx("CNV.xlsx",sheet = "Sheet2")

rownames(mut) <- mut$sample
mat <- mut[,-1]
mat[is.na(mat)]<-""
mat[1:6,1:6]
mat=t(mat)
mat[1:6,1:6]

######画图
col <- c( "Deep deletion" = "blue",
          "Amplification" = "red")

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7, 
              gp = gpar(fill = "#e8e7e3", col = NA))
  },
  
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Amplification"], col = NA))
  },
  
  "Deep deletion" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.7,
              gp = gpar(fill = col["Deep deletion"], col = NA))
  }
  
 )
heatmap_legend <- list(title = "CNV",
                       at = c("Amplification", "Deep deletion"
                              ),
                       labels = c("Amplification", "Deep deletion"
                                  ))
oncoPrint(mat ,  #数据
          alter_fun = alter_fun, col = col,
          #column_title = column_title,
          remove_empty_columns = T, 
          remove_empty_rows = T, 
          row_names_side = "left", #我们可以把基因信息放左边
          pct_side = "right",
          alter_fun_is_vectorized = FALSE,
          heatmap_legend_param = heatmap_legend,  #设置图例
)
pdata <- cli
head(pdata)
ha<-HeatmapAnnotation(Cancer_type=pdata$Cancer.type,
                      show_annotation_name = TRUE,
                      annotation_name_gp = gpar(fontsize = 7))


oncoPrint(mat ,
          alter_fun = alter_fun, col = col,
          remove_empty_columns = T, 
          remove_empty_rows = T, 
          row_names_side = "left", 
          pct_side = "right",
          alter_fun_is_vectorized = FALSE,
          heatmap_legend_param = heatmap_legend,
          bottom_annotation = ha 
)
#自定义样本顺序
s <- pdata[order(pdata$Cancer.type,pdata$CDKN2A,pdata$DLD,pdata$MTF1),]
sample_order <- as.character(s$sample)

#自定义颜色
#连续性变量设置颜色（外）
ha<-HeatmapAnnotation(Cancer_type=pdata$Cancer.type,
                      col = list(Cancer_type = c("ACC" =  "#FFB6C1", "BLCA" =  "#4E72B8", "BRCA" =  "#FAEBD7"
                                                 ,"CESC" =  "#00FFFF", "CHOL" =  "#7FFFAA", "COAD" =  "#FFE4C4"
                                                 ,"DLBC" =  "#0000FF", "ESCA" =  "#8A2BE2", "GBM" =  "#A52A2A"
                                                 ,"HNSC" =  "#7FFF00", "KICH" =  "#5F9EA0", "KIRC" =  "#DEB887"
                                                 ,"KIRP" =  "#FF7F50", "LAML" =  "#6495ED", "LGG" =  "#FFF8DC"
                                                 ,"LIHC" =  "#00008B", "LUAD" =  "#CD5C5C", "LUSC" =  "#A9A9A9"
                                                 ,"MESO" =  "#8B008B", "OV" =  "#FF8C00", "PAAD" =  "#E9967A"
                                                 ,"PCPG" =  "#FF1493", "PRAD" =  "#00BFFF", "READ" =  "#FFFAF0"
                                                 ,"SARC" =  "#FFD700", "SKCM" =  "#FF6347", "STAD" =  "#FF69B4"
                                                 ,"TGCT" =  "#66CDAA", "THCA" =  "#BA55D3", "THYM" =  "#FF4500"
                                                 ,"UCEC" =  "#DDA0DD", "UCS" =  "#008B8B", "UVM" =  "#FFFF00"
                      )),
                      show_annotation_name = TRUE,
                      annotation_name_gp = gpar(fontsize = 7))


oncoPrint(mat ,
          alter_fun = alter_fun, col = col,
          remove_empty_columns = T, 
          remove_empty_rows = T, 
          column_order = sample_order,
          row_names_side = "left", 
          pct_side = "right",
          alter_fun_is_vectorized = FALSE,
          heatmap_legend_param = heatmap_legend,
          bottom_annotation = ha 
)                     
##9*4.5




