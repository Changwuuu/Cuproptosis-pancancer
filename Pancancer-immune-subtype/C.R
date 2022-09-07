
######combine the pancancer immune subtypes and cuscore. Then using Sangerbox to plot the box plot. Please cite:Shen, et al. 2022. Sangerbox: A comprehensive, interaction-friendly clinical bioinformatics analysis platform. iMeta 1(3): e36. https://doi.org/10.1002/imt2.36
library(openxlsx)
C_PAN <- read.xlsx("C_PAN.xlsx",sheet = "Sheet1")
RNA_Cuscore <- read.xlsx("RNA_Cuscore.xlsx",sheet = "Sheet1")

C_PAN_CU=merge(C_PAN,RNA_Cuscore,by="sample")
write.table(C_PAN_CU,file = "C_PAN_CU.txt", quote=F,sep = "\t", col.names = T, row.names = F)
