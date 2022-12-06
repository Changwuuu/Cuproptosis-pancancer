# Cuproptosis-pancancer
These codes are for cuproptosis pancancer analyses. Please find "README.md" in every folder. 

Partial analysis and visualization are performed by Sangerbox (http://vip.sangerbox.com/home.html).

The html file of this study can be found in: https://pan.baidu.com/s/1yH7bxmmQTYTb-s5EnCWrKg 
keyword：p3it

# Please cite:

1.Changwu Wu, Jun Tan, Xiangyu Wang, Chaoying Qin, Wenyong Long, Yimin Pan, Yuzhe Li, Qing Liu. 2022. Pan‐cancer analyses reveal molecular and clinical characteristics of cuproptosis regulators. iMeta e68. https://doi.org/10.1002/imt2.68

2.Shen, W. et al. Sangerbox: A comprehensive, interaction-friendly clinical bioinformatics analysis platform. iMeta 1, e36 (2022).

3.Gu, Z. Complex heatmap visualization. iMeta 1, e43 (2022).

# File description

1.Project.Rmd: This is the R markdown file of all relevant codes in this study.  Ordered as each figure appears in the manuscript.

2.SNV_CNV: This folder contains files used on SNV and CNV analysis. It is related to Figure 2.

3.mRNA: This folder contains files used on mRNA difference analysis between paired tumor and normal tissue, relationship analysis between CNV and mRNA expression and the relationship between cuproptosis activity and hallmark score. This folder is related to Figure 3A and 7.

4.GSVA: This folder contains files used on GSVA algorithm. Including the generation of cuproptosis scores and hallmark score in all datasets. This folder is related to Figure 3,6-9.

5.Cluster: This folder contains files used on cuproptosis-related cluster analysis. Including the difference of cuproptosis scores between clusters. This folder is related to Figure 3B-F.

6.cuproptosis_related_genes: This folder contains files used on cuproptosis-related genes analysis. Including the identification of cuproptosis-related genes and enrichment analysis. This folder is related to Figure S6.

7.Methylation: This folder contains files used on methylation analysis of cuproptosis regulators.  Including the methylation difference analysis between paired tumor and normal tissue, the relationship between methylation value and mRNA expression, survival analysis of methylation.This folder is related to Figure 4.

8.miRNA: This folder contains files used on the identification of potential miRNAs regulated cuproptosis regulators.  It is related to Figure 5A.

9.lncRNA: This folder contains files used on the identification of potential lncRNAs regulated cuproptosis regulators.  It is related to Figure 5B.

10.TF: This folder contains files used on the identification of potential transcription factors regulated cuproptosis regulators.  It is related to Figure 5C.

11.survival:  This folder contains files used on the survival analysis of cuproptosis in TCGA cohort and the validation cohort (CGGA).It is related to Figure 6.

12.Cell_senescence_and_cell_cycle: This folder contains files used on the correlation analysis between cell senescence/cell cycle markers and cuproptosis activity. It is related to Figure S9.

13.immune_analysis:This folder contains files used on the immune-related analysis of cuproptosis. Including the correlation between cuproptosis activity and ImmuneScore, the correlation between cuproptosis activity and immune regulators, the correlation between cuproptosis activity and immune cell infiltration and the diffrenence of cuproptosis activity among six immune subtypes. This folder is related to Figure 8.

14.drug_sensitivity: This folder contains files used on drug sensitivity analysis of cuproptosis. Including the correlation between cuproptosis activity and IC50 values of 189 drugs in all tumors and the drug sensitivity of individual tumor samples in ACC. This folder is related to Figure 9A.

15.immunotherapy: This folder contains files used on the relationship analysis of immunotherapy outcome and cuproptosis activity. This folder is related to Figure 9B-F.




