# Maor_Nof_2020

[TOC]

#### ATAC-seq
ATAC-seq data at WashU genome browser : http://epigenomegateway.wustl.edu/browser/?genome=mm9&session=ZiyQfvwMSC&statusId=990979562

- Peaks:
    - files/peaks/All_peaks_fixed_len.bed - 321788 peaks, fixed length
    - files/peaks/All_peaks_variable_len.bed - 321788 peaks, variable length
- Counts:
    - files/counts/All_peaks_fixed_len.counts - Counts in fixed length peaks
    - files/counts/All_peaks_variable_len.counts - Counts in variable length peaks
    - files/counts/Variable_peaks_kmeans.counts - 3981 variable peaks (CPM)

To get CPM from count files
```R
library(edgeR)
#Get CPM
sizeFactorsTable = read.table("sizeFactors/All.sizes", col.names=c("dataset", "size"))

cpm = edgeR::cpm(as.matrix(counts[,7:29]),
    lib.size=sizeFactorsTable$size,log=T,prior.count=5)
```

#### RNA-seq

RNA-seq DEseq2 pairwise comparison file :

- DEseq2:
    - files/diffexp/GFP_none_day1.deseq.csv - GFP vs Ctrl Day 1
    - files/diffexp/GFP_none_day2.deseq.csv - GFP vs Ctrl Day 2
    - files/diffexp/GFP_none_day3.deseq.csv - GFP vs Ctrl Day 3
    - files/diffexp/GFP_pr_day1.deseq.csv - GFP vs PR Day 1
    - files/diffexp/GFP_pr_day2.deseq.csv - GFP vs PR Day 2
    - files/diffexp/GFP_pr_day3.deseq.csv - GFP vs PR Day 3
    - files/diffexp/GFP_tdp_day1.deseq.csv - GFP vs TDP Day 1
    - files/diffexp/GFP_tdp_day2.deseq.csv - GFP vs TDP Day 2
    - files/diffexp/GFP_tdp_day3.deseq.csv - GFP vs TDP Day 3
    - files/diffexp/PR_day1_day3.deseq.csv - PR Day 1 vs Day 3
    - files/diffexp/TDP_day1_day3.deseq.csv - TDP Day 1 vs Day 3
    - files/diffexp/ctrl_Day1_Day3.deseq.csv - Ctrl Day 1 vs Day 3
    - files/diffexp/ctrl_pr_day1.deseq.csv - Ctrl vs PR Day 1
    - files/diffexp/ctrl_pr_day2.deseq.csv - Ctrl vs PR Day 2
    - files/diffexp/ctrl_pr_day3.deseq.csv - Ctrl vs PR Day 3
    - files/diffexp/ctrl_tdp_day1.deseq.csv - Ctrl vs TDP Day 1
    - files/diffexp/ctrl_tdp_day2.deseq.csv - Ctrl vs TDP Day 2
    - files/diffexp/ctrl_tdp_day3.deseq.csv - Ctrl vs TDP Day 3
