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
To be done
