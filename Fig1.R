library(tglkmeans)
library(edgeR)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggfortify)
library(chromVAR)
library(motifmatchr)
library(chromVARmotifs)
library(JASPAR2018)
library(BSgenome.Mmusculus.UCSC.mm9)
library(ggpubr)
library(ggrepel)
library(viridis)

getMotifs <- function (species = "Homo sapiens", collection = "CORE", ...)
{
      opts <- list()
  opts["species"] <- species
    opts["collection"] <- collection
    opts <- c(opts, list(...))
      out <- TFBSTools::getMatrixSet(JASPAR2018::JASPAR2018, opts)
      if (!isTRUE(all.equal(TFBSTools::name(out), names(out))))
              names(out) <- paste(names(out), TFBSTools::name(out),
                                                          sep = "_")
        return(out)
}


counts = read_tsv("../git_site/files/counts/All_peaks_fixed_len.counts")

#Get CPM
sizeFactorsTable = read.table("sizeFactors/All.sizes", col.names=c("dataset", "size"))
#remove empty peaks
max_c = apply(counts[,7:29],1,max)
coints = counts[max_c>0,]

cpm = edgeR::cpm(as.matrix(counts[,7:29]),
                     lib.size=sizeFactorsTable$size,log=T,prior.count=5)
meta = counts[1:6,]

counts_var = read_tsv("../git_site/files/counts/Variable_peaks_kmeans.counts")

to_clust = counts_var[,12:34]
km = TGL_kmeans(to_clust,12,"euclid",reorder_func = "hclust",id_column=FALSE)
image(t(as.matrix(to_clust[order(km$cluster),])),,xaxt="n",yaxt="n",col=colorRampPalette(brewer.pal(n=9,name="Blues"))(8))

#ChromVAR
genome = BSgenome.Mmusculus.UCSC.mm9
motifs = getMotifs()

rowRanges = makeGRangesFromDataFrame(counts[,1:3])
colData = DataFrame(Treatment = names(counts[,7:29]),row.names=colnames(counts[,7:29]))
rse <- SummarizedExperiment(assays=SimpleList(counts=counts[,7:29]), rowRanges=rowRanges, colData=colData)

counts_filtered <- addGCBias(rse, genome = genome)
motif_ix <- matchMotifs(motifs, counts_filtered, genome = genome)
dev_TF <- computeDeviations(object = counts_filtered, annotations = motif_ix)
bg <- getBackgroundPeaks(object = counts_filtered)
expected <- computeExpectations(counts_filtered)
variability_TF <- computeVariability(dev_TF)

Highest_var_TF <- cbind(variability_TF, rownames(variability_TF))
Highest_var_TF <- Highest_var_TF %>% group_by(name) %>% dplyr::slice(which.max(variability))
Highest_var_TF <- Highest_var_TF[order(Highest_var_TF$variability, decreasing = TRUE),]
to_out = Highest_var_TF
Highest_var_TF <- Highest_var_TF[c(1:100),]
Dev_Highest_var_TF <- assays(dev_TF)[[1]]
Dev_Highest_var_TF <- as.data.frame(Dev_Highest_var_TF)
Dev_Highest_var_TF <- cbind(rownames(variability_TF), variability_TF$name, Dev_Highest_var_TF)
Dev_Highest_var_TF <- (subset(Dev_Highest_var_TF, Dev_Highest_var_TF$`rownames(variability_TF)` %in% Highest_var_TF$`rownames(variability_TF)`))
rownames(Dev_Highest_var_TF) <- Dev_Highest_var_TF$`variability_TF$name`
Dev_Highest_var_TF <- Dev_Highest_var_TF[,c(3:(dim(data)[2]+2))]
pheatmap(Dev_Highest_var_TF, cluster_rows = T, cluster_cols = FALSE, scale = "row")



