BiocManager::install(c("org.Hs.eg.db",
                       "BSgenome.Hsapiens.UCSC.hg38",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene"))
install.packages("WGCNA")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GRaNIE")
install.packages("devtools")

devtools::install_gitlab("grp-zaugg/GRaNIE", host =
                           "git.embl.de", subdir = "src/GRaNIE", force = TRUE)

BiocManager::install(c("ChIPseeker"))

library(WGCNA)
library(tidyverse)
library(GRaNIE)
library(biomaRt)
library(data.table)
library(R.utils)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ChIPseeker)


atac <- fread("atac.pseudobulkFromClusters_cell_type_mean.tsv.gz", header = T)
rna_seq <- fread("rna.pseudobulkFromClusters_cell_type_mean.tsv.gz", header = T)
meta_data <- fread("metadata_cell_type.tsv.gz")

meta.l = list(name = "tutorial", date = "13.02.24")

GRN <- initializeGRN(objectMetadata = meta.l, genomeAssembly = "hg38")

grn <- addData(GRN, sampleMetadata = meta_data, counts_peaks = atac, normalization_peaks = "none", 
        counts_rna = rna_seq, normalization_rna = "none")

####QC
grn <- plotPCA_all(grn)
grn <- addTFBS(grn, motifFolder = 'PWMScan_HOCOMOCOv12/H12INVIVO/', translationTable = 'translationTable.csv',
               translationTable_sep = ' ', fileEnding = '.bed.gz')

grn <- overlapPeaksAndTFBS(grn, nCores = 1)
grn <- addConnections_TF_peak(grn)

saveRDS(grn, 'grn_REAL_intermediary.rds')

grn <- readRDS('grn_REAL_intermediary.rds')

grn <- addConnections_peak_gene(grn, nCores = 1)

#change FDR threshold
grn <- filterGRNAndConnectGenes(grn, TF_peak.fdr.threshold = 0.05, forceRerun = TRUE)

conections.all <- getGRNConnections(grn)
grn <- generateStatsSummary(grn)

grn  <-  build_eGRN_graph(grn)

grn <- visualizeGRN(grn, maxEdgesToPlot = 5300)
saveRDS(grn, 'grn_0.05_semifinal.rds')
#grn <- readRDS('grn_0.05_semifinal.rds')

#network and enrichment analyses for filtered connections

grn <- performAllNetworkAnalyses(grn)


