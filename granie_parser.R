#!/usr/bin/env -S Rscript --vanilla

suppressPackageStartupMessages({
library(WGCNA)
library(tidyverse)
library(GRaNIE)
library(biomaRt)
library(data.table)
library(R.utils)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(argparser)
library(here)
library(lubridate)
})

here() #takes parameters to join into a path

# Create a parser
p <- arg_parser("Run GRaNIE")

# Add command line arguments
p <- add_argument(p, "--meta_data", help = "Path to TSV with metadata")
p <- add_argument(p, "--atac", help = "Path to TSV with atac data")
p <- add_argument(p, "--rna_seq", help = "Path to RNA seq counts data")
p <- add_argument(p, "--output_folder", help = "Folder name for output")
p <- add_argument(p, "--branch", help = "branch name")
p <- add_argument(p, "--hocomoco_path", help = "Path to hocomoco folder", default = here("PWMScan_HOCOMOCOv12/H12INVIVO"))
p <- add_argument(p, "--correlation_method",  help = "Correlation method in addConnections_peak_gene", default = "pearson")
p <- add_argument(p, "--promoter_range", help = "Promoter range in addConnections_peak_gene, bp", default = 250000)
p <- add_argument(p, "--TF_peak_FDR", help = "FDR threshold in filterGRNAndConnectGenes", default = 0.2)
p <- add_argument(p, "--n_cores", help = "Number of cores in addConnections_peak_gene and overlapPeaksAndTFBS", default = 2)

# Parse the command line arguments
args <- parse_args(p)

meta_data <- fread(args$meta_data)
atac <- fread(args$atac, header = T)
rna_seq <- fread(args$rna_seq, header = T)

output_folder <- args$output_folder
hocomoco_path <- args$hocomoco_path
correlation_method <- args$correlation_method
promoter_range <- args$promoter_range
TF_peak_FDR <- args$TF_peak_FDR
branch <- args$branch
n_cores <- args$n_cores

meta.l = list(name = "GRaNIE", date = now()) 

GRN <- initializeGRN(objectMetadata = meta.l, outputFolder = output_folder, genomeAssembly = "hg38")

grn <- addData(GRN, sampleMetadata = meta_data, counts_peaks = atac, normalization_peaks = "none", 
        counts_rna = rna_seq, normalization_rna = "none")

####QC
grn <- plotPCA_all(grn)
grn <- addTFBS(grn, motifFolder = hocomoco_path, translationTable = 'translationTable.csv',
               translationTable_sep = ' ', fileEnding = '.bed.gz')

grn <- overlapPeaksAndTFBS(grn, nCores = n_cores)
grn <- addConnections_TF_peak(grn)

grn <- addConnections_peak_gene(grn, corMethod = correlation_method, promoterRange = promoter_range, nCores = n_cores) 


#change FDR threshold
grn <- filterGRNAndConnectGenes(grn, TF_peak.fdr.threshold = TF_peak_FDR, forceRerun = FALSE)

conections.all <- getGRNConnections(grn)
grn <- generateStatsSummary(grn)

grn  <-  build_eGRN_graph(grn)

grn <- visualizeGRN(grn, maxEdgesToPlot = 5300)

#network and enrichment analyses for filtered connections

grn <- performAllNetworkAnalyses(grn)

###add granpa and other validations as a separate script!
###current output folder will serve as input for that

saveRDS(grn, here(output_folder, paste(branch, 'grn.rds', sep='_')))
