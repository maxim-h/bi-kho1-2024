
suppressPackageStartupMessages({
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
  library(argparser)
  library(here)
  library(lubridate)
  library(GRaNPA)
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(tibble)
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  library("qs")
})

here()

p <- add_argument(p, "--grn_rds", help = "Path to GRaNIE rds", default = here("granie_output"))
p <- add_argument(p, "--de_data", help = "Path to diff expr data for validation", default = here("folder_name"))
p <- add_argument(p, "--de_pval_th", help = "P-value threshold for DE", default = 0.2)
p <- add_argument(p, "--logFC_th", help = "LogFC threshold for DE", default = 0)
p <- add_argument(p, "--n_cores", help = "Number of cores", default = 10)
p <- add_argument(p, "--importance_tf", help = "Algorithm for finding most important TFs", default = permutation)
p <- add_argument(p, "--ml_type", help = "Regression or classification", default = regression)
p <- add_argument(p, "--output_folder_granpa", help = "Output folder name")

# Parse the command line arguments
args <- parse_args(p)

grn_rds <- readRDS(args$grn_rds)
de_data <- qread(args$de_data)
de_pval_th <- args$de_pval_th
logFC_th <- args$logFC_th
n_cores <- args$n_cores
importance_tf <- args$importance_tf
ml_type <- args$ml_type
output_folder_granpa <- args$output_folder_granpa

granpa_result = GRaNPA::GRaNPA_main_function(DE_data = de_data, 
                                                GRN_matrix_filtered = grn_rds,
                                                DE_pvalue_th = de_pval_th,
                                                logFC_th = logFC_th,
                                                num_run = 3,
                                                num_run_CR = 2,
                                                num_run_random = 3,
                                                cores = n_cores,
                                                importance = importance_tf,
                                                ML_type = ml_type,
                                                control = "cv",
                                                train_part = 1)

GRaNPA::plot_GRaNPA_density(GRaNPA.object = granpa_result, plot_name = "density.pdf", outputFolder = ".", width = 4, height = 4)
GRaNPA::plot_GRaNPA_scatter(GRaNPA.object = granpa_result, plot_name = "scatter.pdf", outputFolder = ".", width = 4, height = 4) 
GRaNPA::plot_GRaNPA_TF_imp(GRaNPA.object = ggranpa_result, plot_name = "TF_imp.pdf", outputFolder = ".", width = 4, height = 4) 