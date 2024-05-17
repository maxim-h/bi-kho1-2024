#install.packages('Seurat')
#install.packages('Signac')

library(tidyverse)
library(Seurat)
library(Signac)
library(data.table)
library(dplyr)
library(tibble)
library("AnnotationDbi")
library("org.Hs.eg.db")

#import seurat object
seurat <- readRDS("seurat_v4.rds")

#plot 
seurat$branch <- seurat[[]] %>% mutate(branch = case_when(
  branch_1 ~ "Myeloid",
  branch_2 ~ "Lymphoid",
  branch_3 ~ "Erythroid",
)) %>% pull(branch)
DimPlot(seurat, group.by = "branch")
DimPlot(seurat, group.by = "scenic_clusters")

#this is how we need the extracted data to look
rna_seq <- fread("rna.pseudobulkFromClusters_cell_type_mean.tsv.gz", header = T)

aggr_RNA_br1 <- AggregateExpression(subset(seurat, subset = branch_1), assays = c("RNA"), group.by = "scenic_clusters")
aggr_RNA_br1_tsv <- as.data.frame(aggr_RNA_br1)


#add 'SYMBOL' as a name for the first column
aggr_RNA_br1_tsv <- aggr_RNA_br1_tsv %>% 
  rownames_to_column(var = "SYMBOL")

#convert HGNC symbols to ENSEMBL
aggr_RNA_br1_tsv$SYMBOL = mapIds(org.Hs.eg.db,
                  keys=aggr_RNA_br1_tsv$SYMBOL, 
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")

#rename the first column to "ENSEMBL"
names(aggr_RNA_br1_tsv)[names(aggr_RNA_br1_tsv) == "SYMBOL"] <- "ENSEMBL"

#remove prefix 'RNA' from the column names
colnames(aggr_RNA_br1_tsv) <- gsub("RNA\\.", "", colnames(aggr_RNA_br1_tsv))

#how many NA IDs
sum(is.na(aggr_RNA_br1_tsv$ENSEMBL))

#save the resulting table
fwrite(x = aggr_RNA_br1_tsv, file = "RNA_br1.tsv")

#############RNA_br2
aggr_RNA_br2 <- AggregateExpression(subset(seurat, subset = branch_2), assays = c("RNA"), group.by = "scenic_clusters")
aggr_RNA_br2_tsv <- as.data.frame(aggr_RNA_br2)

#add 'SYMBOL' as a name for the first column
aggr_RNA_br2_tsv <- aggr_RNA_br2_tsv %>% 
  rownames_to_column(var = "SYMBOL")

#convert HGNC symbols to ENSEMBL
aggr_RNA_br2_tsv$SYMBOL = mapIds(org.Hs.eg.db,
                                 keys=aggr_RNA_br2_tsv$SYMBOL, 
                                 column="ENSEMBL",
                                 keytype="SYMBOL",
                                 multiVals="first")

#rename the first column to "ENSEMBL"
names(aggr_RNA_br2_tsv)[names(aggr_RNA_br2_tsv) == "SYMBOL"] <- "ENSEMBL"

#remove prefix 'RNA' and full stops from the column names
colnames(aggr_RNA_br2_tsv) <- gsub("RNA\\.", "", colnames(aggr_RNA_br2_tsv))

#how many NA IDs
sum(is.na(aggr_RNA_br2_tsv$ENSEMBL))

#save the resulting table
fwrite(x = aggr_RNA_br2_tsv, file = "RNA_br2.tsv")

############RNA_br3

aggr_RNA_br3 <- AggregateExpression(subset(seurat, subset = branch_3), assays = c("RNA"), group.by = "scenic_clusters")
aggr_RNA_br3_tsv <- as.data.frame(aggr_RNA_br3)

#add 'SYMBOL' as a name for the first column
aggr_RNA_br3_tsv <- aggr_RNA_br3_tsv %>% 
  rownames_to_column(var = "SYMBOL")

#convert HGNC symbols to ENSEMBL
aggr_RNA_br3_tsv$SYMBOL = mapIds(org.Hs.eg.db,
                                 keys=aggr_RNA_br3_tsv$SYMBOL, 
                                 column="ENSEMBL",
                                 keytype="SYMBOL",
                                 multiVals="first")

#rename the first column to "ENSEMBL"
names(aggr_RNA_br3_tsv)[names(aggr_RNA_br3_tsv) == "SYMBOL"] <- "ENSEMBL"

#remove prefix 'RNA' and full stops from the column names
colnames(aggr_RNA_br3_tsv) <- gsub("RNA\\.", "", colnames(aggr_RNA_br3_tsv))

#how many NA IDs
sum(is.na(aggr_RNA_br3_tsv$ENSEMBL))

#save the resulting table
fwrite(x = aggr_RNA_br3_tsv, file = "RNA_br3.tsv")

############ATAC_br1

aggr_atac_br1 <- AggregateExpression(subset(seurat, subset = branch_1), assays = c("ATAC"), group.by = "scenic_clusters")
aggr_atac_br1_tsv <- as.data.frame(aggr_atac_br1)

#add 'SYMBOL' as a name for the first column
aggr_atac_br1_tsv <- aggr_atac_br1_tsv %>% 
  rownames_to_column(var = "SYMBOL")

#rename the first column to "ENSEMBL"
names(aggr_atac_br1_tsv)[names(aggr_atac_br1_tsv) == "SYMBOL"] <- "peakID"

#remove prefix 'RNA' and full stops from the column names
colnames(aggr_atac_br1_tsv) <- gsub("ATAC\\.", "", colnames(aggr_atac_br1_tsv))

#how many NA IDs
sum(is.na(aggr_atac_br1_tsv$peakID))

#save the resulting table
fwrite(x = aggr_atac_br1_tsv, file = "ATAC_br1.tsv")


############ATAC_br2

aggr_atac_br2 <- AggregateExpression(subset(seurat, subset = branch_2), assays = c("ATAC"), group.by = "scenic_clusters")
aggr_atac_br2_tsv <- as.data.frame(aggr_atac_br2)

#add 'SYMBOL' as a name for the first column
aggr_atac_br2_tsv <- aggr_atac_br2_tsv %>% 
  rownames_to_column(var = "SYMBOL")

#rename the first column to "ENSEMBL"
names(aggr_atac_br2_tsv)[names(aggr_atac_br2_tsv) == "SYMBOL"] <- "peakID"

#remove prefix 'RNA' and full stops from the column names
colnames(aggr_atac_br2_tsv) <- gsub("ATAC\\.", "", colnames(aggr_atac_br2_tsv))

#how many NA IDs
sum(is.na(aggr_atac_br2_tsv$peakID))

#save the resulting table
fwrite(x = aggr_atac_br2_tsv, file = "ATAC_br2.tsv")


############ATAC_br3

aggr_atac_br3 <- AggregateExpression(subset(seurat, subset = branch_3), assays = c("ATAC"), group.by = "scenic_clusters")
aggr_atac_br3_tsv <- as.data.frame(aggr_atac_br3)

#add 'SYMBOL' as a name for the first column
aggr_atac_br3_tsv <- aggr_atac_br3_tsv %>% 
  rownames_to_column(var = "SYMBOL")

#rename the first column to "ENSEMBL"
names(aggr_atac_br3_tsv)[names(aggr_atac_br3_tsv) == "SYMBOL"] <- "peakID"

#remove prefix 'RNA' and full stops from the column names
colnames(aggr_atac_br3_tsv) <- gsub("ATAC\\.", "", colnames(aggr_atac_br3_tsv))

#how many NA IDs
sum(is.na(aggr_atac_br3_tsv$peakID))

#save the resulting table
fwrite(x = aggr_atac_br3_tsv, file = "ATAC_br3.tsv")



#this is how we need the extracted data to look
atac_seq <- fread("atac.pseudobulkFromClusters_cell_type_mean.tsv.gz", header = T)

#metadata table
metadata <- seurat[[]] %>% 
  group_by(scenic_clusters) %>%
  summarise(nCells = n())

#save metadata as tsv
write.table(metadata, file = "metadata.tsv", sep="\t")

#meta <- read.table("metadata.tsv", header = T)


###################################################ANALYSIS br1

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


atac_br1 <- fread("ATAC_br1.tsv", header = T)
rna_br1 <- fread("RNA_br1.tsv", header = T)
meta_data <- read.table("metadata.tsv",  header = T)

meta.l = list(name = "branch1", date = "14.03.24")

GRN_br1 <- initializeGRN(objectMetadata = meta.l, genomeAssembly = "hg38")

grn_br1 <- addData(GRN_br1, sampleMetadata = meta_data, counts_peaks = atac_br1, normalization_peaks = "none", 
               counts_rna = rna_br1, normalization_rna = "none")

####QC
grn_br1 <- plotPCA_all(grn_br1)
grn_br1 <- addTFBS(grn_br1, motifFolder = 'PWMScan_HOCOMOCOv12/H12INVIVO/', translationTable = 'translationTable.csv',
               translationTable_sep = ' ', fileEnding = '.bed.gz')

grn_br1 <- overlapPeaksAndTFBS(grn_br1, nCores = 1)
grn_br1 <- addConnections_TF_peak(grn_br1)

saveRDS(grn_br1, 'grn_br1.rds')

grn_br1 <- readRDS('grn_br1.rds')

grn_br1 <- addConnections_peak_gene(grn_br1, nCores = 1)

#change FDR threshold
grn_br1 <- filterGRNAndConnectGenes(grn_br1, TF_peak.fdr.threshold = 0.05, forceRerun = TRUE)

conections.all <- getGRNConnections(grn_br1)
grn_br1 <- generateStatsSummary(grn_br1)

grn_br1  <-  build_eGRN_graph(grn_br1)

grn_br1 <- visualizeGRN(grn_br1, maxEdgesToPlot = 5300)
saveRDS(grn_br1, 'grn_0.05_br1.rds')
#grn <- readRDS('grn_0.05_br1.rds')

#network and enrichment analyses for filtered connections

grn_br1 <- performAllNetworkAnalyses(grn_br1)


###################################################ANALYSIS br2

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


atac_br2 <- fread("ATAC_br2.tsv", header = T)
rna_br2 <- fread("RNA_br2.tsv", header = T)
meta_data <- read.table("metadata.tsv",  header = T)

meta.l = list(name = "branch2", date = "14.03.24")

GRN_br2 <- initializeGRN(objectMetadata = meta.l, genomeAssembly = "hg38")

grn_br2 <- addData(GRN_br2, sampleMetadata = meta_data, counts_peaks = atac_br2, normalization_peaks = "none", 
                   counts_rna = rna_br2, normalization_rna = "none")

####QC
grn_br2 <- plotPCA_all(grn_br2)
grn_br2 <- addTFBS(grn_br2, motifFolder = 'PWMScan_HOCOMOCOv12/H12INVIVO/', translationTable = 'translationTable.csv',
                   translationTable_sep = ' ', fileEnding = '.bed.gz')

grn_br2 <- overlapPeaksAndTFBS(grn_br2, nCores = 1)
grn_br2 <- addConnections_TF_peak(grn_br2)

saveRDS(grn_br2, 'grn_br2.rds')

grn_br2 <- readRDS('grn_br2.rds')

grn_br2 <- addConnections_peak_gene(grn_br2, nCores = 1)

#change FDR threshold
#No connections passed the filter steps. Rerun the function and be less stringent (0.05)
grn_br2 <- filterGRNAndConnectGenes(grn_br2, TF_peak.fdr.threshold = 0.08, forceRerun = TRUE)

conections.all <- getGRNConnections(grn_br2)
grn_br2 <- generateStatsSummary(grn_br2)

grn_br2  <-  build_eGRN_graph(grn_br2)

grn_br2 <- visualizeGRN(grn_br2, maxEdgesToPlot = 5300)
saveRDS(grn_br2, 'grn_br2.rds')
#grn <- readRDS('grn_0.05_semifinal.rds')

#network and enrichment analyses for filtered connections

grn_br2 <- performAllNetworkAnalyses(grn_br2)


###################################################ANALYSIS br3

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


atac_br3 <- fread("ATAC_br3.tsv", header = T)
rna_br3 <- fread("RNA_br3.tsv", header = T)
meta_data <- read.table("metadata.tsv",  header = T)

meta.l = list(name = "branch3", date = "14.03.24")

GRN_br3 <- initializeGRN(objectMetadata = meta.l, genomeAssembly = "hg38")

grn_br3 <- addData(GRN_br3, sampleMetadata = meta_data, counts_peaks = atac_br3, normalization_peaks = "none", 
                   counts_rna = rna_br3, normalization_rna = "none")

####QC
grn_br3 <- plotPCA_all(grn_br3)
grn_br3 <- addTFBS(grn_br3, motifFolder = 'PWMScan_HOCOMOCOv12/H12INVIVO/', translationTable = 'translationTable.csv',
                   translationTable_sep = ' ', fileEnding = '.bed.gz')

grn_br3 <- overlapPeaksAndTFBS(grn_br3, nCores = 1)

###stopped here
grn_br3 <- addConnections_TF_peak(grn_br3)

saveRDS(grn_br3, 'grn_br3.rds')

grn_br3 <- readRDS('grn_REAL_intermediary.rds')

grn_br3 <- addConnections_peak_gene(grn_br3, nCores = 1)

#change FDR threshold
grn_br3 <- filterGRNAndConnectGenes(grn_br3, TF_peak.fdr.threshold = 0.05, forceRerun = TRUE)

conections.all <- getGRNConnections(grn_br3)
grn_br3 <- generateStatsSummary(grn_br3)

grn_br3  <-  build_eGRN_graph(grn_br3)

grn_br3 <- visualizeGRN(grn_br3, maxEdgesToPlot = 5300)
saveRDS(grn_br3, 'grn_0.05_semifinal.rds')
#grn <- readRDS('grn_0.05_semifinal.rds')

#network and enrichment analyses for filtered connections

grn_br3 <- performAllNetworkAnalyses(grn_br3)

