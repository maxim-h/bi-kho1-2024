#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(GRaNPA)
  library(GRaNIE)
  library(argparser)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(here)
})

#create a parser
p <- arg_parser("Results summary")

#add command line arguments
p <- add_argument(p, "--granpa_rds_list", help = "Comma-separated list of granpa rds files")
p <- add_argument(p, "--TF_peak_FDR_list", help = "Comma-separated list of corresponding TF peak FDRs")
p <- add_argument(p, "--peak_gene_FDR_list", help = "Comma-separated list of corresponding peak gene FDRs")
p <- add_argument(p, "--correlation_list", help = "Comma-separated list of corresponding correlation method")
p <- add_argument(p, "--promoter_range_list", help = "Comma-separated list of corresponding promoter range")
p <- add_argument(p, "--tf_database_list", help = "Comma-separated list of TF database")
p <- add_argument(p, "--output_folder", help = "Folder name for output")

#parse the command line arguments
args <- parse_args(p)

granpa_rds_list <- strsplit(args$granpa_rds_list, ", ")[[1]]
TF_peak_FDR_list <- strsplit(args$TF_peak_FDR_list, ",")[[1]]
peak_gene_FDR_list <- strsplit(args$peak_gene_FDR_list, ",")[[1]]
correlation_list <- strsplit(args$correlation_list, ",")[[1]]
promoter_range_list <- strsplit(args$promoter_range_list, ",")[[1]]
tf_database_list <- strsplit(args$tf_database_list, ",")[[1]]

output_folder <- args$output_folder
dir.create(output_folder, showWarnings = FALSE)

#read each file individually
data_list <- lapply(granpa_rds_list, readRDS)

mean_rsquared_list <- lapply(data_list, function(granpa_rds) {
  rsquared_values <- granpa_rds[["normal_models"]][[1]][["results"]][["Rsquared"]]
  mean(rsquared_values)
})

mean_rsquared_list <- lapply(data_list, function(granpa_rds) {
  if (is.null(granpa_rds)) {
    NA  
  } else {
    rsquared_values <- granpa_rds[["normal_models"]][[1]][["results"]][["Rsquared"]]
    mean(rsquared_values)
  }
})

mean_rsquared_vector <- unlist(mean_rsquared_list)

filenames <- basename(granpa_rds_list)

###creating a data frame

data <- data.frame(
  sample=filenames,
  TF_FDR=TF_peak_FDR_list,
  gene_FDR=as.factor(peak_gene_FDR_list),
  correlation=correlation_list,
  promoter_range=promoter_range_list,
  tf_database=tf_database_list,
  Rsq=mean_rsquared_vector
)

write.csv(data, here(output_folder, "result_summary.csv"), row.names=FALSE)
data <- read.csv("result_summary.csv")

####ggplot2 heatmap

plot <- ggplot(data, aes(x=gene_FDR, y=TF_FDR, fill= Rsq))+ 
  geom_tile()+ 
  labs(x="gene_FDR", y="TF_FDR", fill="Rsq")+
  facet_grid(vars(promoter_range), vars(correlation))

ggsave(here(output_folder,"parameters.pdf"), plot = plot, width = 8, height = 6)
  
###number of unique TFs

#unique_TFs <- Map(function(data_list, filenames) {
#  unique_TFs <- length(unique(data_list[["GRN"]][["TF.name"]]))
#  })
#print(paste("File:", filenames, "- Unique TFs:", unique_TFs))
