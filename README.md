This repository contains scripts for GRaNIE and GRaNPA (https://www.bioconductor.org/packages/release/bioc/html/GRaNIE.html) snakemake pipeline.

```
run_granie.R
```

This is a GRaNIE parser script that takes ATAC-seq, RNA-seq and metadata files as input allowing to vary the parameters:

`--meta_data`          Path to TSV with metadata

`--atac`               Path to TSV with ATAC data

`--rna_seq`            Path to RNA seq counts data

`--output_folder`      Folder name for output

`--branch`             Branch name

`--hocomoco_path`      Path to HOCOMOCO folder (Default: PWMScan_HOCOMOCOv12/H12INVIVO)

`--correlation_method` Correlation method in addConnections_peak_gene (Default: pearson)

`--promoter_range`     Promoter range in addConnections_peak_gene, bp (Default: 250000)

`--TF_peak_FDR`        FDR threshold in filterGRNAndConnectGenes (Default: 0.2)

`--peak_gene_FDR`      Peak gene FDR threshold in filterGRNAndConnectGenes (Default: 0.2)

`--n_cores`            Number of cores in addConnections_peak_gene and overlapPeaksAndTFBS (Default: 36)

`--res`                Clustering resolution

### Example

```bash
Rscript run_granie.R --atac data/atac.tsv --rna_seq data/rna_seq_counts.tsv --output_folder results_granie --branch B_cell --res 0.5
```

```
run_granpa.R
```

This is a GRaNPA parser script that takes the `.rds` file from `run_granie.R` as input. For more info on GRaNPA: https://grp-zaugg.embl-community.io/GRaNPA/

`--grn_rds`              Path to GRaNIE rds (Default: granie_output)

`--de_data`              Path to differential expression data for validation (Default: folder_name)

`--de_pval_th`           P-value threshold for differential expression (Default: 0.2)

`--logFC_th`             LogFC threshold for differential expression (Default: 0)

`--n_cores`              Number of cores (Default: 10)

`--importance_tf`        Algorithm for finding most important TFs (Default: permutation)

`--ml_type`              Regression or classification (Default: regression)

`--output_folder_granpa` Output folder name

`--branch`               Branch name

### Example

```bash
Rscript run_granpa.R --grn_rds results_granie/granie_output.rds --de_data data/differential_expression.tsv --output_folder_granpa results --branch B_cell
```


# Snakemake pipeline

This snakemake pipeline includes `run_granie.R`, `run_granpa.R` and `summarise_results.R` to process multiple files with different parameters and creates a heatmap that shows resulting R^2 values. The parameters are provided in `config2.yaml` file.
 
### Example config file

```yaml
branches:
traj1
de:
de1
de3
res:
0.5
2
5
10
15
20'
tf_peak_fdr:
0.05
0.1
0.2
0.3
prom_range:
50000
100000
1000000
cormethod:
pearson
spearman
bicor
```



