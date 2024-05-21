configfile: 'config.yaml'

branch_list = config['branches'].split(' ')
de_list = config['de'].split(' ')
res_list = config['res'].split(' ')
cormethod_list = config['cormethod'].split(' ')
tf_peak_fdr_list = config['tf_peak_fdr'].split(' ')
prom_range_list = config['prom_range'].split(' ')
de_pval_th_list = config['de_pval_th'].split(' ')
peak_gene_fdr_list = config['peak_gene_fdr'].split(' ')

rule all:
    input:
        expand('{branch}/{branch}_res{res}_cormethod_{cormethod}_tf-peak-fdr_{tf_peak_fdr}_prom-range_{prom_range}_de_{de}.qs_de_pval_th_{de_pval_th}_peak_gene_fdr_{peak_gene_fdr}granpa.rds', branch=branch_list, res=res_list, cormethod=cormethod_list, tf_peak_fdr=tf_peak_fdr_list, prom_range=prom_range_list, de=de_list, de_pval_th=de_pval_th_list)

rule get_grn_from_tsv:
    input:
        '{branch}_atac.pseudobulkFromClusters_res{res}_mean.tsv.gz',
        '{branch}_rna.pseudobulkFromClusters_res{res}_mean.tsv.gz',
        '{branch}_metadata_res{res}_mean.tsv.gz'
    conda:
        'r_many_packages'
    output:
        '{branch}/{branch}_res{res}_cormethod_{cormethod}_tf-peak-fdr_{tf_peak_fdr}_peak_gene_fdr_{peak_gene_fdr}_prom-range_{prom_range}_grn.rds'
    shell:
        './granie_parser.R --meta_data {wildcards.branch}_metadata_res{wildcards.res}_mean.tsv.gz --atac {wildcards.branch}_atac.pseudobulkFromClusters_res{wildcards.res}_mean.tsv.gz --rna_seq {wildcards.branch}_rna.pseudobulkFromClusters_res{wildcards.res}_mean.tsv.gz --output_folder {wildcards.branch} --branch {wildcards.branch} --correlation_method {wildcards.cormethod} --promoter_range {wildcards.prom_range} --TF_peak_FDR {wildcards.tf_peak_fdr} --res {wildcards.res} --n_cores 12'

rule get_granpa_from_grn_and_de:
    input:
        '{branch}/{branch}_res{res}_cormethod_{cormethod}_tf-peak-fdr_{tf_peak_fdr}_peak_gene_fdr_{peak_gene_fdr}_prom-range_{prom_range}_grn.rds',
        '{de}.qs'
    conda:
        'r_many_packages'
    output:
        '{branch}/{branch}_res{res}_cormethod_{cormethod}_tf_peak_fdr_{tf_peak_fdr}_peak_gene_fdr_{peak_gene_fdr}_prom-range_{prom_range}_de_{de}.qs_de_pval_th_{de_pval_th}_granpa.rds'
    shell:
        './granpa_parser.R --grn_rds {wildcards.branch}/{wildcards.branch}_res{wildcards.res}_cormethod_{wildcards.cormethod}_tf_peak_fdr_{wildcards.tf_peak_fdr}_peak_gene_fdr_{wildcards.peak_gene_fdr}_prom_range_{wildcards.prom_range}_grn.rds --de_data {wildcards.de}.qs --de_pval_th {wildcards.de_pval_th} --logFC_th 0 --branch {wildcards.branch} --n_cores 12 --importance_tf "permutation" --ml_type "regression" --output_folder_granpa {wildcards.branch}_de_{wildcards.de}_granpa'

rule summarise_results:
    input:
        expand('{branch}_de_{de}_granpa/{branch}_de_{de}_granpa.rds', branch=branch_list, de=de_list)
    output:
        summary_csv = '{branch}_de_{de}_granpa/result_summary.csv',
        parameters_pdf = '{branch}_de_{de}_granpa/parameters.pdf'
    conda:
        'r_many_packages'
    params:
        granpa_rds_list=lambda wildcards: expand('{branch}_de_{de}_granpa/{branch}_de_{de}_granpa.rds', branch=wildcards.branch, de=wildcards.de),
        TF_peak_FDR_list=tf_peak_fdr_list,
        peak_gene_FDR_list=peak_gene_fdr_list,
        correlation_list=cormethod_list,
        promoter_range_list=prom_range_list,
        output_folder=lambda wildcards: f"{wildcards.branch}_de_{wildcards.de}_granpa"
    shell:
        './summarise_results.R --granpa_rds_list {params.granpa_rds_list} --TF_peak_FDR_list {params.TF_peak_FDR_list} '
        '--peak_gene_FDR_list {params.peak_gene_fdr_list} --correlation_list {params.correlation_list} '
        '--promoter_range_list {params.promoter_range_list} '
        '--output_folder {params.output_folder}'
