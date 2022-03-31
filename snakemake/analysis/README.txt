# to get top 1000 most_simialar:
nohup snakemake --use-conda --cores 48 --configfiles config.yaml > 2022-01-21_all.out & 

# to get txt embedding:
nohup snakemake --use-conda --cores 1 > 2022-01-10_txt_emb.out &

# to get robust genes:
nohup snakemake --use-conda --cores 1 > 2022-01-11_robust_genes.out &

# to validate tissues:
nohup snakemake --use-conda --cores 1 > 2022-01-11_tissues.out &

# to plot embedding:
nohup snakemake --use-conda --cores 1 > 2022-01-11_plt_embedding.out &

# to plot gene_ranking:
nohup snakemake --use-conda --cores 1 --configfile config.yaml > 2022-01-18_gene_ranking.out &
