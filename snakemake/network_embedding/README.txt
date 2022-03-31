# prepare input
nohup snakemake --use-conda --cores 2 --configfile config.yml  > 2022-01-16_prepare.out &

nohup snakemake --use-conda --cores 2 --configfile config.yml  > 2022-01-16_dc_filt_genewalk_nohup.out &
