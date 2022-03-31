#!/bin/bash
# Author:Ghalia Rehawi
# Date 02.09.2021
# This script calls the PRC-CS tool on selected traits/diseases 

#GWAS_summaries: a summery of all GWAS ss that we have, it includes information
#like the name of the disease/trait, number of samples, polygenecity and sig number of SNPs
traits=( $(cut -f8 /GWAS/GWAS_summaries.tsv | tail -n +2) )
nsamples=( $(cut -f3 /GWAS/GWAS_summaries.tsv | tail -n +2) )
phicode=( $(cut -f10 /GWAS/GWAS_summaries.tsv | tail -n +2) )

selected_traits=( "covid19hg" "covid19hghospitalized" "covid19hgsevere" )

#create directory for storing results of each trait
for (( n=0; n<${#traits[@]}; n++ )); do

if [[ " ${selected_traits[*]} " == *" ${traits[$n]} "* ]]; then

mkdir -p "/calculated_PRS/${traits[$n]}"

#run PRS-CS
if [ ${phicode[$n]} == 2 ]; then

	python3 /PRScs-master/PRScs.py --ref_dir=/ldblk_1kg_eur --bim_prefix="mapped_GTExSNPs.nodup" --sst_file="/GWAS/${traits[$n]}/${traits[$n]}.sumstats.tsv" --n_gwas=${nsamples[$n]} --out_dir="/calculated_PRS/${traits[$n]}/${traits[$n]}"

elif [ ${phicode[$n]} == 1 ]; then

	python3 /PRScs-master/PRScs.py --ref_dir=/ldblk_1kg_eur --bim_prefix="mapped_GTExSNPs.nodup" --sst_file="/GWAS/${traits[$n]}/${traits[$n]}.sumstats.tsv" --n_gwas=${nsamples[$n]}  --phi=1e-2 --out_dir="/calculated_PRS/${traits[$n]}/${traits[$n]}"

else
	python3 /PRScs-master/PRScs.py --ref_dir=/ldblk_1kg_eur --bim_prefix="mapped_GTExSNPs.nodup" --sst_file="/GWAS/${traits[$n]}/${traits[$n]}.sumstats.tsv" --n_gwas=${nsamples[$n]} --phi=1e-4  --out_dir="/calculated_PRS/${traits[$n]}/${traits[$n]}"
fi
fi
done
