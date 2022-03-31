#!/bin/bash
# Author: Nathalie Gerstner, Edited by Ghalia Rehawi
# Date: 09.09.2021
# This script is to run PLINK2 to calculate overall risk of each individual in the GTEx cohort for the different diseases/traits

traits=( "covid19hg" "covid19hghospitalized" "covid19hgsevere" ) 

for n in "${traits[@]}"
do

	plink2 --score "/calculated_PRS/${n}/mapped_${n}_all_chromosomes.txt" 2 4 6 --allow-extra-chr --extract "SNPs_list.txt" --out "/calculated_PRS/${n}/${n}_plink_out"  --bfile "/GTEX-lifted" 
done

