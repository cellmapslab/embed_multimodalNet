#!/bin/sh
# Author Nathalie Gerstner, Edited by Ghalia Rehawi
# Date 06.09.2021
# This script concatenates files of all chromosomes for each trait into a single file

#list of diseases/traits from summary file
traits=( "covid19hg" "covid19hghospitalized" "covid19hgsevere" )

for i in "${traits[@]}"
do

	cat /calculated_PRS/${i}/* > /calculated_PRS/${i}/${i}_all_chromosomes.txt

done
