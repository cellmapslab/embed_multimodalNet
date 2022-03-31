#Author: Ghalia Rehawi
#Date: 10.09.2021
#This script maps the SNPs IDs in the bim file from format RSXXXXX to CHR_POS 

#Load libraries
library(dplyr)
#library(tidyr)
library(data.table)

setwd('...')
#Read mapped SNP file 
mapped_snps = fread('mapped_GTExSNPs.nodup.bim')
mapped_snps = mapped_snps[, V1:=as.character(V1)]

#Read the lifted bim file containing the old format of SNPs
old_snps = fread('GTEX-lifted.bim')

joined_dt = merge(mapped_snps, old_snps, by.x = c("V1", "V4"), by.y = c("V1", "V4"), all.x = FALSE, all.y = FALSE)
print(joined_dt[1:10,])
print(dim(joined_dt))
new_bim = joined_dt[ ,.(V1, V2.x, V2.y)]
fwrite(new_bim, 'remapping_file.bim', sep= "\t")

