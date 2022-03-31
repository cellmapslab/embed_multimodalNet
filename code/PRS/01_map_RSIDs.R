#Author: Ghalia Rehawi
#Date: 27.08.2021
#This script maps the SNPs IDs from format CHR_POS to RSXXXXX 

#Load libraries
library(dplyr)
#library(tidyr)
library(data.table)

setwd('...')
#Read target data SNP file 
snps = fread('GTEX-lifted.bim')

#Read the LD reference snp info file
LD_ref = fread('ldblk_1kg_eur/snpinfo_1kg_hm3')
LD_ref = LD_ref[, CHR:=as.character(CHR)]

joined_dt = merge(snps, LD_ref, by.x = c("V1", "V4"), by.y = c("CHR", "BP"), all.x = FALSE, all.y = FALSE)
print(joined_dt[1:10,])
print(dim(joined_dt))
new_bim = joined_dt[ ,.(V1, SNP, V3, V4, V5, V6)]
fwrite(new_bim, 'mapped_GTEx_SNPs.bim', sep= "\t")
