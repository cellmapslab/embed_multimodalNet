# Author: Nathalie Gerstner, Edited by:Ghalia Rehawi
# Identify GWAS with less than 250000 samples
# and evaluate if highly polygenic
# Count number of SNPs with p <= 5e-08 (This info is to be addes in a GWAS summery file GWAS_summaries.tsv)

base_dir <- "/GWAS/"
#Repeat the following for all disease/traits summery statitics files
cov19hgsev <- read.table(paste0(base_dir, "covid19hgsevere/covid19hgsevere.sumstats.tsv"), sep = "\t", header = TRUE)
n_sig <- length(which(cov19hgsev$P <= 5e-08))

print(n_sig)
