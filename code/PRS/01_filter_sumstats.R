##################################################
## Date: 11.10.2021
## Author: Nathalie Gerstner, edited by Ghalia Rehawi
##################################################
# Filter SNPs with NAs out of sumstats 
# Flip SNPs with negative BETA

### INFO: flip step is useless, can be skipped
dir <- "directory/to/GWAS-ss/"


### FUNCTIONS ############################

filterRisk <- function(risks) # filter out SNPs with NA in BETA or P
  {
  risks <- na.omit(risks, cols = c("BETA", "P"))
  return(risks)
  }

filterID <- function(risks) # remove variants without rs ID and cut suffices like :G:A
  {
  snp_id <- as.character(risks$SNP)
  missing_rs <- which(!startsWith(snp_id, "rs"))
  if (length(missing_rs)>0) risks <- risks[-missing_rs,]
  risks$SNP <- gsub(":[A-Z]:[A-Z]", "", risks$SNP)
  return(risks)
  }


### FILTER #############################
#This is an example on one summery statistic of a disease
#It should be carried out on all available GWAS SS
psoriasis <- read.table(paste0(dir, "psoriasis.sumstats.tsv"), header = TRUE)

psoriasis <- filterRisk(psoriasis)
psoriasis <- filterID(psoriasis)

write.table(psoriasis,  paste0(dir, "GWAS/psoriasis/", "psoriasis.sumstats.tsv"),r=F,qu=F,sep="\t")

