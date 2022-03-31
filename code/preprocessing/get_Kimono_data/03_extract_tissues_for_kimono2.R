#!/usr/bin/env Rscript

###############################################
library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Takes expression matrix and phenotypes matrix to subset accoridng to tissue group')
parser$add_argument( '-expr', '--expr_file', metavar='<expression_332_SMTSD.csv>', 
                     help='expression_data')
parser$add_argument( '-pheno', '--pheno_file', metavar='<phenotypes_332_SMTSD.csv>',
                     help='phenotypes data')
parser$add_argument( '-prs', '--prs_file', metavar='<PRS.csv>',
                     help='prs data')
parser$add_argument( '-como', '--como_file', metavar='<como.csv>',
                     help='como data')
parser$add_argument('-o', '--out_prefix')
args <- parser$parse_args()

parser$print_help()

##############################################
### 1 libraries
###############################################

libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify",
               "magrittr", "foreach", "doParallel", "ggthemes", "dplyr", "furrr", "purrr")
lapply(libraries, require, character.only = TRUE)

print("Libraries loaded")

###############################################
### 2 Data
###############################################
expr <- fread(args$expr_file)
pheno <- fread(args$pheno_file)
prs <- fread(args$prs_file)
como <- fread(args$como_file)

print("Files reads")
system(paste0("mkdir -p ", args$out_prefix))
###############################################
### 3 Subset
###############################################
#mytissue="Lung" 
alltissues <- pheno$SMTS %>% unique


make_subset_pheno_expr <- function(mytissue, df_pheno=pheno, df_expr=expr, df_prs=prs, df_como=como){
  # subset phenotypes
  subset_pheno <- df_pheno[SMTS == mytissue,]
  # remove dummy coded columns of SMTSD
  tiss_cols <- colnames(subset_pheno)[grep("SMTSD_",colnames(subset_pheno))]
  to_be_dropped <- !(colnames(subset_pheno) %in% tiss_cols)
  subset_pheno <- subset_pheno[,..to_be_dropped]

  # samples to keep
  sampid <- subset_pheno$SAMPID


  #subset expression
  subset_expr <- df_expr[SAMPID %in% sampid,]
  subset_expr <- subset_expr[match(sampid, subset_expr$SAMPID),]

  #subset como
  subset_prs <- df_prs[SAMPID %in% sampid,]
  subset_prs <- subset_prs[match(sampid, subset_prs$SAMPID),]
  
  #subset prs
  subset_como <- df_como[SAMPID %in% sampid,]
  subset_como <- subset_como[match(sampid, subset_como$SAMPID),]
  
  
  # write files
  print(paste0(args$out_prefix, "expr_biogrid", "_SMTS_",  mytissue, ".csv"))
  fwrite(subset_pheno, file=paste0(args$out_prefix, "pheno_biogrid", "_SMTS_",  mytissue, ".csv"))
  fwrite(subset_expr, file=paste0(args$out_prefix, "expr_biogrid", "_SMTS_",  mytissue, ".csv"))
  fwrite(subset_prs, file=paste0(args$out_prefix, "prs_biogrid", "_SMTS_",  mytissue, ".csv"))
  fwrite(subset_como, file=paste0(args$out_prefix, "como_biogrid", "_SMTS_",  mytissue, ".csv"))
}

print("Starting tissue extraction for phenotypes and expression/n Saving files...")

walk(alltissues,make_subset_pheno_expr)
print("Files written")



