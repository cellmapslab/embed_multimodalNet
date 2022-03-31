#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Creates biogrid gene-gene mapping')
parser$add_argument('-i', '--input')
parser$add_argument('-o', '--outprefix')
args <- parser$parse_args()

parser$print_help()

resultsdir <- args$outprefix
##############################################################################################
##############################################################################################
print("workspace and libraries")

# libraries
libraries <- c("data.table", "tidyr", "ggplot2", 
               "reshape2", "oem", "ramify", "magrittr",
               "grpreg", "foreach", "doParallel", "ggthemes",
               "dplyr","rjson")
lapply(libraries, require, character.only = TRUE)



##############################################################################################
##############################################################################################

# MAPPING EXPR-EXPR - biogrid


load(args$input)
prior_expression_biogrid <- biogrid@data %>% setDT
prior_expression_biogrid<-prior_expression_biogrid %>% unique

tmp <- prior_expression_biogrid[,.(gene.b, gene.a)]
colnames(tmp) <- c("gene.a", "gene.b")
prior_expression_biogrid <- rbind(prior_expression_biogrid,tmp) %>% unique

prior_expression_biogrid[1:5,]

fwrite(prior_expression_biogrid, file=paste0(resultsdir, "prior_expression_biogrid.csv"))

