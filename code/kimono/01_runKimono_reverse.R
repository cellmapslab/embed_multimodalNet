#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Takes expression matrix, phenotype file, priors to generate network')

parser$add_argument( '-expr', '--expr_file', metavar='<expression.csv>', 
                     help='expression_data')
parser$add_argument( '-pheno', '--pheno_file', metavar='<phenotypes.csv>',
                     help='phenotypes data')
parser$add_argument( '-prs', '--prs_file', metavar='<prs.csv>',
                     help='prs data')
parser$add_argument( '-como', '--como_file', metavar='<como.csv>',
                     help='comorbidities data')
parser$add_argument( '-network', '--network_file', metavar='<network.csv>',
                     help='comorbidities data')

parser$add_argument('-o', '--output')
args <- parser$parse_args()

parser$print_help()

###############################################
### 1 libraries
###############################################

libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify",
               "magrittr", "foreach", "doParallel", "ggthemes", "dplyr", "furrr", "purrr", "stringr")
lapply(libraries, require, character.only = TRUE)

print("libraries done")

###############################################

source("../../code/kimono/kimono/R/infer_sgl_model.R")
source('../../code/kimono/kimono/R/utility_functions.R')
source('../../code/kimono/kimono/R/kimono.R')
print("moni main done")

###############################################
### 2 Data
###############################################

# mrna <- fread(paste0(inputdir,'expression332.csv')) %>% as.data.frame
expr <- fread(paste0(args$expr)) %>% as.data.frame

myvars <- colnames(expr) %in% c("SAMPID", "subid", "SUBJID")
expr <- as.data.table(expr[,!myvars])
expr[1:5,1:5]
print("expr done")
dim(expr)

###############################################
# pheno <- fread(paste0(inputdir,'phenotypes332_SMTS.csv')) %>% as.data.frame
pheno <- fread(args$pheno) %>% as.data.frame
rownames(pheno) <- pheno[,"SAMPID"]
pheno <- as.data.table(pheno[,! (colnames(pheno) %in% c("SAMPID", "SMTS", "SMTSD", "subid", "SUBJID"))])
pheno[1:5,]
dim(pheno)

###############################################
# prs <- fread(paste0(inputdir,'phenotypes332_SMTS.csv')) %>% as.data.frame
prs <- fread(args$prs) %>% as.data.frame
rownames(prs) <- prs[,"SAMPID"]
prs <- as.data.table(prs[,! (colnames(prs) %in% c("SAMPID", "SMTS", "SMTSD", "subid", "SUBJID"))])
prs[1:5,]
dim(prs)

###############################################
#como <- fread("/home/icb/yue.hu/proj/data/kimono/pre/como_biogrid.csv") %>% as.data.frame
como <- fread(args$como) %>% as.data.frame
rownames(como) <- como[,"SAMPID"]
como <- as.data.table(como[,! (colnames(como) %in% c("SAMPID", "SMTS", "SMTSD", "subid", "SUBJID"))])
como[1:5,]
dim(como)


###############################################
# read in network 
network <- fread(args$network)



# run
options(future.globals.maxSize = 15000 * 1024^2)
plan(multisession, workers = 12)
start_time <- Sys.time()
print(paste0("start time:", start_time))

###############################################
###############################################
###############################################
###############################################

###############################################
###############################################
###############################################
# PRS
###############################################
###############################################
###############################################

start_time <- Sys.time()
#set input parameters
input_list <- list(
  as.data.table(prs),
  as.data.table(expr)
)
names(input_list) <- c('prs',
                       'expr')

#########################
mapping_list <- list(
  network[relation=="expr_prs"][,.(target,predictor)]
)
node_list <- network[relation=="expr_prs", predictor] %>% unique

#########################
metainfo <-data.frame('ID'   = c('prs_expr'),
                      'main_to'   =  c(2)
)

print("data created for prs")

###############################################

results_prs <- future_map(node_list, run_kimono_para, sel_iterations = 30, 
.options = furrr::future_options(seed = TRUE), .progress = T) #add one more layer of parallelization in kimono.R (seeds)
results_prs <- do.call(rbind, results_prs) # make one big data frame

###############################################


print("kimono done"); end_time <- Sys.time();end_time - start_time
print(paste0("end time:", Sys.time()))

###############################################
###############################################
# PHENOPTYPE
###############################################
###############################################
###############################################


start_time <- Sys.time()
#set input parameters
input_list <- list(
  as.data.table(pheno),
  #as.data.table(expr[,1:500,with=FALSE])
  as.data.table(expr)
)

names(input_list) <- c('pheno',
                       'expr')

#########################
mapping_list <- list(
  network[relation=="expr_bio"][,.(target,predictor)]
)

foi <- fread("../../data/lookup/genelist.csv", 
  header = FALSE)
node_list <- network[relation=="expr_bio", predictor] %>% unique

node_list <- node_list[node_list%in% foi$V1]
node_list <- node_list[4:length(node_list)]
#########################
metainfo <-data.frame('ID'   = c('bio_expr'),
                      'main_to'   =  c(2)
)
node_list
print("data created for bio")

###############################################

results_bio <- future_map(node_list, run_kimono_para, sel_iterations = 30, 
.options = furrr::future_options(seed = TRUE), .progress = T) #add one more layer of parallelization in kimono.R (seeds)
results_bio <- do.call(rbind, results_bio) # make one big data frame

###############################################

print("kimono done"); end_time <- Sys.time();end_time - start_time
print(paste0("end time:", Sys.time()))

# ###############################################
# ###############################################
# ###############################################
# # COMO
# ###############################################
# ###############################################
# ###############################################

start_time <- Sys.time()
#set input parameters
input_list <- list(
  as.data.table(como),
  as.data.table(expr)
)
names(input_list) <- c('como',
                       'expr')

#########################
mapping_list <- list(
  network[relation=="expr_como"][,.(target,predictor)]
)
node_list <- network[relation=="expr_como", predictor] %>% unique

#########################
metainfo <-data.frame('ID'   = c('como_expr'),
                      'main_to'   =  c(2)
)
head(node_list)

print("data created for como")

###############################################

results_como <- future_map(node_list, run_kimono_para, sel_iterations = 30, 
.options = furrr::future_options(seed = TRUE), .progress = T) #add one more layer of parallelization in kimono.R (seeds)
results_como <- do.call(rbind, results_como) # make one big data frame

###############################################


print("kimono done"); end_time <- Sys.time();end_time - start_time
print(paste0("end time:", Sys.time()))



# ###############################################
# ###############################################
# ###############################################

results_all <- rbind(results_prs, results_como, results_bio )
#end
fwrite(results_all,file=args$output)
print("table written")


