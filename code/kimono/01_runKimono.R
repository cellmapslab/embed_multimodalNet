#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Takes expression matrix, phenotype file, priors to generate network')
parser$add_argument( '-startnode', '--startnode', type="integer",
                     help='the node to start inferring the network, the next 500 gene models will be calculated')

parser$add_argument( '-expr', '--expr_file', metavar='<expression.csv>', 
                     help='expression_data')
parser$add_argument( '-pheno', '--pheno_file', metavar='<phenotypes.csv>',
                     help='phenotypes data')
parser$add_argument( '-prs', '--prs_file', metavar='<prs.csv>',
                     help='prs data')
parser$add_argument( '-como', '--como_file', metavar='<como.csv>',
                     help='comorbidities data')

parser$add_argument( '-pr_pheno', '--prior_pheno', metavar='<prior_expr_bio.csv>', 
                     help='prior phenotypes - gene covariates mapping')
parser$add_argument( '-pr_prs', '--prior_prs', metavar='<prior_expr_bio.csv>', 
                     help='prior PRS - gene covariates mapping')
parser$add_argument( '-pr_como', '--prior_como', metavar='<prior_expr_como.csv>', 
                     help='prior comorbidities - gene covariates mapping')

parser$add_argument( '-pr_expr', '--prior_expr', metavar='<prior_expr.csv>',
                     help='prior expression - gene mapping')
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
# prior_pheno <- fread(paste0(inputdir,'prior_expr_bio.csv'));prior_pheno[1:5,]
prior_pheno <- fread(args$prior_pheno);tail(prior_pheno)
prior_expr <- fread(args$prior_expr); prior_expr[1:5,]
prior_prs <- fread(args$prior_prs); prior_prs[1:5,]
prior_como <- fread(args$prior_como); prior_como[1:5,]


print("read mapping")

###############################################
# 3 Assemble into lists
###############################################

#set input parameters
input_list <- list(
  # as.data.table(expr[,1:5,with=FALSE]),
  as.data.table(expr),
  as.data.table(pheno),
  as.data.table(prs),
   as.data.table(como)

)
names(input_list) <- c('expr',
                       'pheno', 
                       'prs',
                       'como')

#########################
mapping_list <- list(
  as.data.table(prior_expr),
  as.data.table(prior_pheno),
  as.data.table(prior_prs),
  as.data.table(prior_como)
  
)
#########################
metainfo <-data.frame('ID'   = c( 'prior_expr', 'expr_bio', 'expr_prs', 'expr_como'),
                      'main_to'   =  c(1,2,3,4)
)
print("data created")
###############################################
# 4 Run MONI
###############################################


# run
options(future.globals.maxSize = 20000 * 1024^2)
plan(multisession, workers = 46)
start_time <- Sys.time()
print(paste0("start time:", start_time))


startnode=args$startnode
endnode <- min(startnode+498, length(colnames(input_list$expr)))
print(paste0("startnode=", startnode,"; endnode=", endnode))
node_list <- colnames(input_list$expr)[startnode:endnode]
results <- future_map(node_list, run_kimono_para, sel_iterations = 30, .options = furrr::future_options(seed = TRUE), .progress = T) #add one more layer of parallelization in kimono.R (seeds)



###############################################
results <- do.call(rbind, results) # make one big data frame

print("kimono done"); end_time <- Sys.time();end_time - start_time
print(paste0("end time:", Sys.time()))



#end
fwrite(results,file=args$output)
print("table written")


