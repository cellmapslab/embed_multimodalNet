#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Takes expression matrix, phenotype file, priors to generate network')
parser$add_argument( '-expr', '--expr_file', metavar='<expression.gct>', 
                     help='expression_data')
parser$add_argument( '-anno', '--anno_file', metavar='<annotation.csv>',
                     help='phenotypes data')
parser$add_argument('-o', '--output')
args <- parser$parse_args()

parser$print_help()

###############################################
###############################################
###############################################

libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify",
               "magrittr", "foreach", "doParallel", "ggthemes", "dplyr", "furrr", "purrr")
lapply(libraries, require, character.only = TRUE)

###############################################
###############################################
###############################################

expr <- fread(args$expr_file)
print("read in expression data")

a <- fread(args$anno_file)
print("read in annotation data")


###############################################
###############################################
###############################################



# vector with only protein coding genes not on chr X, Y and M
tokeep <- a[type=="gene"][
  gene_type=="protein_coding"][
  !chr %in% c("chrM", "chrY", "chrX"),][,
                                        gene_id]


# remove MXY from expression data
expr <- expr[Name%in%tokeep]
print("removed chr x,y,m")
cat("\t", "nrow: ", nrow(expr),", ncol: ", ncol(expr), "\n")



###############################################
###############################################
###############################################

# remove sex tissues c("Ovary", "Uterus", "Vagina", "Prostate", "Testis", "Cervix_Ectocervix", "Fallopian_Tube", "Cervix_Endocervix", "Cervix_Ectocervix")
nosextissue_samples <- fread('../../../data/phenotypes/input/nosextissue_samples.csv')
tokeep2 <- c("Name", "Description", nosextissue_samples$SAMPID %>% unique)

tokeep2 <- gsub("\\.","-",tokeep2)

expr <- expr[,..tokeep2]
print("removing sex tissues")
cat("\t", "these are:", "Ovary", "Uterus", "Vagina", "Prostate", "Testis",
      "Cervix_Ectocervix", "Fallopian_Tube", "Cervix_Endocervix", "Cervix_Ectocervix", "\n")
cat("\t", "nrow: ", nrow(expr),", ncol: ", ncol(expr), "\n")

###############################################
###############################################
###############################################

# remove low expression genes
# at least expression of 0.1 in 80 % of the samples

remove_low_expr <- function(threshold=0.1, percentage=0.8, dt=expr_keeps){
  count <-apply( select(dt,-c("Name", "Description"))>0.1, 1, sum) # 0.1 expression threshold
  crit <- nrow(dt)*0.8 # in 80 % of the samples
  return(dt[count >= crit])
}

print("removing low expr")
expr <- remove_low_expr(dt=expr)
cat("\t", "nrow: ", nrow(expr),", ncol: ", ncol(expr), "\n")
print("done")

fwrite(expr, 
       file=args$output)


