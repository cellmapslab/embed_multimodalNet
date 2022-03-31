#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Takes expression matrix, phenotype file, priors to generate network')
parser$add_argument( '--node_anno', 
                     help='node annotation')

parser$add_argument( '--raw_data',
                     help='most similar and least similar from 100 runs')

parser$add_argument('--outputfile')
args <- parser$parse_args()

parser$print_help()

##################################
# libraries
libraries <- c("data.table", "tidyr", "ggplot2", "magrittr", "ggthemes", "dplyr", "igraph", "RColorBrewer", "gridExtra", "grid", "purrr")
lapply(libraries, require, character.only = TRUE)

##################################
## data

node_anno <-  fread(args$node_anno)


all <- fread(args$raw_data)

colnames(all) <- c("V1","node1","node2","sim","neighbour","run")
all[,V1:=NULL]
all$node1 <- gsub("GO:", "", all$node1)
all$node2 <- gsub("GO:", "", all$node2)
library(stringr)
colnames(node_anno) <- c("type", "node")

head(all)
##################################
## merge with node annotation, count how many times gene appears for trait of interest in the top/bottom 1000, calculate max and mean and sd of similarity score 



z <- merge(all, node_anno, by.x="node2", by.y="node", all.x=T) %>% dplyr::rename(type_node2=type)
res_full <- merge(z, node_anno, by.x="node1", by.y="node", all.x=T) %>% dplyr::rename(type_node1=type)
res_full <- res_full[,.(node1, node2, type_node1, type_node2, sim, neighbour, run)]
rm(z)

head(res_full)

# get max, mean and sd
res_full <- res_full[order(run,-sim)]
max_sim_score <- res_full[, max(sim), by = .(node2, node1)] %>% dplyr::rename(max_sim := V1)
mean_sim_score <-  res_full[, mean(sim), by = .(node2, node1)] %>% dplyr::rename(mean_sim := V1)
var_sim <- res_full[, var(sim), by = .(node2, node1)] %>% dplyr::rename(var_sim := V1)

head(max_sim_score);head(mean_sim_score);head(var_sim)

y <- res_full[, .N, by = .(node1, node2, neighbour, type_node1, type_node2)]
y <- merge(y, max_sim_score, by = c("node2", "node1"))
y <- merge(y, mean_sim_score, by = c("node2", "node1"))
y <- merge(y, var_sim,  by = c("node2", "node1"))
result <- y[order(node1,-mean_sim)]

head(y)

##################################
# ## robust results - found over 80 % of the runs


num_robust <- 0.8 * max(result$N)
result_stable <- result[N>=num_robust,]

##################################
## results & save

head(result_stable)
fwrite(result_stable, file=args$outputfile)