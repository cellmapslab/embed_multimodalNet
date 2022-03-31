#!/usr/bin/env Rscript

###############################################
library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Assemble one big network out of the many files; filter abs(beta) >0.01 and mean_rsq > 0.1')
parser$add_argument( '-i', '--input_dir',
                     help='input directory where network_GTex_startnode1.tsv files are stored')
args <- parser$parse_args()

parser$print_help()

##############################################
###############################################
### 1 libraries
###############################################

libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify",
               "magrittr", "foreach", "doParallel", "ggthemes", "dplyr", "furrr", "purrr")
lapply(libraries, require, character.only = TRUE)
print("Libraries loaded")

# ###############################################
### 2 libraries
###############################################

abs_files <- system(paste0("ls -d ", args$input_dir, "/*.tsv"), intern=T)

# create empty data table
y <- data.table(
  target = character(),
  predictor = character(),
  relation = character(),
  sel_freq = numeric(),
  mean_value = numeric(),
  sd_value = integer(),
  mean_rsq = numeric(),
  sd_rsq = numeric(),
  mean_mse = numeric(),
  sd_mse = numeric()
)


# read in all over loop and concat to y
for (i in 1:length(abs_files)){
  x <- fread(abs_files[i])
  y <- rbind(y,x)
}
y <- y %>% unique
print(paste0("number of rows: ",nrow(y)))


# remove beta=0
network_raw <- y[mean_value!=0,]

print(paste0("number of rows after removing beta = 0: ",nrow(network_raw)))

# get a normalized beta - corrected for the number of features that each gene model has
network_raw[,id:=1:nrow(network_raw)]
network_raw[,features:=.N, by=target]
network_raw[, beta:=(mean_value*features)]
fwrite(network_raw, file=paste0(args$input_dir, "/full_network_nofilter.csv"))

# filtering on selected freq > 0.7; abs(beta) > 0.01; rsq > 0.1
network <- network_raw %>% 
  filter(sel_freq > 0.7) %>%
  filter((beta > 0.01) | (beta < (-0.01))) %>%
  filter(mean_rsq > 0.1) %>%
  filter(predictor != '(Intercept)') %>% setDT

print(paste0("number of rows after filtering: ",nrow(network)))

print(paste0("number of models: ", length(unique(y$target))))
print(paste0("number of models before filtering: ", length(unique(network_raw$target))))
print(paste0("number of models aftr filtering: ",
  length(unique(network$target))))
fwrite(network, file=paste0(args$input_dir, "/full_network_filter.csv"))
print("written network file")

per_fil_sparse <- nrow(network_raw)/nrow(y)
per_fil <- nrow(network)/nrow(network_raw)
per_fil_all <- nrow(network)/nrow(y)

paste0("Percentage of edges retained from introduction of sparsity: ", 100*round(per_fil_sparse, digit=3), " %")
paste0("Percentage of edges retained after filtering for performance > 0.1 and normalized for number of retained features: ", 100*round(per_fil, digit=3), " %")
paste0("Percentage of edges retained from all possible edges: ", 100*round(per_fil_all, digit=3), " %")
