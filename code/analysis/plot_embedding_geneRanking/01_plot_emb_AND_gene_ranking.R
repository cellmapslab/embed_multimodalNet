#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Takes expression matrix, phenotype file, priors to generate network')

parser$add_argument( '--most_similar', 
                     help='node annotation')
parser$add_argument( '--robust_genes',
                     help='robust genes')
parser$add_argument( '--emb_vec',
                     help='embedding vectors')

parser$add_argument('--outputfile')
args <- parser$parse_args()

parser$print_help()

##################################
# libraries
libraries <- c("data.table", "tidyr", "ggplot2", "magrittr", 
               "ggthemes", "dplyr", "igraph", "RColorBrewer", "gridExtra", "grid", "purrr",
               "ggpubr", "ggrepel", "FactoMineR", "factoextra")
lapply(libraries, require, character.only = TRUE)


print("library loaded")


##################################


# code
source('../../../code/analysis/plot_embedding_geneRanking/01_add_fun.R')
source('../../../code/general/add_fun.R')

print("additional functions loaded")
# data
colors <- fread( file="../../../data/lookup/color_code.txt")
node_anno <-  fread( file="../../../data/lookup/node_anno.txt")
traits <- fread("../../../data/lookup/genelist.csv", header=F)

all <- fread(args$most_similar)
results <- fread(args$robust_genes)

result_genes <- results[type_node2=="gene"]

# ##################################
# # PCA
i <- sub(".*_(\\d+).*", "\\1", args$emb_vec)
vectors <- fread(args$emb_vec, skip=1)
vectors$V1 <- gsub("GO:","",vectors$V1)

print("all data read")
##################################
vec_type <- merge(node_anno, vectors, by.x="name", by.y="V1", all.y=T)

##################################
PCA
print("starting with PCA")

col_num <- which(colnames(vec_type) %in% c("name", "type"))
res.pca2 <- PCA(X=vec_type, scale.unit = F, quali.sup =col_num, graph=F)
fviz_pca_ind(X=res.pca2, 
             habillage="type",
             label = "none", addEllipses = T)
eigenvalues <- res.pca2$eig
head(eigenvalues[, 1:2])

pca_data <- as.data.table(cbind(as.character(res.pca2$call$X$name),res.pca2$ind$coord))
pca_data_sub <- pca_data[,.(V1, Dim.1, Dim.2, Dim.3)]
pca_data_sub[,Dim.1:=as.numeric(Dim.1)]
pca_data_sub[,Dim.2:=as.numeric(Dim.2)]
pca_data_sub[,Dim.3:=as.numeric(Dim.3)]
pca_data_sub$V1 <- gsub("name_", "", pca_data_sub$V1)
head(pca_data_sub)

#pca_data_sub <- fread(file="../../../plots/data_dir/pca.txt")
print("plot pca for all traits")

##################################
# plot


pdf(file=paste0(args$outputfile), height=6, width=12)

map(traits$V1, get_plot)

dev.off()
