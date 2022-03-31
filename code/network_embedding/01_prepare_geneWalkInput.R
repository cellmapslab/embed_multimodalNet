#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Takes output of kimono (merged file) and returns the adjlist')
parser$add_argument('-i', '--input', help='output of Kimono')
parser$add_argument('-o', '--outputdir')
args <- parser$parse_args()

parser$print_help()

###############################################
### 1 libraries
###############################################

libraries <- c("data.table", "tidyr", "ggplot2", "reshape2", "oem", "ramify",
               "magrittr", "foreach", "doParallel", "ggthemes", "dplyr", "furrr", "purrr", "stringr", "igraph")
lapply(libraries, require, character.only = TRUE)

print("libraries done")

###############################################
### 2 make undirected, simple graph
###############################################

network <- fread(args$input)


graph <- graph_from_edgelist(as.matrix(network[,.(target, predictor)])); 
graph <- as.undirected(graph)
graph <- simplify(graph)

edgelist <- as.data.table(as_edgelist(graph))
edgelist[,relation:="some_relation"]

###############################################
### 3 get features of interest
###############################################


phenotypes <- network[relation=="bio_expr", predictor] %>% unique
prs <- network[relation=="prs_expr", predictor] %>% unique
como <- network[relation=="como_expr", predictor] %>% unique
tiss <- network[relation=="tiss_expr", predictor] %>% unique

f_oI <- c(phenotypes, prs, como, tiss)




###############################################
### 4 write files
###############################################
edgelist[,V3:="{}"]
fwrite(edgelist[,.(V1, V2,V3)],paste0(args$outputdir, "/network.edgelist"), col.names = F, sep=" ")


# write custom network file
edgelist[!(V2%in%f_oI), V2:=paste0("GO:", V2)]
edgelist[!(V1%in%f_oI),V1:=paste0("GO:", V1)]
fwrite(edgelist[,.(V1,relation, V2)],paste0(args$outputdir, "/network_sif_full.sif"), col.names = F)

# write gene list
genes <- data.frame(genes=f_oI) %>% setDT
fwrite(genes,paste0(args$outputdir, "/genelist.csv"), col.names = F)

# write my costum go, consisting of only the go terms (aka normal genes)
allnodes <- unique(c(edgelist$V1, edgelist$V2))
mygoterms <- allnodes[grep("GO:",allnodes)]

setwd(args$outputdir)
logFile = "go.obo"
for(i in 1:length(mygoterms)){
  cat("[Term]",
      "\nid:", as.character(mygoterms[i]), 
      "\nname:",gsub("GO:", "", as.character(mygoterms[i])),
      "\nnamespace: gene\n\n",
      file=logFile, append=T)
}

###############################################
### 5 fin
###############################################

