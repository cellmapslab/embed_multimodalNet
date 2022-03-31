#!/usr/bin/env Rscript

###############################################

library(argparse)

# read in arguments
args = commandArgs(trailingOnly=TRUE)

parser = ArgumentParser(description='Creates data for Kimono; get expr, phenotypes into right format')

parser$add_argument( '-s', '--which_subset', default="biogrid",type ="character",
                     help='')

parser$add_argument('-o', '--outprefix')
args <- parser$parse_args()

parser$print_help()



which_subset <- gsub(" ", "", as.character(args$which_subset))
resultsdir <- args$outprefix

print(which_subset)


which_subset="biogrid"
##############################################################################################
##############################################################################################
print("workspace and libraries")


# libraries
libraries <- c("data.table", "tidyr", "ggplot2",
               "reshape2", "oem", "ramify", "magrittr", "grpreg", 
               "foreach", "doParallel", "ggthemes", "plyr", "dplyr", "purrr", "fastDummies")
lapply(libraries, require, character.only = TRUE)

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
print("Phenotypes from openly available resources")
sample_attributes <- fread("scratch/data/2019_GTEx_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
sample_phenotype <- fread("scratch/data/2019_GTEx_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# add subject ID to sample_attributes
sub <- strsplit(sample_attributes$SAMPID, "-")
subid <- paste0("GTEX-",map(sub,2))
sample_attributes <- cbind(subid, sample_attributes)

##############################################################################################

# merge the two attribute files
sample_attributes <- as.data.table(sample_attributes)
colstoget <- c("subid", "SAMPID", "SMTS", "SMTSD")
y <- merge(sample_attributes[,..colstoget],
               sample_phenotype, by.x="subid", by.y="SUBJID")
y$SMTS <- y$SMTS %>% gsub(" - ", "_",.) %>% gsub(" ", "_",.) %>% gsub("\\(", "",.) %>% gsub("\\)", "",.)
y$SMTSD <- y$SMTSD %>% gsub(" - ", "_",.) %>% gsub(" ", "_",.) %>% gsub("\\(", "",.) %>% gsub("\\)", "",.)

# remove AGE + SEX (as already exist in phenotype file (age with better resolution))
y[,SEX:=NULL];y[,AGE:=NULL]


y[1:5]
##############################################################################################
##############################################################################################
print("Add in extended phenotypes")

y2 <- fread("../../..//data/phenotypes/subsetPhenotypes_ext.csv")
pheno_raw <- merge(y2, y, by.y="subid", by.x="SUBJID")


# exclude sex tissues
tiss_toexclude_SMTSD <- c("Testis", "Cervix_Ectocervix",
                      "Cervix_Endocervix", "Kidney_Medulla",
                      "Bladder", "Fallopian_Tube",
                      "Cells_Leukemia_cell_line_CML", 
                      "Ovary", "Uterus", "Vagina", "Prostate")

pheno_raw <- pheno_raw[! (SMTSD%in%tiss_toexclude_SMTSD),]
print("Sex tissues excluded")


pheno_raw[1:5]
pheno_raw <- as.data.table(as.data.frame(apply(pheno_raw,2,function(x)gsub('\\s+', '',x))))

pheno_col <- c("SAMPID","SUBJID","SEX", "AGE", "BMI", "RACE", "COHORT", "DTHHRDY","SMTS", "SMTSD")
como_col <- c("SAMPID","SUBJID","MHRNLFLR", "MHCOPD", "MHCVD","MHHRTDIS", "MHHTN", "MHLVRDIS",  "MHT1D", "MHT2D", "MHPNMNIA","MHPNMIAB", "MHASTHMA", "MHDPRSSN")
pheno <- pheno_raw[,..pheno_col]
como <- pheno_raw[,..como_col]



##############################################################################################
##############################################################################################

# add the covariates to PHENOTYPE

print("Adding covariates to phenotype")


covariates_dir = "scratch/data/2019_GTEx_v8/GTEx_Analysis_v8_eQTL_covariates"
sample_attributes <- sample_attributes %>% dplyr::rename(SUBJID=subid)


abs_files <- system(paste0("ls -d ", covariates_dir, "/*"), intern=T)


x=abs_files[1]
create_phenotypes_with_covariates <- function(x){
  print(x)
  x_tmp <- strsplit(x, "\\.")  %>%  map(., 1) %>% unlist %>% strsplit(., "/")
  tissue <-  map(x_tmp, length(unlist(x_tmp))) %>% unlist
  print(tissue)
  if(!tissue %in% tiss_toexclude_SMTSD){
    cova <- fread(x)
    mysamples <- colnames(cova)
    # mysamples <- c("GTEX-1117F", "GTEX-ZZPU")
    newdf <- merge(pheno, sample_attributes[,.(SUBJID, SAMPID)], by=c("SAMPID", "SUBJID"), all.x=T)

    newdf <- newdf[SUBJID %in% mysamples][SMTSD%in%tissue,]
    
    samples_cova <- colnames(cova)[-1]
    cols_to_keep <- c("ID", unique(newdf$SUBJID))
    new_cova <- cova[,..cols_to_keep]
    # row_to_keep <- c("PC1", "PC2", "PC3", "PC4", "PC5", "pcr", "platform", "sex")
    # result_df <- new_cova[ID%in%row_to_keep]
    result_df <- new_cova
    cova_interest <- result_df$ID
    result_df_t <- t(result_df[,-1]) %>% data.frame %>% data.table(keep.rownames = T)
    
    colnames(result_df_t)  <- c("SUBJID",cova_interest)
    
    result_df_t_m <- merge(newdf,result_df_t,by="SUBJID")
    return(result_df_t_m)
  }
}

#create_phenotypes_with_covariates(abs_files[1])
pheno_cova <- map(abs_files,create_phenotypes_with_covariates)
# some tissues have 60, 45, 30 and 15 PEER factors
## fill: fills up the empty columns
## set NA to 0


phenotypes <- rbind.fill(pheno_cova) %>% setDT
phenotypes[is.na(phenotypes)] <- 0
phenotypes[,sex:=NULL]
phenotypes[1:5,1:5]
################################################################################
################################################################################
# PHENOTYPES
# dummy code the tissue
phenotypes <- dummy_cols(phenotypes, select_columns = 'SMTSD')
phenotypes[,SMTSD_Cells_Cultured_fibroblasts:=NULL]

# dummy code DTHHRDY
phenotypes <- dummy_cols(phenotypes, select_columns = 'DTHHRDY')
phenotypes[,DTHHRDY_0:=NULL]
phenotypes[,DTHHRDY:=NULL]

# dummy code DTHHRDY
phenotypes <- dummy_cols(phenotypes, select_columns = 'RACE')
phenotypes[,RACE_0:=NULL]
phenotypes[,RACE:=NULL]

# dummy code COHORT
phenotypes <- dummy_cols(phenotypes, select_columns = 'COHORT')
phenotypes[,COHORT_0:=NULL]
phenotypes[,COHORT:=NULL]


##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

# PREPARE EXPRESSION

print("Starting with expression preparation")


expr_raw <- fread("scratch/data/2019_GTEx_v8/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm_filtered.csv")

if(which_subset=="pcnet"){
  # all genes
  prior_expression_sym <- fread(paste0(resultsdir, "prior_expression_pcnet.csv"))
  mappable_genes <- c("CTSL", "ACE2", "TMPRSS2", prior_expression_sym$ILMN.A, 
                      prior_expression_sym$ILMN.B) %>% unique
  expr <- expr_raw[Description %in% mappable_genes,]
  
}else if(which_subset=="krogan"){
  interactions <- fread(paste0("SARS-CoV-2 Host-Pathogen Interaction Map (Fig. 3)).nw_NM20200421.txt"),
                        skip=16) %>% rename(humanprot=`Human Protein`, virusprot=`Virus Protein`)
  mygenes <- c("CTSL", "ACE2", "TMPRSS2", interactions$humanprot)
  expr <- expr_raw[Description %in% mygenes]
  
}else if(which_subset=="biogrid"){
  load("../../../data/phenotypes/input/biogrid.RData")
  # load("/storage/groups/ccm01/workspace/Benchmark_MONI_MDD/Benchmarking/biogrid.RData")
  prior_expression_biogrid <- biogrid@data %>% setDT
  prior_expression_biogrid <- prior_expression_biogrid %>% unique
  mappable_genes <- c(prior_expression_biogrid$gene.a, prior_expression_biogrid$gene.b) %>% unique
  expr <- expr_raw[Description %in% mappable_genes,]
}


# there should be no duplicated genes (genes with multiple ENS ID- Remove them for now)
expr <- expr[!duplicated(expr$Description),]

### Transpose the data frame
# first remember the names
genes <- expr$Description

# expr <- expr[1:10,1:10]
# genes <- genes[1:10]
# transpose all but the first column (name)
t_expr <- as.data.frame(t(expr[,-c(1:2)]))
colnames(t_expr) <- genes
t_expr <- cbind(SAMPID=as.character(row.names(t_expr)), t_expr) %>% setDT
t_expr[, SAMPID := as.character(SAMPID)]


################################################################################
################################################################################
# take intersect of all data
all_sampid <- intersect(intersect(t_expr$SAMPID, phenotypes$SAMPID), como$SAMPID)
t_expr <- t_expr[SAMPID %in% all_sampid,]
como <- como[SAMPID %in% all_sampid,]
phenotypes <- phenotypes[SAMPID %in% all_sampid,]

################################################################################
# order according to idorder
idorder <- as.character(phenotypes$SAMPID)
t_expr <- t_expr[match(idorder, t_expr$SAMPID),]
como <- como[match(idorder, como$SAMPID),]
phenotypes <- phenotypes[match(idorder, phenotypes$SAMPID),]
################################################################################
# order sample_attributes too
sample_attributes <- sample_attributes[match(idorder, sample_attributes$SAMPID),]
#sample_attributes <- sample_attributes %>% dplyr::rename(SUBJID=subid)
fwrite(sample_attributes, file=paste0(resultsdir, "sample_attributes_", which_subset,".csv"))

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
print("Add in PRS")

# PRS directory
prs_dir = "scratch/data/GTEx_PRS/calculated_PRS"
prs_files <- system(paste0("find ", prs_dir," -type f -name *.sscore"), intern=T)

# phenotype directory
p <- phenotypes[,.(SUBJID, SAMPID)]#

p[, names(p) := lapply(.SD, setattr, "label", NULL)]
p[, names(p) := lapply(.SD, setattr, "sas format", NULL)]
idorder <- p$SAMPID

create_PRS_data <- function(x){
  x_tmp <- strsplit(x, "\\.")  %>%  map(., 1) %>% unlist %>% strsplit(., "/")
  trait <-  map(x_tmp, length(unlist(x_tmp))) %>% unlist
  trait <- gsub("_plink_out", "", trait)
  print(paste0(trait, ": ", x))

  y <- fread(x)
  res <- merge(p, y[,.(IID, SCORE1_AVG)], by.x="SUBJID", by.y="IID")
  res <- res[match(idorder, res$SAMPID),]
  res <- res[,.(SAMPID,SCORE1_AVG)]
  res <- res %>% dplyr::rename(!!trait := SCORE1_AVG)
  return(res)
}

prs_data <- map(prs_files,create_PRS_data)
prs <- prs_data%>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="SAMPID"), .)
toxclude <- c("suicide", "EA", "SESA", "CP", "ASD", "AD",
  "age", "gender", "covid19", "covid19severe", "covid19hospital")
tokeep <- colnames(prs)[!(colnames(prs)%in%toxclude)]

prs <- prs[,..tokeep]
#prs <- prs %>% dplyr::select(sort(names(.)))
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
# order according to idorder (both phenotypes and expression AGAIN because of subset for phenotypes)
# order prs too

idorder <- as.character(phenotypes$SAMPID)
t_expr <- t_expr[match(idorder, t_expr$SAMPID),]
como <- como[match(idorder, como$SAMPID),]
phenotypes <- phenotypes[match(idorder, phenotypes$SAMPID),]
prs <- prs[match(idorder, prs$SAMPID),]

t_expr[1:5,1:5];dim(t_expr)
phenotypes[1:5,1:5];dim(phenotypes)
como[1:5,1:5]; dim(como)
prs[1:5,1:5];dim(prs)

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
# # identically ordered?
identical(t_expr$SAMPID,phenotypes$SAMPID)
identical(t_expr$SAMPID,prs$SAMPID)

# SAVE EXPRESSION, PRS AND PHENOTYPES DATA

print("Saving expr data")
fwrite(t_expr, file=paste0(resultsdir, "expression", "_" , which_subset, ".csv"))


# save phenotype data
print("Saving phenotypes data")
fwrite(phenotypes, file=paste0(resultsdir, "phenotypes", "_" , which_subset,".csv"))

# save comorbidities data
print("Saving phenotypes data")
fwrite(como, file=paste0(resultsdir, "como", "_" , which_subset,".csv"))

# save PRS 
print("Saving PRS data")
fwrite(prs, file=paste0(resultsdir, "prs", "_" , which_subset, ".csv"))

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################

# MAPPING EXPR-BIO
print("Preparing mapping expr pheno")

prior_expr_pheno <- merge(colnames(phenotypes)[!(colnames(phenotypes) %in% c("SUBJID", "subid", "SAMPID", "SMTS", "SMTSD"))],
  colnames(t_expr)[!(colnames(t_expr) %in% c("SUBJID","subid", "SAMPID", "SMTS", "SMTSD", "DTHHRDY"))]) %>% 
  dplyr::rename(bio=x, gene=y) %>% setDT %>% .[, .(gene, bio)]

fwrite(prior_expr_pheno, file=paste0(resultsdir, "prior_expr_pheno_", which_subset,".csv"))

##############################################################################################
##############################################################################################

# MAPPING EXPR-PRS
print("Preparing mapping expr prs")

prior_expr_prs <- merge(colnames(prs)[!(colnames(prs) %in% c("SAMPID"))],
  colnames(t_expr)[!(colnames(t_expr) %in% c("SUBJID","subid", "SAMPID", "SMTS", "SMTSD", "DTHHRDY"))]) %>% 
  dplyr::rename(prs=x, gene=y) %>% setDT %>% .[, .(gene, prs)]

fwrite(prior_expr_prs, file=paste0(resultsdir, "prior_expr_prs_", which_subset,".csv"))

##############################################################################################
##############################################################################################


# MAPPING EXPR-BIO
print("Preparing mapping expr como")

prior_expr_como <- merge(colnames(como)[!(colnames(como) %in% c("SUBJID","SAMPID"))],
  colnames(t_expr)[!(colnames(t_expr) %in% c("SUBJID","subid", "SAMPID", "SMTS", "SMTSD", "DTHHRDY"))]) %>% 
  dplyr::rename(como=x, gene=y) %>% setDT %>% .[, .(gene, como)]

fwrite(prior_expr_como, file=paste0(resultsdir, "prior_expr_como_", which_subset,".csv"))

##############################################################################################
##############################################################################################

print("All done")
