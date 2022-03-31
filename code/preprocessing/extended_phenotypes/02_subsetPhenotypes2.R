#!/usr/bin/Rscript
library(data.table)
library(magrittr)
library(expss)
pheno_raw <- fread("scratch/data/77854/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/PhenotypeFiles/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt", skip=10)

legend <- fread("../../../data/phenotypes/Extended_phenotypes.csv")

print("Select phenotypes")




name_subset <- c("dbGaP_Subject_ID","SUBJID", "AGE","BMI", "COHORT", "SEX", "RACE", "MHRNLFLR", "MHCOPD", "MHCVD","MHHRTDIS", "MHHTN", "MHLVRDIS",  "MHT1D", "MHT2D", "MHPNMNIA","MHPNMIAB", "MHASTHMA", "MHDPRSSN")



pheno <- pheno_raw[,..name_subset]
# recode SEX
pheno[SEX==1, SEX:=100] # 1 = males; make them the reference as there are more samples
pheno[SEX==2, SEX:=200]  
pheno[SEX==100, SEX:=0]
pheno[SEX==200, SEX:=1]

# recode RACE
pheno[RACE==3, RACE:=0] # white
pheno[RACE==2, RACE:=200] # african american
pheno[RACE==1, RACE:=100] # asian
pheno[RACE==4, RACE:=400] # american indian

pheno[RACE==200, RACE:=1]
pheno[RACE==100, RACE:=2]
pheno[RACE==400, RACE:=3]

# recode COHORT
pheno[COHORT=="Postmortem", COHORT:=0]
pheno[COHORT=="Organ Donor (OPO)", COHORT:=1]
pheno[COHORT=="Surgical", COHORT:=2]
pheno$COHORT <- as.integer(pheno$COHORT)

pheno <- na_if(pheno, 98, with_labels = FALSE)
pheno <- na_if(pheno, 99, with_labels = FALSE)
pheno <- na.omit(pheno)
colSums(is.na(pheno))
rowSums(is.na(pheno))
sapply(pheno, table)


fwrite(pheno, file="../../../data/phenotypes/subsetPhenotypes_ext.csv")


# having only the essential without the comorbidites
name_subset <- c("SUBJID","SEX", "AGE", "BMI")
pheno <- pheno[,..name_subset]
sapply(pheno, table)

fwrite(pheno, file="../../../data/phenotypes/subsetPhenotypes.csv")
