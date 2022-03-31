#!/usr/bin/Rscript

library(XML)
library(purrr)
library(data.table)

doc = xmlParse("/scratch/data/77854/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/PhenotypeFiles/phs000424.v8.pht002742.v8.GTEx_Subject_Phenotypes.data_dict.xml")


mylist=xmlToList(doc)
mylist=mylist[-1]
mylist=mylist[-190]
library(data.table)
dt =rbindlist(mylist,fill=T, idcol=NULL)
colnames(dt) <- c("name", "description", "type", "comment", "attrs", "value1", "value2", "value3", "unit", "value4", "value5", "value6")



add_leg <- function(x){
dt[x,val1:=paste(dt[x]$value1,dt[(x+1)]$value1, sep="_")] 
dt[x,val2:=paste(dt[x]$value2,dt[(x+1)]$value2, sep="_")] 
dt[x,val3:=paste(dt[x]$value3,dt[(x+1)]$value3, sep="_")] 
dt[x,val4:=paste(dt[x]$value4,dt[(x+1)]$value4, sep="_")] 
dt[x,val5:=paste(dt[x]$value5,dt[(x+1)]$value5, sep="_")] 
dt[x,val6:=paste(dt[x]$value6,dt[(x+1)]$value6, sep="_")] 
return(dt)
}

myvec = (which(duplicated(dt$description))) -1
map(myvec, add_leg)

dt_order = dt[!which(duplicated(dt$description))]
dt_order = dt_order[,.(name,description,attrs,val1, val2, val3, val4, val5,val6)]
dt_order
fwrite(dt_order, file="../data/Extended_phenotypes.csv")
