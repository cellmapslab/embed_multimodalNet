node_name <- fread("../../../data/lookup/node_name_mapping.txt")

change_names <- function(dt, lookup, mycol="name"){
  # function to replace names of nodes, using a lookup table
  ## dt: data frame to change
  ## mycol: column in dt to change
  ## lookup table with 'old' and 'new' column, representing old and new names
  col_order <- as.character(colnames(dt))
  dt[, myid:=1:nrow(dt)]
  dt <- merge(dt, lookup[,.(old, old2)], by.x=mycol, by.y="old", all.x=T) %>% dplyr::rename(new=old2)
  dt <- dt[!is.na(new), as.character(mycol):=new][,new:=NULL]
  dt <- dt[order(myid)]
  dt <- dt[,..col_order]
  dt[,myid:=NULL]
  
  return(dt)
}
order_type <- function(dt){
  preferred.order <- c("pheno","como", "prs", "tissue", "cova", "genes")
  dt[, xFac := factor(type, levels=preferred.order)]
  dt <- setorder(dt, xFac)
  dt[,xFac:=NULL]
  return(dt)
}



