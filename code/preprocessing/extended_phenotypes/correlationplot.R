#!/usr/bin/Rscript
library(data.table)
library(magrittr)
library(ggplot2)

source('../../../code/general/add_fun.R')
colors <- fread( file="../../../data/lookup/color_code.txt")
colors[node_types=="prs", node_colors:="#355C7D"]
colors[node_types=="pheno", node_colors:="#F8B195"]

node_anno <-  fread( file="../../../data/lookup/node_anno.txt")
node_color <- merge(colors, node_anno, by.x="node_types", by.y="type")
#####################################
sam_att <- fread("../../../data/kimono/pre/sample_attributes_biogrid.csv")
sam_att <- sam_att[,.(SUBJID, SAMPID)]
#####################################
como <- fread(file="../../../data/kimono/pre/como_biogrid.csv")
como <- unique(como[,SAMPID:=NULL])

#####################################
pheno <- fread("../../../data/kimono/pre/phenotypes_biogrid.csv")
pheno <-pheno[,.(SUBJID, SEX, AGE, BMI)]
pheno <- unique(pheno)

#####################################
prs <- fread("../../../data/kimono/pre/prs_biogrid.csv")
prs <- merge(prs, sam_att, by="SAMPID")
prs <- unique(prs[, SAMPID:=NULL])
#####################################
# merge data frames
dt1 <- merge(pheno, prs, by="SUBJID")
dt2 <- merge(dt1, como, by="SUBJID")
dt3 <- setnames(dt2, node_name$old, node_name$new,skip_absent=TRUE)
dt3[, SUBJID:=NULL]
plt_bar <- function(x){
	plt <-ggplot(plt_data, aes_string(as.character(x), fill="type")) + geom_bar() + geom_bar(position = "dodge2") + theme_minimal()
	return(plt)
}

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }



reorder_cormat <- function(cormat){
	dd <- as.dist((1-cormat)/2)
	hc <- hclust(dd)
	cormat <-cormat[hc$order, hc$order]
	return(cormat)
}



cormat <- round(cor(dt3),2)
# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- reshape2::melt(upper_tri, na.rm = TRUE)

node_color <- change_names(node_color, node_name, mycol="name")

melted_cormat1 <- merge(melted_cormat, node_color[,.(node_colors, name)],
  by.x="Var1", by.y="name", all.x=T) %>% dplyr::rename(color_Var1=node_colors)
melted_cormat2 <- merge(melted_cormat1, node_color[,.(node_colors, name)],
 by.x="Var2", by.y="name", all.x=T) %>% dplyr::rename(color_Var2=node_colors)

color_lookup <- unique(as.data.table(melted_cormat2)[,.(Var1, color_Var1)] ) %>% dplyr::rename(node=Var1, color=color_Var1)
color_lookup <- color_lookup[match( as.character(unique(melted_cormat$Var2)),color_lookup$node)]

mycolors <- setNames(color_lookup$color, color_lookup$node)


# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(x=Var2, y=Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#C7CEEA", high = "#F67280", mid = "#FFFFFF", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1, color=mycolors),
    axis.text.y=element_text(color=mycolors))+
 coord_fixed()  + 
geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.justification = c(1, 0),
  legend.position = c(0.6, 0.7),
  legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                title.position = "top", title.hjust = 0.5))



pdf(height=10, width=10, 
  file=paste0("results/data_cor/", Sys.Date(), "_cor_all.pdf"))
ggheatmap
dev.off()

