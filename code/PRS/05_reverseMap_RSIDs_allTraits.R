#Author: Ghalia Rehawi
#Date: 10.09.2021
#This script maps the SNPs IDs from format RSXXXXX to CHR_POS in the concatenated chromosome files of each trait 

#Load libraries
library(dplyr)
#library(tidyr)
library(data.table)

# Read the mapping file 
mapping_file = fread("remapping_file.nodup.bim")

base_dir = "/calculated_PRS"
files_to_read = list.files(
  path = base_dir,        # directory to search within
  pattern = "*_all_chromosomes.txt$", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)

for(f in files_to_read){
	trait = fread(f)
	joined_dt = merge(trait, mapping_file, by.x = c("V2"), by.y = c("V2"), all.x = FALSE, all.y = FALSE)
	new_dt = joined_dt[ ,.(V1.x, V3.y, V3.x, V4, V5, V6)]
	fwrite(new_dt, paste0(dirname(f), "/mapped_", basename(f)), col.names =FALSE, sep= "\t")

	
}
