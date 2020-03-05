## ---------------------------
##
## Script name: Combine_BED.R
##
## Purpose of script: Combine the bed files used for target capture and assing a uniq ID for each region
##
## Author: Chathura J. Gunasekara, PhD
##
## Date Created: 2020-02-26
##
## Copyright (c) Chathura J Gunasekara, 2020
## Email: gunaseka@bcm.edu/cjgunase@mtu.edu
##
## ---------------------------
##
## Notes: Every tine the script runs, it will assin a different uniq ID so, do not overwrite an exiting file.
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("~/Documents/")          # set to paren directory where the bed files resides


## ---------------------------

options(scipen = 6, digits = 4) # I prefer to view outputs in non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

require(tidyverse)
require(data.table)
require(bedr)
library(tidyr)
library(ids)
## ---------------------------


if (check.binary("bedtools")) {

    index <- list()
    for (file in dir("./beds",full.names = T)){
        print(file)
        file_base = basename(file)
        temp_bed <- read.table(file)
        bed_string <- paste0(paste0(temp_bed$V1,":",temp_bed$V2),"-",temp_bed$V3)
        is.valid.region(bed_string)
        
        if(!is.sorted.region(bed_string)){
            bed_string <- bedr.sort.region(bed_string)
            index[[file_base]] <- bed_string
        }else{
            index[[file_base]] <- bed_string
        }
        
        
    }
    combine_all_bed <- bedr.join.multiple.region(index)
    
    
}

index <- bedr.merge.region(combine_all_bed$index)
uniq_id <- ids::random_id(n=length(index),bytes = 4)
combine_all_bed <- data.frame(index,uniq_id)


temp <- separate(data = combine_all_bed,col = "index",into = c("chr","start-end"),sep = ":")
final <- separate(data = temp,col = "start-end",into = c("start","end"),sep = "-")
final$UCSC_Coord <- combine_all_bed$index
write.table(final,file = "./combined_corsivs.bed",col.names = F,row.names = F,quote = F,sep = "\t")


