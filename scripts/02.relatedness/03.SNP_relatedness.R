setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/")
library(data.table)
library(dplyr)
library(tidyr)
library(R.utils)

args <- commandArgs(trailingOnly = T)

data <- fread(paste0(args[1],"/WB.couples_",args[2],".raw"))
if(args[1] == "randompairs"){
  couples <- fread(paste0(args[1],"/WBcouples.sample"))
  couples <- couples %>% arrange(V3)
  colnames(couples) <- c("FID", "IID", "uniqueCode")
  couples <- couples %>% inner_join(data, by=c("IID"="IID"))
}
if(args[1] == "realcouples"){
  couples <- readRDS("couples_step1.rds")
  couples <- couples %>% arrange(uniqueCode)
  couples <- couples %>% inner_join(data, by=c("f.eid"="IID"))
}

assign_value <- function(list_object){
  x <- as.numeric(list_object[1])
  y <- as.numeric(list_object[2])
  if(is.na(x) | is.na(y)){
    return(NA)
  }
  if(x == 0 & y == 0){
    return(1)
  }
  if(x == 2 & y == 2){
    return(1)
  }
  if(x == 0 & y == 2){
    return(0)
  }
  if(x == 2 & y == 0){
    return(0)
  }
  else{
    return(0.5)
    }
}

couples <- data.frame(couples)

SNP_diversity <- couples %>% group_by(uniqueCode) %>%
  mutate_at(.vars=vars(colnames(couples)[grepl("rs|Aff", colnames(couples))]),
            .funs=funs(assign_value(.))) %>%
  distinct(uniqueCode, .keep_all=TRUE)

SNP_diversity <- data.frame(SNP_diversity)
SNP_diversity <- SNP_diversity %>% 
  mutate(Qc = rowMeans(.[, c(colnames(couples)[grepl("rs|Aff", colnames(couples))])], na.rm=T),
         Nsnp = length(colnames(couples)[grepl("rs|Aff", colnames(couples))])) %>% select(uniqueCode, Qc, Nsnp)

saveRDS(SNP_diversity, file=paste0(args[1],"/SNP_diversity_WBcouples_",args[2],".rds"))

#data <- data %>% select(!colnames(data)[grepl("rs", colnames(data))])

#Affx-
