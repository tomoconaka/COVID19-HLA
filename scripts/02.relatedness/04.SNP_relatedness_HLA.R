setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/")
library(data.table)
library(dplyr)
library(tidyr)
library(R.utils)

couples <- readRDS("couples_step1.rds")
couples <- couples %>% arrange(uniqueCode)

pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))
couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

#saveRDS(couples, file="relatedness/couples_WB.rds")

#write.table(cbind(couples$f.eid, couples$f.eid), file="relatedness/WB.couples", sep="\t", quote=F, col.names = F, row.names = F)

data <- fread("relatedness/WB.couples_")
couples <- couples %>% inner_join(data, by=c("f.eid"="IID"))

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
  if(x == 2 & y == 1){
    return(0)
  }
  else{
    return(0.5)
    }
}

SNP_diversity <- couples %>% group_by(uniqueCode) %>%
  mutate_at(.vars=vars(colnames(couples)[grepl("rs|Affi", colnames(couples))]),
            .funs=funs(assign_value(.))) %>%
  distinct(uniqueCode, .keep_all=TRUE)

SNP_diversity <- data.frame(SNP_diversity)
SNP_diversity <- SNP_diversity %>% 
  mutate(Qc = rowMeans(.[, c(colnames(couples)[grepl("rs|Affi", colnames(couples))])], na.rm=T)) %>% select(uniqueCode, Qc)

saveRDS(SNP_diversity, file="SNP_diversity_WBcouples_HLA.rds")


#Affx-
