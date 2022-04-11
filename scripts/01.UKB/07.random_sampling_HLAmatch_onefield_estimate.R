setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(stringr)
library(data.table)
library(dplyr)
library(tidyr)

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

##HLA
batch10 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_10.tsv.gz") 
batch11 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_11.tsv.gz") 
batch12 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_12.tsv.gz") 
batch13 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_13.tsv.gz") 
batch14 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_14.tsv.gz") 
batch15 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_15.tsv.gz") 
batch16 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_16.tsv.gz") 
batch17 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_17.tsv.gz") 
batch18 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_18.tsv.gz") 
batch19 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_19.tsv.gz") 
batch20 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_20.tsv.gz") 
batch21 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_21.tsv.gz") 
batch22 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_22.tsv.gz") 
batch23 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_23.tsv.gz") 
batch24 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_24.tsv.gz") 
batch25 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_25.tsv.gz") 
batch26 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_26.tsv.gz") 
batch27 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_27.tsv.gz") 
batch28 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_28.tsv.gz") 
batch29 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_29.tsv.gz") 
batch30 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_30.tsv.gz") 
batch31 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_31.tsv.gz") 
batch32 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_32.tsv.gz") 
batch33 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_33.tsv.gz") 
batch34 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_34.tsv.gz") 
batch35 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_35.tsv.gz") 
batch36 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_36.tsv.gz") 
batch37 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_37.tsv.gz") 
batch38 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_38.tsv.gz") 
batch39 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_39.tsv.gz") 
batch40 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_40.tsv.gz") 
batch41 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_41.tsv.gz") 
batch42 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_42.tsv.gz") 
batch43 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_43.tsv.gz") 
batch44 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_44.tsv.gz") 
batch45 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_45.tsv.gz") 
batch46 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_46.tsv.gz") 
batch47 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_47.tsv.gz") 
batch48 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_48.tsv.gz") 
batch49 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_49.tsv.gz") 
batch50 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_50.tsv.gz") 
batch51 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_51.tsv.gz") 
batch52 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_52.tsv.gz") 
batch53 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_53.tsv.gz") 
batch54 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_54.tsv.gz") 
batch55 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_55.tsv.gz") 
batch56 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_56.tsv.gz") 
batch57 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_57.tsv.gz") 
batch58 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_58.tsv.gz") 
batch59 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_59.tsv.gz") 
batch60 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_60.tsv.gz") 

hla_all <- bind_rows(batch10, batch11, batch12, batch13, batch14, batch15, batch16, batch17, batch18, batch19, 
                     batch20, batch21, batch22, batch23, batch24, batch25, batch26, batch27, batch28, batch29,
                     batch30, batch31, batch32, batch33, batch34, batch35, batch36, batch37, batch38, batch39,
                     batch40, batch41, batch42, batch43, batch44, batch45, batch46, batch47, batch48, batch49,
                     batch50, batch51, batch52, batch53, batch54, batch55, batch56, batch57, batch58, batch59, batch60)

hla_all <- hla_all %>% inner_join(pc, by="ID")
household <- fread("ukb27449_20688_household.tab.gz") %>% dplyr::select(f.eid, f.31.0.0, f.21022.0.0)

hla_all <- hla_all %>% inner_join(household, by=c("ID"="f.eid"))
hla_all <- hla_all %>% group_by(ID) %>% 
  mutate_at(.vars=vars(colnames(hla_all)[2:63]), .funs=funs(paste0(unlist(strsplit(.,":"))[1])))
hla_all <- hla_all %>% group_by(ID) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQA1 = list(c(DQA1_1, DQA1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

## select random 40000 pairs from UKB WB.

HLAmatch_estimate_onefield <- function(hla_all, f, m){
  set.seed(f)
  dat_female_index <- sample(hla_all$ID[hla_all$f.31.0.0 == 0], size = 40000,  replace = FALSE, prob = NULL)
  set.seed(m)
  dat_male_index <- sample(hla_all$ID[hla_all$f.31.0.0 == 1], size = 40000,  replace = FALSE, prob = NULL)
  dat_female <- hla_all %>% filter(ID %in% dat_female_index) %>% mutate(group = seq(1,40000))
  dat_male <- hla_all %>% filter(ID %in% dat_male_index) %>% mutate(group = seq(1,40000))
  dat <- bind_rows(dat_female, dat_male)
  dat <- dat %>% group_by(group) %>% 
    mutate(HLAAmatch = case_when(length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))) == 2 ~ 1,
                                 length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))) == 1 & length(unique(unlist(HLAA[1]))) == 1 & length(unique(unlist(HLAA[2]))) == 1 ~ 1,
                                 TRUE ~ 0),
           HLABmatch = case_when(length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))) == 2 ~ 1,
                                 length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))) == 1 & length(unique(unlist(HLAB[1]))) == 1 & length(unique(unlist(HLAB[2]))) == 1 ~ 1,
                                 TRUE ~ 0),
           HLACmatch = case_when(length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))) == 2 ~ 1,
                                 length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))) == 1 & length(unique(unlist(HLAC[1]))) == 1 & length(unique(unlist(HLAC[2]))) == 1 ~ 1,
                                 TRUE ~ 0),
           HLADRB1match = case_when(length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))) == 2 ~ 1,
                                    length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))) == 1 & length(unique(unlist(HLADRB1[1]))) == 1 & length(unique(unlist(HLADRB1[2]))) == 1 ~ 1,
                                    TRUE ~ 0),
           HLADQA1match = case_when(length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))) == 2 ~ 1,
                                    length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))) == 1 & length(unique(unlist(HLADQA1[1]))) == 1 & length(unique(unlist(HLADQA1[2]))) == 1 ~ 1,
                                    TRUE ~ 0),
           HLADQB1match = case_when(length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))) == 2 ~ 1,
                                    length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))) == 1 & length(unique(unlist(HLADQB1[1]))) == 1 & length(unique(unlist(HLADQB1[2]))) == 1 ~ 1,
                                    TRUE ~ 0)
    )
  dat <- dat %>% group_by(group) %>% 
    mutate(PC_distance = sqrt(abs(PC1[1] - PC1[2])^2 + abs(PC2[1] - PC2[2])^2 +
                                abs(PC3[1] - PC3[2])^2 + abs(PC3[1] - PC3[2])^2 +
                                abs(PC4[1] - PC4[2])^2 + abs(PC4[1] - PC4[2])^2 + 
                                abs(PC5[1] - PC5[2])^2 + abs(PC5[1] - PC5[2])^2 + 
                                abs(PC6[1] - PC6[2])^2 + abs(PC7[1] - PC7[2])^2 + 
                                abs(PC8[1] - PC8[2])^2 + abs(PC8[1] - PC8[2])^2 + 
                                abs(PC9[1] - PC9[2])^2 + abs(PC9[1] - PC9[2])^2),
           age_diff = abs(f.21022.0.0[1]- f.21022.0.0[2])) %>% ungroup()
  tmp <- dat %>% dplyr::select(group, PC_distance, age_diff, colnames(dat)[grepl("match", colnames(dat))]) %>% unique()
  return(c(i, median(tmp$age_diff), median(tmp$PC_distance),sum(tmp$HLAAmatch)/dim(tmp)[1], sum(tmp$HLABmatch)/dim(tmp)[1], sum(tmp$HLACmatch)/dim(tmp)[1], 
           sum(tmp$HLADRB1match)/dim(tmp)[1], sum(tmp$HLADQA1match)/dim(tmp)[1], sum(tmp$HLADQB1match)/dim(tmp)[1])) 
}

out <- data.frame(matrix(0, 100, 9))
colnames(out) <- c("iter", "median_age_diff","median_PC_distance","HLAAmatch", "HLABmatch", "HLACmatch", "HLADRB1match",
                   "HLADQA1match", "HLADQB1match")
for(i in seq(1, 100)){
  print(i)
  out[i,] <- HLAmatch_estimate_onefield(hla_all, i, (100-i))
}

saveRDS(out, file="bootstrap_randomsample_HLAmatch_onefield.rds")

couplesmatch <- fread("Table_baseline_couples_onefield.csv")
