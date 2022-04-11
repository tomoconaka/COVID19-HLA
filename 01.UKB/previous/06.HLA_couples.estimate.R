setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(stringr)
couples <- readRDS(file="couples_step1.rds")
##HLA
batch10 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_10.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch11 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_11.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch12 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_12.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch13 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_13.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch14 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_14.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch15 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_15.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch16 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_16.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch17 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_17.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch18 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_18.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch19 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_19.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch20 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_20.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch21 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_21.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch22 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_22.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch23 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_23.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch24 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_24.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch25 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_25.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch26 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_26.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch27 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_27.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch28 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_28.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch29 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_29.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch30 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_30.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch31 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_31.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch32 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_32.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch33 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_33.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch34 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_34.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch35 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_35.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch36 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_36.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch37 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_37.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch38 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_38.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch39 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_39.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch40 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_40.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch41 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_41.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch42 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_42.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch43 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_43.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch44 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_44.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch45 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_45.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch46 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_46.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch47 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_47.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch48 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_48.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch49 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_49.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch50 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_50.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch51 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_51.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch52 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_52.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch53 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_53.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch54 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_54.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch55 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_55.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch56 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_56.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch57 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_57.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch58 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_58.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch59 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_59.tsv.gz") %>% filter(ID %in% couples$f.eid)
batch60 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_60.tsv.gz") %>% filter(ID %in% couples$f.eid)

hla_all <- bind_rows(batch10, batch11, batch12, batch13, batch14, batch15, batch16, batch17, batch18, batch19, 
                     batch20, batch21, batch22, batch23, batch24, batch25, batch26, batch27, batch28, batch29,
                     batch30, batch31, batch32, batch33, batch34, batch35, batch36, batch37, batch38, batch39,
                     batch40, batch41, batch42, batch43, batch44, batch45, batch46, batch47, batch48, batch49,
                     batch50, batch51, batch52, batch53, batch54, batch55, batch56, batch57, batch58, batch59, batch60)

couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

## two digits
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[5:66]), .funs=funs(paste0(unlist(strsplit(.,":"))[1])))

couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQA1 = list(c(DQA1_1, DQA1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% group_by(uniqueCode) %>% 
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



# tmp <- couples %>% filter(uniqueCode == "1_278000_665000_11004_1_1_21_2_3") %>% dplyr::select(HLACmatch)

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(PC_distance = sqrt(abs(PC1[1] - PC1[2])^2 + abs(PC2[1] - PC2[2])^2 +
                              abs(PC3[1] - PC3[2])^2 + abs(PC3[1] - PC3[2])^2 +
                              abs(PC4[1] - PC4[2])^2 + abs(PC4[1] - PC4[2])^2 + 
                              abs(PC5[1] - PC5[2])^2 + abs(PC5[1] - PC5[2])^2 + 
                              abs(PC6[1] - PC6[2])^2 + abs(PC7[1] - PC7[2])^2 + 
                              abs(PC8[1] - PC8[2])^2 + abs(PC8[1] - PC8[2])^2 + 
                              abs(PC9[1] - PC9[2])^2 + abs(PC9[1] - PC9[2])^2),
         age_diff = abs(f.21022.0.0[1]- f.21022.0.0[2])) %>% ungroup()


# couples <- couples %>% group_by(uniqueCode) %>% 
#   mutate(HLAAmatch = length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))), 
#          HLABmatch = length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))),
#          HLACmatch = length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))),
#          HLADRB1match = length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))),
#          HLADQA1match = length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))),
#          HLADQB1match = length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))))
# 

couples <- couples %>% 
  mutate(HLAclass1match = case_when(HLAAmatch == 1 & HLABmatch == 1 & HLACmatch == 1 ~ 1,
                                    TRUE ~ 0),
         HLAclass2match = case_when(HLADRB1match == 1 & HLADQA1match == 1 & HLADQB1match == 1 ~ 1,
                                    TRUE ~ 0))

tmp <- couples %>% dplyr::select(uniqueCode, PC_distance, age_diff, colnames(couples)[grepl("match", colnames(couples))]) %>% unique()

library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
Vars <- c(catVars, "PC_distance", "age_diff")
tableOne <- CreateTableOne(vars = Vars, data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
write.csv(tab3Mat, file = "Table_baseline_couples_onefield.csv")

couples <- readRDS(file="couples_step1.rds")
couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

## four digits
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[5:66]), .funs=funs(paste0(unlist(strsplit(.,":"))[1], ":", unlist(strsplit(.,":"))[2])))


couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQA1 = list(c(DQA1_1, DQA1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% group_by(uniqueCode) %>% 
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

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(PC_distance = sqrt(abs(PC1[1] - PC1[2])^2 + abs(PC2[1] - PC2[2])^2 +
                              abs(PC3[1] - PC3[2])^2 + abs(PC3[1] - PC3[2])^2 +
                              abs(PC4[1] - PC4[2])^2 + abs(PC4[1] - PC4[2])^2 + 
                              abs(PC5[1] - PC5[2])^2 + abs(PC5[1] - PC5[2])^2 + 
                              abs(PC6[1] - PC6[2])^2 + abs(PC7[1] - PC7[2])^2 + 
                              abs(PC8[1] - PC8[2])^2 + abs(PC8[1] - PC8[2])^2 + 
                              abs(PC9[1] - PC9[2])^2 + abs(PC9[1] - PC9[2])^2),
         age_diff = abs(f.21022.0.0[1]- f.21022.0.0[2])) %>% ungroup()


# couples <- couples %>% group_by(uniqueCode) %>% 
#   mutate(HLAAmatch = length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))), 
#          HLABmatch = length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))),
#          HLACmatch = length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))),
#          HLADRB1match = length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))),
#          HLADQA1match = length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))),
#          HLADQB1match = length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))))
# 

couples <- couples %>% 
  mutate(HLAclass1match = case_when(HLAAmatch == 1 & HLABmatch == 1 & HLACmatch == 1 ~ 1,
                                    TRUE ~ 0),
         HLAclass2match = case_when(HLADRB1match == 1 & HLADQA1match == 1 & HLADQB1match == 1 ~ 1,
                                    TRUE ~ 0))

tmp <- couples %>% dplyr::select(uniqueCode, PC_distance, age_diff, colnames(couples)[grepl("match", colnames(couples))]) %>% unique()

library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
Vars <- c(catVars, "PC_distance", "age_diff")
tableOne <- CreateTableOne(vars = Vars, data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
write.csv(tab3Mat, file = "Table_baseline_couples_twofield.csv")


