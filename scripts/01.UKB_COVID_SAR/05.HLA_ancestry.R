setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(stringr)
couples <- readRDS(file="couples_step4_incl_samesex_whichfirst_20230414.rds")

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))
pc <- pc %>% dplyr::select(ID)

p <- read.table("/project/richards/tomoko.nakanishi/14.PanUKBB/data/UKBfullcPCrev.list")
gp <- fread("/project/richards/tomoko.nakanishi/14.PanUKBB/data/EMC-1KGUKB-GENO-PC125-C12-HM3.class")
colnames(p) <- colnames(gp)
pc <- pc %>% inner_join(p, by=c("ID"="FID")) %>% dplyr::select(ID, PC1:PC10)

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_WB_20230414.rds")

write.table(cbind(couples$f.eid, couples$f.eid), file="WB.couples", sep="\t", quote=F, col.names = F, row.names = F)

batch10 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_10.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch11 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_11.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch12 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_12.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch13 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_13.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch14 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_14.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch15 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_15.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch16 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_16.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch17 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_17.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch18 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_18.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch19 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_19.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch20 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_20.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch21 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_21.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch22 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_22.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch23 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_23.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch24 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_24.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch25 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_25.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch26 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_26.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch27 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_27.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch28 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_28.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch29 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_29.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch30 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_30.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch31 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_31.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch32 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_32.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch33 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_33.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch34 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_34.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch35 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_35.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch36 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_36.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch37 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_37.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch38 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_38.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch39 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_39.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch40 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_40.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch41 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_41.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch42 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_42.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch43 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_43.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch44 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_44.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch45 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_45.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch46 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_46.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch47 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_47.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch48 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_48.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch49 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_49.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch50 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_50.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch51 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_51.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch52 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_52.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch53 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_53.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch54 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_54.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch55 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_55.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch56 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_56.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch57 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_57.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch58 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_58.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch59 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_59.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))
batch60 <- fread("/project/richards/guillaume.butler-laporte/HLA/ukb_wes/calls/hla_df_batch_qced_60.tsv.gz") %>% mutate_all( .funs=funs(as.character(.)))


hla_all <- bind_rows(batch10, batch11, batch12, batch13, batch14, batch15, batch16, batch17, batch18, batch19,
                     batch20, batch21, batch22, batch23, batch24, batch25, batch26, batch27, batch28, batch29,
                     batch30, batch31, batch32, batch33, batch34, batch35, batch36, batch37, batch38, batch39,
                     batch40, batch41, batch42, batch43, batch44, batch45, batch46, batch47, batch48, batch49,
                     batch50, batch51, batch52, batch53, batch54, batch55, batch56, batch57, batch58, batch59, batch60)

couples$f.eid <- as.character(couples$f.eid)
couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_WB_withHLA_20230414.rds")

##SAS
couples <- readRDS(file="couples_step4_incl_samesex_whichfirst_20230414.rds")
pc <- fread("/home/richards/tomoko.nakanishi/my_project/14.PanUKBB/data/SAS.pc")
pc <- pc %>% dplyr::select(FID)
pc <- pc %>% inner_join(p, by=c("FID"="FID")) %>% dplyr::select(FID, PC1:PC10)

couples <- couples %>% inner_join(pc, by=c("f.eid" = "FID"))
couples$f.eid <- as.character(couples$f.eid)
couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_SAS_withHLA_20230414.rds")

##EAS
# couples <- readRDS(file="couples_step4_incl_samesex_whichfirst_20230414.rds")
# pc <- fread("/home/richards/tomoko.nakanishi/my_project/14.PanUKBB/data/EAS.pc")
# 
# couples <- couples %>% inner_join(pc, by=c("f.eid" = "FID"))
# couples$f.eid <- as.character(couples$f.eid)
# couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
# couples <- couples %>% group_by(uniqueCode) %>%
#   filter(length(f.eid) == 2) %>% ungroup()
# dim(couples)
# 
# #saveRDS(couples, file="couples_step4_EAS_withHLA.rds")
# 
##AFR
# couples <- readRDS(file="couples_step4_incl_samesex_whichfirst_20230414.rds")
# pc <- fread("/home/richards/tomoko.nakanishi/my_project/14.PanUKBB/data/AFR.pc")
# 
# couples <- couples %>% inner_join(pc, by=c("f.eid" = "FID"))
# couples$f.eid <- as.character(couples$f.eid)
# couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
# couples <- couples %>% group_by(uniqueCode) %>%
#   filter(length(f.eid) == 2) %>% ungroup()
# dim(couples)
