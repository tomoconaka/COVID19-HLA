setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(stringr)
couples <- readRDS(file="couples_step4_incl_samesex_whichfirst.rds")

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

saveRDS(couples, file="couples_step4_WB.rds")

write.table(cbind(couples$f.eid, couples$f.eid), file="WB.couples", sep="\t", quote=F, col.names = F, row.names = F)

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

couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_WB_withHLA.rds")

##SAS
couples <- readRDS(file="couples_step4_incl_samesex_whichfirst.rds")
pc <- fread("/home/richards/tomoko.nakanishi/my_project/14.PanUKBB/data/SAS.pc")
pc <- pc %>% dplyr::select(FID)
pc <- pc %>% inner_join(p, by=c("FID"="FID")) %>% dplyr::select(FID, PC1:PC10)

couples <- couples %>% inner_join(pc, by=c("f.eid" = "FID"))
couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_SAS_withHLA.rds")

# ##EAS
# couples <- readRDS(file="couples_step4_incl_samesex_whichfirst.rds")
# pc <- fread("/home/richards/tomoko.nakanishi/my_project/14.PanUKBB/data/EAS.pc")
# 
# couples <- couples %>% inner_join(pc, by=c("f.eid" = "FID"))
# couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
# couples <- couples %>% group_by(uniqueCode) %>%
#   filter(length(f.eid) == 2) %>% ungroup()
# 
# #saveRDS(couples, file="couples_step4_EAS_withHLA.rds")
# 
# ##AFR
# couples <- readRDS(file="couples_step4_incl_samesex_whichfirst.rds")
# pc <- fread("/home/richards/tomoko.nakanishi/my_project/14.PanUKBB/data/AFR.pc")
# 
# couples <- couples %>% inner_join(pc, by=c("f.eid" = "FID"))
# couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
# couples <- couples %>% group_by(uniqueCode) %>%
#   filter(length(f.eid) == 2) %>% ungroup()
