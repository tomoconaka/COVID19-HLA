setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(stringr)
couples <- readRDS(file="couples_step4_whichfirst.rds")

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
death_data <- fread("/scratch/richards/tomoko.nakanishi/DATA/UKB/death_20210930.txt.gz")
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_WB.rds")

write.table(cbind(couples$f.eid, couples$f.eid), file="WB.couples", sep="\t", quote=F, col.names = F, row.names = F)

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
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_WB_withHLA.rds")


## two digits
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[34:95]), .funs=funs(paste0(unlist(strsplit(.,":"))[1])))

couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQA1 = list(c(DQA1_1, DQA1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% arrange(desc(first))

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(HLAAmatch = case_when(all(unlist(HLAA[1]) %in% unlist(HLAA[2])) ~ 1,
                               TRUE ~ 0),
         HLABmatch = case_when(all(unlist(HLAB[1]) %in% unlist(HLAB[2])) ~ 1,
                               TRUE ~ 0),
         HLACmatch = case_when(all(unlist(HLAC[1]) %in% unlist(HLAC[2])) ~ 1,
                               TRUE ~ 0),
         HLADRB1match = case_when(all(unlist(HLADRB1[1]) %in% unlist(HLADRB1[2])) ~ 1,
                                  TRUE ~ 0),
         HLADQA1match = case_when(all(unlist(HLADQA1[1]) %in% unlist(HLADQA1[2])) ~ 1,
                                  TRUE ~ 0),
         HLADQB1match = case_when(all(unlist(HLADQB1[1]) %in% unlist(HLADQB1[2])) ~ 1,
                                  TRUE ~ 0)
  )

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(PC1dif = abs(PC1[1] - PC1[2]),
         PC2dif = abs(PC2[1] - PC2[2]),
         PC3dif = abs(PC3[1] - PC3[2]),
         PC4dif = abs(PC4[1] - PC4[2]),
         PC5dif = abs(PC5[1] - PC5[2]),
         PC6dif = abs(PC6[1] - PC6[2]),
         PC7dif = abs(PC7[1] - PC7[2]),
         PC8dif = abs(PC8[1] - PC8[2]),
         PC9dif = abs(PC9[1] - PC9[2]),
         PC10dif = abs(PC10[1] - PC10[2]),
         age_diff = abs(f.21022.0.0[1]- f.21022.0.0[2]),
         age_sum = f.21022.0.0[1] + f.21022.0.0[2]) %>% ungroup()

couples <- couples %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                                  TRUE ~ 0))

tmp <- couples %>% dplyr::select(uniqueCode, age_diff, age_sum, f.31.0.0, colnames(couples)[grepl("match", colnames(couples))], coinfection,
                                 colnames(couples)[grepl("dif$", colnames(couples))]) %>% unique()
tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))

library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
tableOne <- CreateTableOne(vars = catVars, strata = c("outcome"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
write.csv(tab3Mat, file = "Table_onefield.csv")


out <- data.frame(matrix(0, 6, 4))
colnames(out) <- c("HLAgene", "beta", "se", "pval")

HLAgene <- colnames(couples)[grepl("match", colnames(couples))]

for(i in seq(1, 6)){
  LM <- glm(paste0("outcome ~ ",HLAgene[i]," + age_diff + age_sum + f.31.0.0 + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif"), data=tmp, family="binomial")
  out$HLAgene[i] <- HLAgene[i]
  out$beta[i] <- summary(LM)$coefficients[2, 1]
  out$se[i] <- summary(LM)$coefficients[2, 2]
  out$pval[i] <- summary(LM)$coefficients[2, 4]
}

write.table(out, file="onefield_regression.tsv", col.names = T, row.names = F, sep="\t", quote = F)

tmp_cor <- tmp %>% dplyr::select(colnames(couples)[grepl("match", colnames(couples))])
library(corrplot)
M <- rcorr(as.matrix(tmp_cor))
library(Hmisc)
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

saveRDS(M, file="HLA_onefield_correlation.rds")


library(metafor)

U <- M$r
vi <- out$se
yi <- out$beta

V <- U*matrix(vi)%*%t(matrix(vi))

metaresult <- summary(rma.mv(yi,V))

save(tab3Mat, out, M, metaresult, file = "onefield.Rdata")
