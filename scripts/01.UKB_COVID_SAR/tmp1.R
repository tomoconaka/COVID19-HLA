setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")
library(tidyr)
library(dplyr)
## two digits
couples <- readRDS("couples_step4_WB_withHLA.rds")
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

couples <- couples %>% mutate(outcome = case_when(coinfection == TRUE & interval >= 1 & interval <= 7 ~ 1,
                                                  TRUE ~ 0))
couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(second_sex = f.31.0.0[2])

tmp <- couples %>% dplyr::select(uniqueCode, age_diff, age_sum, interval, f.31.0.0, second_sex, colnames(couples)[grepl("match", colnames(couples))], coinfection,
                                 colnames(couples)[grepl("dif$", colnames(couples))]) %>% unique()
tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))
tmp <- tmp[!duplicated(tmp$uniqueCode),]

out <- data.frame(matrix(0, 6, 4))
colnames(out) <- c("HLAgene", "beta", "se", "pval")

HLAgene <- colnames(couples)[grepl("match", colnames(couples))]

tmp <- tmp %>% mutate(interval = ifelse(interval > 7, 7, interval))

library(survival)
hla.cox <- coxph(Surv(interval, outcome) ~ HLAAmatch + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif, data=tmp)
out$HLAgene <- HLAgene
out$beta[1] <- summary(hla.cox)$coefficients[1, 1]
out$se[1] <- summary(hla.cox)$coefficients[1, 3]
out$pval[1] <- summary(hla.cox)$coefficients[1, 5]

hla.cox <- coxph(Surv(interval, outcome) ~ HLABmatch + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif, data=tmp)
out$beta[2] <- summary(hla.cox)$coefficients[1, 1]
out$se[2] <- summary(hla.cox)$coefficients[1, 3]
out$pval[2] <- summary(hla.cox)$coefficients[1, 5]

hla.cox <- coxph(Surv(interval, outcome) ~ HLACmatch + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif, data=tmp)
out$beta[3] <- summary(hla.cox)$coefficients[1, 1]
out$se[3] <- summary(hla.cox)$coefficients[1, 3]
out$pval[3] <- summary(hla.cox)$coefficients[1, 5]

hla.cox <- coxph(Surv(interval, outcome) ~ HLADRB1match + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif, data=tmp)
out$beta[4] <- summary(hla.cox)$coefficients[1, 1]
out$se[4] <- summary(hla.cox)$coefficients[1, 3]
out$pval[4] <- summary(hla.cox)$coefficients[1, 5]

hla.cox <- coxph(Surv(interval, outcome) ~ HLADQA1match + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif, data=tmp)
out$beta[5] <- summary(hla.cox)$coefficients[1, 1]
out$se[5] <- summary(hla.cox)$coefficients[1, 3]
out$pval[5] <- summary(hla.cox)$coefficients[1, 5]

hla.cox <- coxph(Surv(interval, outcome) ~ HLADQB1match + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif, data=tmp)
out$beta[6] <- summary(hla.cox)$coefficients[1, 1]
out$se[6] <- summary(hla.cox)$coefficients[1, 3]
out$pval[6] <- summary(hla.cox)$coefficients[1, 5]

write.table(out, file="onefield_cox_WB.tsv", col.names = T, row.names = F, sep="\t", quote = F)
tmp_cor <- tmp %>% ungroup() %>% dplyr::select(colnames(couples)[grepl("match", colnames(couples))])
library(corrplot)
library(Hmisc)
M <- Hmisc::rcorr(as.matrix(tmp_cor))
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

library(metafor)

U <- M$r
vi <- out$se
yi <- out$beta

V <- U*matrix(vi)%*%t(matrix(vi))

metaresult <- summary(rma.mv(yi,V))

res <- rma.mv(yi,V)
png("onefield.cox.meta_WB.png", width=500, height = 400)
forest(res, atransf=exp, slab = out$HLAgene)
dev.off()
save(tab3Mat, out, M, metaresult, res, file = "onefield_WB.Rdata")



