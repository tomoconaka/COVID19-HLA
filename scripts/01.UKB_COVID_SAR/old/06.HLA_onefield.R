setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")
library(tidyr)
library(dplyr)
## two digits
couples <- readRDS("couples_step4_WB_withHLA_20230414.rds")
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[24:85]), .funs=funs(paste0(unlist(strsplit(.,":"))[1])))

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

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(second_sex = f.31.0.0[2])
  
tmp <- couples %>% dplyr::select(uniqueCode, age_diff, age_sum, f.31.0.0, second_sex, colnames(couples)[grepl("match", colnames(couples))], coinfection,
                                 colnames(couples)[grepl("dif$", colnames(couples))]) %>% unique()
tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))
tmp <- tmp[!duplicated(tmp$uniqueCode),]
library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
tableOne <- CreateTableOne(vars = catVars, strata = c("outcome"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
write.csv(tab3Mat, file = "Table_onefield_WB_2023044.csv")

tab3Mat <- read.csv("Table_onefield_WB_2023044.csv")
tab3Mat

# tmp_long <- gather(tmp, HLAgene, match, HLAAmatch:HLADQB1match, factor_key=TRUE)
# 
# library(lme4)
# #age_diff + age_sum + f.31.0.0 + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
# #PC7dif + PC8dif + PC9dif + PC10dif +
# LM <- glmer(outcome ~ match +  (1 | uniqueCode), data=tmp_long, family="binomial")
# summary(LM)

out <- data.frame(matrix(0, 6, 4))
colnames(out) <- c("HLAgene", "beta", "se", "pval")

HLAgene <- colnames(couples)[grepl("match", colnames(couples))]

for(i in seq(1, 6)){
  LM <- glm(paste0("outcome ~ ",HLAgene[i]," + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif"), data=tmp, family="binomial")
  out$HLAgene[i] <- HLAgene[i]
  out$beta[i] <- summary(LM)$coefficients[2, 1]
  out$se[i] <- summary(LM)$coefficients[2, 2]
  out$pval[i] <- summary(LM)$coefficients[2, 4]
}

write.table(out, file="onefield_regression_WB_20230414.tsv", col.names = T, row.names = F, sep="\t", quote = F)

tmp_cor <- tmp %>% dplyr::select(colnames(couples)[grepl("match", colnames(couples))])
library(corrplot)
library(Hmisc)
M <- Hmisc::rcorr(as.matrix(tmp_cor))
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

saveRDS(M, file="HLA_onefield_correlation_WB_20230414.rds")


library(metafor)

U <- M$r
vi <- out$se
yi <- out$beta

V <- U*matrix(vi)%*%t(matrix(vi))

metaresult <- summary(rma.mv(yi,V))

res <- rma.mv(yi,V)
png("onefield.meta_WB_20230414.png", width=500, height = 400)
forest(res, atransf=exp, slab = out$HLAgene)
dev.off()
save(tab3Mat, out, M, metaresult, res, file = "onefield_WB_20230414.Rdata")

# for(i in seq(1, 6)){
#   LM <- glm(paste0("outcome ~ ",HLAgene[i]," + age_diff + age_sum + f.31.0.0"), data=tmp, family="binomial")
#   out$HLAgene[i] <- HLAgene[i]
#   out$beta[i] <- summary(LM)$coefficients[2, 1]
#   out$se[i] <- summary(LM)$coefficients[2, 2]
#   out$pval[i] <- summary(LM)$coefficients[2, 4]
# }
# 
# write.table(out, file="onefield_regression_woPCdif.tsv", col.names = T, row.names = F, sep="\t", quote = F)

#SAS
couples <- readRDS("couples_step4_SAS_withHLA_20230414.rds")
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[24:85]), .funs=funs(paste0(unlist(strsplit(.,":"))[1])))

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
couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(second_sex = f.31.0.0[2])

tmp <- couples %>% dplyr::select(uniqueCode, age_diff, age_sum, f.31.0.0, second_sex, colnames(couples)[grepl("match", colnames(couples))], coinfection,
                                 colnames(couples)[grepl("dif$", colnames(couples))]) %>% unique()

tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))
tmp <- tmp[!duplicated(tmp$uniqueCode),]
library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
tableOne <- CreateTableOne(vars = catVars, strata = c("outcome"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
write.csv(tab3Mat, file = "Table_onefield_SAS_20230414.csv")
tab3Mat <- read.csv("Table_onefield_SAS_20230414.csv")
tab3Mat


out <- data.frame(matrix(0, 6, 4))
colnames(out) <- c("HLAgene", "beta", "se", "pval")

HLAgene <- colnames(couples)[grepl("match", colnames(couples))]

for(i in seq(1, 6)){
  LM <- glm(paste0("outcome ~ ",HLAgene[i]," + age_diff + age_sum + f.31.0.0 + second_sex + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif"), data=tmp, family="binomial")
  out$HLAgene[i] <- HLAgene[i]
  out$beta[i] <- summary(LM)$coefficients[2, 1]
  out$se[i] <- summary(LM)$coefficients[2, 2]
  out$pval[i] <- summary(LM)$coefficients[2, 4]
}

write.table(out, file="onefield_regression_SAS_2023044.tsv", col.names = T, row.names = F, sep="\t", quote = F)

tmp_cor <- tmp %>% dplyr::select(colnames(tmp)[grepl("match", colnames(tmp))])
tmp_cor <- as.matrix(tmp_cor)
tmp_cor <- tmp_cor %>% dplyr::select(colnames(tmp_cor)[grepl("match", colnames(tmp_cor))])

library(corrplot)
library(Hmisc)
tmp_cor <- tmp_cor[!is.na(tmp_cor)]
M <- Hmisc::rcorr(as.matrix(tmp_cor))
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

saveRDS(M, file="HLA_onefield_correlation_SAS_2023044.rds")


library(metafor)

U <- M$r
vi <- out$se
yi <- out$beta

V <- U*matrix(vi)%*%t(matrix(vi))

metaresult <- summary(rma.mv(yi,V))

res <- rma.mv(yi,V)
png("onefield.meta_SAS_2023044.png", width=500, height = 400)
forest(res, atransf=exp, slab = out$HLAgene)
dev.off()
save(tab3Mat, out, M, metaresult, res, file = "onefield_SAS_2023044.Rdata")

#forest(res1, atransf=exp, slab = out$HLAgene)
##meta
load("onefield_WB_20230414.Rdata")
res_WB <- res
load("onefield_SAS_2023044.Rdata")
res_SAS <- res

res_WB <- data.frame(cbind(res_WB$b, res_WB$se, res_WB$pval))
colnames(res_WB) <- c("beta", "se", "pval")
res_WB <- res_WB %>% mutate(ancestry = "European")
res_SAS <- data.frame(cbind(res_SAS$b, res_SAS$se, res_SAS$pval))
colnames(res_SAS) <- c("beta", "se", "pval")
res_SAS <- res_SAS %>% mutate(ancestry = "South Asian")

res <- bind_rows(res_WB, res_SAS)
library(meta)
m1 <- metagen(beta,
              se,
              data=res,
              studlab=paste(ancestry),
              fixed = TRUE,
              random = TRUE,
              prediction=FALSE,
              sm="OR")

forest(m1, smlab="",leftcols=c("studlab"),
       rightcols=c("effect", "ci"),print.I2.ci = TRUE,print.tau2 = FALSE,
       leftlabs = c("Ancestry"),zero.pval = TRUE)
