setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

# family
couples <- readRDS(file="couples_step4_WB_withHLA.rds")

## family
mapping1 <- read.xlsx("Supertypes HLA 1.xlsx", sheet="Sheet2")
mapping1 <- mapping1 %>% dplyr::select(Allele, Family)
mapping2 <- read.xlsx("misssingAllele.xlsx", sheet="Sheet 1")
mapping2 <- mapping2 %>% dplyr::select(Supertype_1, X3) %>% 
  rename(Allele = Supertype_1,
         Family = X3)
mapping <- bind_rows(mapping1, mapping2)

couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[34:95]), .funs=funs(paste0(unlist(strsplit(.,":"))[1], ":", unlist(strsplit(.,":"))[2])))
couples <- couples %>% group_by(f.eid) %>% mutate_at(.vars=vars(colnames(couples)[34:95]), .funs=funs(gsub("HLA-", "", .)))

couples1 <- couples %>% dplyr::select(f.eid, colnames(couples)[34:95])

setDT(mapping)
setDT(couples1)

couples1 <- couples1[ , melt(.SD, id.vars = 'f.eid')     
                      # merging
][mapping, new_value := i.Family, on = c(value = 'Allele') 
  #reform back to original shape
][ , dcast(.SD, f.eid ~ variable, value.var = 'new_value')]

# tmp <- c(couples$A_1, couples$A_2,
#                  couples$B_1, couples$B_2,
#                  couples$C_1, couples$C_1,
#                  couples$DRB1_1, couples$DRB1_2,
#                  couples$DQB1_1, couples$DQB1_2)
# tmp <- unique(tmp)
# 
# tmp <- data.frame(tmp)
# tmp <- tmp %>% filter(!(tmp %in% mapping$Allele))
# write.xlsx(tmp, "misssingAllele.xlsx")

couples2 <- couples %>% dplyr::select(f.eid, !(colnames(couples)[34:95]))

couples <- inner_join(couples1, couples2, by="f.eid")

couples <- couples %>% drop_na(A_1, A_2, B_1, B_2, C_1, C_2, DRB1_1, DRB1_2, DQB1_1, DQB1_2)

couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

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
write.csv(tab3Mat, file = "Table_supertypeHLA1.csv")

out <- data.frame(matrix(0, 5, 4))
colnames(out) <- c("HLAgene", "beta", "se", "pval")

HLAgene <- colnames(tmp)[grepl("match", colnames(tmp))]

for(i in seq(1, 5)){
  LM <- glm(paste0("outcome ~ ",HLAgene[i]," + age_diff + age_sum + f.31.0.0 + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif"), data=tmp, family="binomial")
  out$HLAgene[i] <- HLAgene[i]
  out$beta[i] <- summary(LM)$coefficients[2, 1]
  out$se[i] <- summary(LM)$coefficients[2, 2]
  out$pval[i] <- summary(LM)$coefficients[2, 4]
}

out

library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
tableOne <- CreateTableOne(vars = c(catVars, "age_diff"), strata = c("outcome"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
write.csv(tab3Mat, file = "Table_supertypeHLA1.csv")

tmp <- tmp %>% dplyr::select(colnames(tmp)[c(5:8,10)])

library(corrplot)
M <- rcorr(as.matrix(tmp))
library(Hmisc)
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

saveRDS(M, file="HLA_supertypeHLA1_correlation.rds")

library(metafor)

U <- M$r
vi <- out$se
yi <- out$beta

V <- U*matrix(vi)%*%t(matrix(vi))

metaresult <- summary(rma.mv(yi,V))

save(tab3Mat, out, M, metaresult, file = "supertypeHLA1.Rdata")


# family
couples <- readRDS(file="couples_step4_WB_withHLA.rds")

## family
mapping1 <- read.xlsx("Supertypes HLA 2.xlsx", sheet="Sheet2")
mapping1 <- mapping1 %>% dplyr::select(Allele, Family)
mapping2 <- read.xlsx("misssingAllele.xlsx", sheet="Sheet 1")
mapping2 <- mapping2 %>% dplyr::select(Supertype_2, X5) %>% 
  rename(Allele = Supertype_2,
         Family = X5)
mapping <- bind_rows(mapping1, mapping2)

couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[34:95]), .funs=funs(paste0(unlist(strsplit(.,":"))[1], ":", unlist(strsplit(.,":"))[2])))
couples <- couples %>% group_by(f.eid) %>% mutate_at(.vars=vars(colnames(couples)[34:95]), .funs=funs(gsub("HLA-", "", .)))

couples1 <- couples %>% dplyr::select(f.eid, colnames(couples)[34:95])

setDT(mapping)
setDT(couples1)

couples1 <- couples1[ , melt(.SD, id.vars = 'f.eid')     
                      # merging
][mapping, new_value := i.Family, on = c(value = 'Allele') 
  #reform back to original shape
][ , dcast(.SD, f.eid ~ variable, value.var = 'new_value')]

# tmp <- c(couples$A_1, couples$A_2,
#                  couples$B_1, couples$B_2,
#                  couples$C_1, couples$C_1,
#                  couples$DRB1_1, couples$DRB1_2,
#                  couples$DQB1_1, couples$DQB1_2)
# tmp <- unique(tmp)
# 
# tmp <- data.frame(tmp)
# tmp <- tmp %>% filter(!(tmp %in% mapping$Allele))
# write.xlsx(tmp, "misssingAllele.xlsx")

couples2 <- couples %>% dplyr::select(f.eid, !(colnames(couples)[34:95]))

couples <- inner_join(couples1, couples2, by="f.eid")

couples <- couples %>% drop_na(A_1, A_2, B_1, B_2, C_1, C_2, DRB1_1, DRB1_2, DQB1_1, DQB1_2)

couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

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
write.csv(tab3Mat, file = "Table_supertypeHLA2.csv")

out <- data.frame(matrix(0, 5, 4))
colnames(out) <- c("HLAgene", "beta", "se", "pval")

HLAgene <- colnames(tmp)[grepl("match", colnames(tmp))]

for(i in seq(1, 5)){
  LM <- glm(paste0("outcome ~ ",HLAgene[i]," + age_diff + age_sum + f.31.0.0 + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif"), data=tmp, family="binomial")
  out$HLAgene[i] <- HLAgene[i]
  out$beta[i] <- summary(LM)$coefficients[2, 1]
  out$se[i] <- summary(LM)$coefficients[2, 2]
  out$pval[i] <- summary(LM)$coefficients[2, 4]
}

out

tmp <- tmp %>% dplyr::select(colnames(tmp)[c(5:9)])

library(corrplot)
M <- rcorr(as.matrix(tmp))
library(Hmisc)
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

saveRDS(M, file="HLA_supertypeHLA2_correlation.rds")

library(metafor)

U <- M$r
vi <- out$se
yi <- out$beta

V <- U*matrix(vi)%*%t(matrix(vi))

metaresult <- summary(rma.mv(yi,V))

save(tab3Mat, out, M, metaresult, file = "supertypeHLA2.Rdata")

