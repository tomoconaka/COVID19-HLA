setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(stringr)
couples <- readRDS(file="couples_step4.rds")

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
death_data <- fread("/scratch/richards/tomoko.nakanishi/DATA/UKB/death_20210930.txt.gz")
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples, file="couples_step4_WB.rds")

write.table(cbind(couples$f.eid, couples$f.eid), file="WB.couples", sep="\t", quote=F, col.names = F, row.names = F)


## two digits
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[12:73]), .funs=funs(paste0(unlist(strsplit(.,":"))[1])))

couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQA1 = list(c(DQA1_1, DQA1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(HLAAmatch = case_when(!any(is.na(unlist(HLAA[1]))) & !any(is.na(unlist(HLAA[2]))) & length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))) == 2 ~ 1,
                               !any(is.na(unlist(HLAA[1]))) & !any(is.na(unlist(HLAA[2]))) & length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))) == 1 & length(unique(unlist(HLAA[1]))) == 1 & length(unique(unlist(HLAA[2]))) == 1 ~ 1,
                               TRUE ~ 0),
         HLABmatch = case_when(!any(is.na(unlist(HLAB[1]))) & !any(is.na(unlist(HLAB[2]))) & length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))) == 2 ~ 1,
                               !any(is.na(unlist(HLAB[1]))) & !any(is.na(unlist(HLAB[2]))) & length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))) == 1 & length(unique(unlist(HLAB[1]))) == 1 & length(unique(unlist(HLAB[2]))) == 1 ~ 1,
                               TRUE ~ 0),
         HLACmatch = case_when(!any(is.na(unlist(HLAC[1]))) & !any(is.na(unlist(HLAC[2]))) & length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))) == 2 ~ 1,
                               !any(is.na(unlist(HLAC[1]))) & !any(is.na(unlist(HLAC[2]))) & length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))) == 1 & length(unique(unlist(HLAC[1]))) == 1 & length(unique(unlist(HLAC[2]))) == 1 ~ 1,
                               TRUE ~ 0),
         HLADRB1match = case_when(!any(is.na(unlist(HLADRB1[1]))) & !any(is.na(unlist(HLADRB1[2]))) & length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))) == 2 ~ 1,
                                  !any(is.na(unlist(HLADRB1[1]))) & !any(is.na(unlist(HLADRB1[2]))) & length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))) == 1 & length(unique(unlist(HLADRB1[1]))) == 1 & length(unique(unlist(HLADRB1[2]))) == 1 ~ 1,
                                  TRUE ~ 0),
         HLADQA1match = case_when(!any(is.na(unlist(HLADQA1[1]))) & !any(is.na(unlist(HLADQA1[2]))) & length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))) == 2 ~ 1,
                                  !any(is.na(unlist(HLADQA1[1]))) & !any(is.na(unlist(HLADQA1[2]))) & length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))) == 1 & length(unique(unlist(HLADQA1[1]))) == 1 & length(unique(unlist(HLADQA1[2]))) == 1 ~ 1,
                                  TRUE ~ 0),
         HLADQB1match = case_when(!any(is.na(unlist(HLADQB1[1]))) & !any(is.na(unlist(HLADQB1[2]))) & length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))) == 2 ~ 1,
                                  !any(is.na(unlist(HLADQB1[1]))) & !any(is.na(unlist(HLADQB1[2]))) & length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))) == 1 & length(unique(unlist(HLADQB1[1]))) == 1 & length(unique(unlist(HLADQB1[2]))) == 1 ~ 1,
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

couples <- couples %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))

tmp <- couples %>% dplyr::select(uniqueCode, PC_distance, age_diff, colnames(couples)[grepl("match", colnames(couples))], coinfection) %>% unique()
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
  LM <- glm(paste0("outcome ~ ",HLAgene[i]," + PC_distance + age_diff"), data=tmp, family="binomial")
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


# family
couples <- readRDS(file="couples_step4.rds")

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
death_data <- fread("/scratch/richards/tomoko.nakanishi/DATA/UKB/death_20210930.txt.gz")
death_data <- death_data %>% filter(as.Date(date_of_death, format="%d/%m/%Y") < as.Date("2020-01-13"))
couples <- couples %>% filter(!(f.eid %in% death_data$eid))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

## family
mapping1 <- read.xlsx("Supertypes HLA 1.xlsx", sheet="Sheet2")
mapping1 <- mapping1 %>% dplyr::select(Allele, Family)
mapping2 <- read.xlsx("misssingAllele.xlsx", sheet="Sheet 1")
mapping2 <- mapping2 %>% dplyr::select(Supertype_1, X3) %>% 
  rename(Allele = Supertype_1,
        Family = X3)
mapping <- bind_rows(mapping1, mapping2)

couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[12:73]), .funs=funs(paste0(unlist(strsplit(.,":"))[1], ":", unlist(strsplit(.,":"))[2])))
couples <- couples %>% group_by(f.eid) %>% mutate_at(.vars=vars(colnames(couples)[12:73]), .funs=funs(gsub("HLA-", "", .)))

couples1 <- couples %>% dplyr::select(f.eid, colnames(couples)[12:73])

setDT(mapping)
setDT(couples)

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

couples2 <- couples %>% dplyr::select(f.eid, !(colnames(couples)[12:73]))

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
# 
# couples <- couples %>% 
#   mutate(HLAclass1match = case_when(HLAAmatch == 1 & HLABmatch == 1 & HLACmatch == 1 ~ 1,
#                                     TRUE ~ 0),
#          HLAclass2match = case_when(HLADRB1match == 1 & HLADQA1match == 1 & HLADQB1match == 1 ~ 1,
#                                     TRUE ~ 0))

couples <- couples %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                                  TRUE ~ 0))

tmp <- couples %>% dplyr::select(uniqueCode, PC_distance, age_diff, colnames(couples)[grepl("match", colnames(couples))], coinfection) %>% unique()
tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))

out <- data.frame(matrix(0, 5, 4))
colnames(out) <- c("HLAgene", "beta", "se", "pval")

HLAgene <- colnames(tmp)[grepl("match", colnames(tmp))]

for(i in seq(1, 5)){
  LM <- glm(paste0("outcome ~ ",HLAgene[i]," + PC_distance + age_diff"), data=tmp, family="binomial")
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

tmp <- tmp %>% dplyr::select(colnames(tmp)[4:8])

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


