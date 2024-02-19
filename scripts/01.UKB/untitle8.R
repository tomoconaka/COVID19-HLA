setwd("/home/richards/tomoko.nakanishi/09.COVID19/scratch/01.UKBB/02.LongCOVID")

library(stringr)
couples <- readRDS(file="couples_step4_2020.rds")

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))

## two digits
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[11:72]), .funs=funs(paste0(unlist(strsplit(.,":"))[1])))

couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQA1 = list(c(DQA1_1, DQA1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(HLAAmatch = case_when(length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))) == 2 ~ 1,
                               TRUE ~ 0),
         HLABmatch = case_when(length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))) == 2 ~ 1,
                               TRUE ~ 0),
         HLACmatch = case_when(length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))) == 2 ~ 1,
                               TRUE ~ 0),
         HLADRB1match = case_when(length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))) == 2 ~ 1,
                                  TRUE ~ 0),
         HLADQA1match = case_when(length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))) == 2 ~ 1,
                                  TRUE ~ 0),
         HLADQB1match = case_when(length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))) == 2 ~ 1,
                                  TRUE ~ 0)
  )

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

couples <- couples %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                                  TRUE ~ 0))

tmp <- couples %>% dplyr::select(uniqueCode, colnames(couples)[grepl("match", colnames(couples))], coinfection) %>% unique()
tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))

library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
tableOne <- CreateTableOne(vars = catVars, strata = c("coinfection"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
# write.csv(tab3Mat, file = "Table_onefield_strict.csv")

tableOne <- CreateTableOne(vars = catVars, strata = c("outcome"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
# write.csv(tab3Mat, file = "Table_onefield_relaxed.csv")

tmp <- tmp %>% dplyr::select(colnames(tmp)[2:7])

library(corrplot)
M <- rcorr(as.matrix(tmp))
library(Hmisc)
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

# saveRDS(M, file="HLA_onefield_correlation.rds")




couples <- readRDS(file="couples_step4_2020.rds")

# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))

## two digits
couples <- couples %>% group_by(f.eid) %>% 
  mutate_at(.vars=vars(colnames(couples)[11:72]), .funs=funs(paste0(unlist(strsplit(.,":"))[1], ":", unlist(strsplit(.,":"))[2])))

couples <- couples %>% group_by(f.eid) %>% 
  mutate(HLAA = list(c(A_1, A_2)),
         HLAB = list(c(B_1, B_2)),
         HLAC = list(c(C_1, C_2)),
         HLADRB1 = list(c(DRB1_1, DRB1_2)),
         HLADQA1 = list(c(DQA1_1, DQA1_2)),
         HLADQB1 = list(c(DQB1_1, DQB1_2))) %>% ungroup()

couples <- couples %>% group_by(uniqueCode) %>% 
  mutate(HLAAmatch = case_when(length(intersect(unlist(HLAA[1]),unlist(HLAA[2]))) == 2 ~ 1,
                               TRUE ~ 0),
         HLABmatch = case_when(length(intersect(unlist(HLAB[1]),unlist(HLAB[2]))) == 2 ~ 1,
                               TRUE ~ 0),
         HLACmatch = case_when(length(intersect(unlist(HLAC[1]),unlist(HLAC[2]))) == 2 ~ 1,
                               TRUE ~ 0),
         HLADRB1match = case_when(length(intersect(unlist(HLADRB1[1]),unlist(HLADRB1[2]))) == 2 ~ 1,
                                  TRUE ~ 0),
         HLADQA1match = case_when(length(intersect(unlist(HLADQA1[1]),unlist(HLADQA1[2]))) == 2 ~ 1,
                                  TRUE ~ 0),
         HLADQB1match = case_when(length(intersect(unlist(HLADQB1[1]),unlist(HLADQB1[2]))) == 2 ~ 1,
                                  TRUE ~ 0)
  )

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

couples <- couples %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                                  TRUE ~ 0))

tmp <- couples %>% dplyr::select(uniqueCode, colnames(couples)[grepl("match", colnames(couples))], coinfection) %>% unique()
tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))

library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
tableOne <- CreateTableOne(vars = catVars, strata = c("coinfection"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
# write.csv(tab3Mat, file = "Table_twofield_strict.csv")

tableOne <- CreateTableOne(vars = catVars, strata = c("outcome"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
# write.csv(tab3Mat, file = "Table_twofield_relaxed.csv")

tmp <- tmp %>% dplyr::select(colnames(tmp)[2:7])

library(corrplot)
M <- rcorr(as.matrix(tmp))
library(Hmisc)
col2 = colorRampPalette(c('#005BB6', 'white', '#C70000'))
corrplot(M$r, method = 'square', tl.col = "black", col = col2(100))

# saveRDS(M, file="HLA_twofield_correlation.rds")



summary(glm(outcome ~ HLAclass2match, data=tmp, family="binomial"))

summary(glm(outcome ~ HLADRB1match, data=tmp, family="binomial"))

summary(glm(outcome ~ HLACmatch, data=tmp, family="binomial"))

summary(glm(outcome ~ HLADRB1match, data=tmp, family="binomial"))

summary(glm(outcome ~ HLADQA1match, data=tmp, family="binomial"))

summary(glm(outcome ~ HLADQB1match, data=tmp, family="binomial"))




library(lme4)
summary(glmer(outcome ~ HLADRB1match + (1 | uniqueCode) + f.31.0.0 + f.21022.0.0, data = couples, family="binomial"))



