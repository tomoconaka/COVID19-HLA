setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")
library(data.table)
library(dplyr)
library(tidyr)

couples_ori <- readRDS("couples_step4_WB_withHLA.rds")
couples <- couples_ori %>% group_by(uniqueCode) %>%
  mutate(ID1 = f.eid[1],
         ID2 = f.eid[2]) %>%
  distinct(uniqueCode, .keep_all=TRUE) %>% select(uniqueCode, ID1, ID2)

# write.table(couples, file="WBcoupleswithHLA", sep="\t", quote=F, col.names = F, row.names = F)

related <- fread("WB.couples-pruned-related-degree13.kin0_filtered")
related <- related %>% select(ID1, ID2, IBS0, Kinship, InfType)
related <- related %>% mutate(Relatedness = case_when(Kinship > 0.177 ~ 1,
                                                      Kinship > 0.0884 ~ 2,
                                                      Kinship > 0.0442 ~ 3,
                                                      Kinship > 0.0442/2 ~ 4,
                                                      Kinship > 0.0442/4 ~ 5,
                                                      Kinship > 0.0442/8 ~ 7,
                                                      Kinship > 0.0442/16 ~ 8,
                                                      Kinship > 0.0442/32 ~ 9,
                                                      Kinship > 0.0442/64 ~ 10,
                                                      Kinship > 0.0442/128 ~ 11,
                                                      Kinship > 0.0442/256 ~ 12,
                                                      TRUE ~ 13
))

library(ggplot2)
ggplot(related, aes(x=IBS0, y=Kinship, col=as.factor(Relatedness))) + geom_point()

couples_ori <- readRDS("couples_step4_WB_withHLA.rds")
couples <- couples_ori %>% arrange(desc(first))
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

couples <- couples %>% group_by(uniqueCode) %>%
  mutate(ID1 = f.eid[1],
         ID2 = f.eid[2]) %>%
  distinct(uniqueCode, .keep_all=TRUE)

couples1 <- couples %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                           TRUE ~ 0))

couples1 <- couples1 %>% left_join(related, by=c("ID1"="ID1", "ID2"="ID2"))
couples1 <- couples1 %>% mutate(Relatedness = replace_na(Relatedness, 14))

out <- data.frame(matrix(0, 8, 8))
colnames(out) <- c("Relatedness", "beta", "se", "pval", "N_match_SA", "N_match_NSA", "N_unmatch_SA", "N_unmatch_NSA")

for(i in seq(6,13)){
  tmp <- couples1 %>% mutate(match = ifelse(Relatedness <= i, 1, 0))
  LM <- glm(paste0("outcome ~ match + age_diff + age_sum + f.31.0.0"), data=tmp, family="binomial")
  out$Relatedness[i-5] <- i
  out$beta[i-5] <- summary(LM)$coefficients[2, 1]
  out$se[i-5] <- summary(LM)$coefficients[2, 2]
  out$pval[i-5] <- summary(LM)$coefficients[2, 4]
  out$N_match_SA[i-5] <- sum(tmp$match[tmp$outcome == 1])
  out$N_match_NSA[i-5] <- sum(tmp$match[tmp$outcome == 0])
  out$N_unmatch_SA[i-5] <- sum(tmp$outcome[tmp$match == 0])
  out$N_unmatch_NSA[i-5] <- dim(tmp)[1] - sum(tmp$match[tmp$outcome == 1]) - sum(tmp$match[tmp$outcome == 0]) - sum(tmp$outcome[tmp$match == 0])
}

library(openxlsx)
out %>% write.xlsx("relatedness.xlsx")

related %>% filter(InfType == "2nd")
related %>% filter(InfType == "3rd")
