setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")
library(data.table)
library(dplyr)
library(tidyr)

couples_ori <- readRDS("couples_step4_WB_withHLA.rds")
couples <- couples_ori %>% arrange(desc(first))
couples1 <- couples %>% group_by(uniqueCode) %>% 
  mutate(PC1dif = abs(PC1[1] - PC1[2]),
         age_diff = abs(f.21022.0.0[1]- f.21022.0.0[2]),
         age_sum = f.21022.0.0[1] + f.21022.0.0[2]) %>% distinct(uniqueCode, .keep_all=TRUE)

couples1 <- couples1 %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                                   TRUE ~ 0))

quantile(couples1$PC1dif, probs = seq(.1, .9, by = .1))

out <- data.frame(matrix(0, 9, 8))
colnames(out) <- c("PC1diff_Cutoff", "beta", "se", "pval", "N_match_SA", "N_match_NSA", "N_unmatch_SA", "N_unmatch_NSA")

for(i in seq(1,9)){
  tmp <- couples1 %>% mutate(match = ifelse(PC1dif <= quantile(couples1$PC1dif, probs = seq(.1, .9, by = .1))[i], 1, 0))
  LM <- glm(paste0("outcome ~ match + age_diff + age_sum + f.31.0.0"), data=tmp, family="binomial")
  out$PC1diff_Cutoff[i] <- quantile(couples1$PC1dif, probs = seq(.1, .9, by = .1))[i]
  out$beta[i] <- summary(LM)$coefficients[2, 1]
  out$se[i] <- summary(LM)$coefficients[2, 2]
  out$pval[i] <- summary(LM)$coefficients[2, 4]
  out$N_match_SA[i] <- sum(tmp$match[tmp$outcome == 1])
  out$N_match_NSA[i] <- sum(tmp$match[tmp$outcome == 0])
  out$N_unmatch_SA[i] <- sum(tmp$outcome[tmp$match == 0])
  out$N_unmatch_NSA[i] <- dim(tmp)[1] - sum(tmp$match[tmp$outcome == 1]) - sum(tmp$match[tmp$outcome == 0]) - sum(tmp$outcome[tmp$match == 0])
}

library(openxlsx)
out %>% write.xlsx("PC1.xlsx")

