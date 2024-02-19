setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

#https://www.frontiersin.org/files/Articles/799519/fmicb-12-799519-HTML/image_m/fmicb-12-799519-g001.jpg
#https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=23165

couples <- readRDS("couples_step4_WB_withHLA.rds")

ABO <- fread("WB.couples.ABO.raw")
ABO <- ABO %>% mutate(ABO = case_when(rs8176719_T == 2 ~ "OO",
                                      rs8176746_G == 0 ~ "BB",
                                      rs8176746_G == 1 & rs8176719_T == 1 ~ "BO",
                                      rs8176746_G == 1 ~ "AB",
                                      rs8176719_T == 1 ~ "AO",
                                      TRUE ~ "AA"))
ABO <- ABO %>% select(IID, ABO)
couples <- couples %>% inner_join(ABO, by=c("f.eid"="IID"))

couples <- couples %>% arrange(desc(first))

couples <- couples %>% group_by(uniqueCode) %>%
  mutate(ABOmatch = case_when(ABO[1] == ABO[2] ~ 1,
                              ABO[1] == "OO" ~ 1,
                              ABO[2] == "AB" ~ 1,
                              TRUE ~ 0
                              ))

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

tmp <- couples %>% dplyr::select(uniqueCode, age_diff, age_sum, f.31.0.0, colnames(couples)[grepl("match", colnames(couples))], coinfection,
                                 colnames(couples)[grepl("dif$", colnames(couples))])
tmp <- tmp %>% mutate(outcome = case_when(coinfection == TRUE ~ 1,
                                          TRUE ~ 0))
tmp <- tmp[!duplicated(tmp$uniqueCode),]

library(tableone)
tmp <- data.frame(tmp)
catVars <- colnames(tmp)[grepl("match", colnames(tmp))]
tableOne <- CreateTableOne(vars = catVars, strata = c("outcome"), data = tmp,
                           factorVars = catVars )
tab3Mat <- print(tableOne, showAllLevels = TRUE, quote=F, noSpaces=TRUE, exact = catVars)
write.csv(tab3Mat, file = "Table_onefield_ABO.csv")

LM <- glm(paste0("outcome ~ ABOmatch + age_diff + age_sum + f.31.0.0 + PC1dif + PC2dif + PC3dif + PC4dif + PC5dif +PC6dif +
                   PC7dif + PC8dif + PC9dif + PC10dif"), data=tmp, family="binomial")
summary(LM)



                                   