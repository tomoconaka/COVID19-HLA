setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")
library(data.table)
library(tidyr)
library(dplyr)
couples <- readRDS("couples_step3.rds")
couples1 <- couples %>% group_by(uniqueCode) %>%
  filter(coinfection == 1 | (sum(infection) >= 1 & death == 0))

get_first_covid19_test_date <- function(x, y){
  tmp1 <- unique(unlist(x))
  tmp2 <- unique(unlist(y))
  tmp3 <- data.frame(expand.grid(tmp1, tmp2))
  tmp3 <- tmp3 %>% mutate(interval = abs(as.Date(Var1, format="%d/%m/%Y") - as.Date(Var2, format="%d/%m/%Y")))
  return(min(as.Date(tmp3$Var1[tmp3$interval == as.numeric(min(tmp3$interval))], format="%d/%m/%Y"), as.Date(tmp3$Var2[tmp3$interval == as.numeric(min(tmp3$interval))], format="%d/%m/%Y")))
}

couples1 <- couples1 %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

couples1 <- couples1 %>% filter(!(bothinfected == TRUE & interval == 0))
couples1 <- couples1 %>% 
  group_by(uniqueCode) %>% 
  mutate(first_covid19_test_date = case_when(bothinfected == TRUE ~ get_first_covid19_test_date(positive_date[1],positive_date[2]),
                                             !is.na(positive_date[1]) ~ min(as.Date(unlist(positive_date[1]), format="%d/%m/%Y")),
                                             !is.na(positive_date[2]) ~ min(as.Date(unlist(positive_date[2]), format="%d/%m/%Y")),
                                             ))

couples1 <- couples1 %>% 
  group_by(f.eid) %>% 
  mutate(first = case_when(bothinfected == FALSE ~ infection,
                           first_covid19_test_date %in% as.Date(unlist(positive_date), format="%d/%m/%Y") ~ 1, 
                           !(first_covid19_test_date %in% as.Date(unlist(positive_date), format="%d/%m/%Y")) ~ 0))

saveRDS(couples1, file="couples_step4_whichfirst.rds")

# female_first <- couples1 %>%
#   filter(f.31.0.0 == 0 & first == 1)
# sum(female_first$coinfection, na.rm=T)/dim(female_first)[1]
# male_first <- couples1  %>%
#   filter(f.31.0.0 == 1 & first == 1)
# sum(male_first$coinfection, na.rm=T)/dim(male_first)[1]

