setwd("/home/richards/tomoko.nakanishi/09.COVID19/scratch/01.UKBB/02.LongCOVID")
library(data.table)
library(tidyr)
library(dplyr)
covid1 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_england_20210930.txt.gz") %>% select(eid, specdate,result)
covid2 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_scotland_20210831.txt.gz") %>% select(eid, specdate, result)
covid3 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_wales_20210831.txt.gz") %>% select(eid, specdate, result)

covid <- bind_rows(covid1, covid2, covid3)
covid <- covid %>% filter(as.Date(specdate, format="%d/%m/%Y") < as.Date("2021-01-01"))
couples <- readRDS(file="couples_step1.rds")
couples <- couples %>% arrange(uniqueCode)
couples <- couples %>% group_by(f.eid) %>%
  mutate(infection = ifelse(any(covid$result[covid$eid %in% f.eid]) == 1, 1, 0),
         positive_date = list(covid$specdate[covid$eid %in% f.eid & covid$result == 1]),
         negative_date = list(covid$specdate[covid$eid %in% f.eid & covid$result == 0]))
saveRDS(couples, file="couples_step2_2020.rds")
couples <- readRDS("couples_step2_2020.rds")

couples1 <- couples %>% group_by(uniqueCode) %>%
  filter(sum(infection) >= 1)
couples1 <- couples1 %>% group_by(uniqueCode) %>% 
  mutate(bothinfected = ifelse(sum(infection) == 2, TRUE, FALSE)) %>% ungroup()


# tmp <- couples1 %>% filter(bothinfected == TRUE) %>% select(positive_date)
# tmp1 <- unique(unlist(tmp$positive_date[1]))
# tmp2 <- unique(unlist(tmp$positive_date[2]))
# tmp3 <- data.frame(expand.grid(tmp1, tmp2))
# tmp3 <- tmp3 %>% mutate(interval = abs(as.Date(Var1, format="%d/%m/%Y") - as.Date(Var2, format="%d/%m/%Y")))

calculate_minimum_interval <- function(x, y){
  tmp1 <- unique(unlist(x))
  tmp2 <- unique(unlist(y))
  tmp3 <- data.frame(expand.grid(tmp1, tmp2))
  tmp3 <- tmp3 %>% mutate(interval = abs(as.Date(Var1, format="%d/%m/%Y") - as.Date(Var2, format="%d/%m/%Y")))
  return(as.numeric(min(tmp3$interval)))
}

couples1 <- couples1 %>% 
  group_by(uniqueCode) %>% 
  mutate(interval = case_when(bothinfected == TRUE ~ calculate_minimum_interval(positive_date[1],positive_date[2]),
                              bothinfected == FALSE & infection[1] == 1 & infection[2] == 0 & length(unlist(negative_date[2])) != 0 ~ calculate_minimum_interval(positive_date[1],negative_date[2]),
                              bothinfected == FALSE & infection[1] == 0 & infection[2] == 1 & length(unlist(negative_date[1])) != 0 ~ calculate_minimum_interval(negative_date[1],positive_date[2]),
                              TRUE ~ as.numeric(NA))) %>% 
  ungroup()

head(couples1)

tmp <- couples1 %>% select(uniqueCode, interval) %>% filter(!is.na(interval) & bothinfected == TRUE) %>% unique()
median(tmp$interval)
mad(tmp$interval)

tmp_rev <- tmp %>% filter(interval <= 365)
ggplot(tmp_rev, aes(y=interval)) + geom_histogram(bins = 365/7) 
tmp_rev <- tmp %>% filter(interval <= 14)
ggplot(tmp_rev, aes(y=interval)) + geom_histogram(bins = 15) 

couples1 <- couples1 %>% mutate(coinfection = case_when(bothinfected == TRUE & interval <= 7 & interval >= 1 ~ TRUE,
                                                        bothinfected == FALSE & interval <= 14 ~ FALSE,
                                                        bothinfected == TRUE & interval > 15 ~ FALSE))
table(couples1$coinfection)
saveRDS(couples1, file="couples_step3_2020.rds")

