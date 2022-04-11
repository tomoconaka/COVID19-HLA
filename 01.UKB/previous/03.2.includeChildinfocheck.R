setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(data.table)
library(tidyr)
library(dplyr)
household <- fread("ukb27449_20688_household.tab.gz")

# Identification of couple pairs
# Participants were assigned to couple pairs on the basis of a shared household identifier. 
# Individuals who shared a household, 
# reported living in a household with two individuals, 
# and who reported living with a husband, wife, or partner were selected. 
# Any couples with an age gap of >10 years were removed, as were couples whose parental ages matched for either parent. After further selecting White British unrelated individuals from this group there were 34,987 opposite-sex pairs available for analysis. A total of 407 same-sex couples were also identified using the above algorithm. Due to the lower number of same-sex pairs genetic associations were not analysed in these individuals although phenotypic associations were estimated.

household <- household %>% mutate(uniqueCode = paste0("1_",f.20074.0.0, "_",f.20075.0.0, "_",f.54.0.0,"_",f.670.0.0, "_",f.680.0.0, "_",f.699.0.0, "_",f.709.0.0, "_",f.728.0.0))

couples1 <- household %>% group_by(uniqueCode) %>%
  filter(f.6141.0.0 == 1 | f.6141.0.1 == 1 | f.6141.0.2 == 1 | f.6141.0.3 == 1 | f.6141.0.4 == 1 ) %>% ungroup()

couples1 %>% group_by(uniqueCode) %>% summarise(length = length(f.eid)) %>% count(length)

couples1 <- couples1 %>% 
  group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% 
  ungroup()

couples1 <- couples1 %>% 
  group_by(uniqueCode) %>%
  filter(max(f.2734.0.0, na.rm=T) != max(f.2405.0.0, na.rm=T)) %>% 
  ungroup()

length(unique(couples1$uniqueCode))

couples1 <- couples1 %>% dplyr::select(f.eid, f.31.0.0, f.21022.0.0, uniqueCode)

household <- household %>% mutate(uniqueCode = paste0("2_",f.20074.1.0, "_",f.20075.1.0, "_",f.54.1.0,"_",f.670.1.0, "_",f.680.1.0, "_",f.699.1.0, "_",f.709.1.0, "_",f.728.1.0))

couples2 <- household %>% group_by(uniqueCode) %>%
  filter((f.6141.1.0 == 1 | f.6141.1.1 == 1 | f.6141.1.2 == 1 | f.6141.1.3 == 1 |
            f.6141.1.4 == 1 ) & !(f.eid %in% couples1$f.eid)) %>% ungroup()

couples2 %>% group_by(uniqueCode) %>% summarise(length = length(f.eid)) %>% count(length)

couples2 <- couples2 %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()
length(unique(couples2$uniqueCode))
couples2 <- couples2 %>% 
  group_by(uniqueCode) %>%
  filter(max(f.2734.1.0, na.rm=T) != max(f.2405.1.0, na.rm=T)) %>% 
  ungroup()

couples2 <- couples2 %>% dplyr::select(f.eid, f.31.0.0, f.21022.0.0, uniqueCode)

household <- household %>% mutate(uniqueCode = paste0("3_",f.20074.2.0, "_",f.20075.2.0, "_",f.54.2.0,"_",f.670.2.0, "_",f.680.2.0, "_",f.699.2.0, "_",f.709.2.0, "_",f.728.2.0))

couples3 <- household %>% group_by(uniqueCode) %>%
  filter((f.6141.2.0 == 1 | f.6141.2.1 == 1 | f.6141.2.2 == 1 | f.6141.2.3 == 1 |
            f.6141.2.4 == 1 ) & !(f.eid %in% couples1$f.eid) & !(f.eid %in% couples2$f.eid)) %>% ungroup()

couples3 %>% group_by(uniqueCode) %>% summarise(length = length(f.eid)) %>% count(length)

couples3 <- couples3 %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()
length(unique(couples3$uniqueCode))
couples3 <- couples3 %>% 
  group_by(uniqueCode) %>%
  filter(max(f.2734.2.0, na.rm=T) != max(f.2405.2.0, na.rm=T)) %>% 
  ungroup()

couples3 <- couples3 %>% dplyr::select(f.eid, f.31.0.0, f.21022.0.0, uniqueCode)

couples <- bind_rows(couples1, couples2, couples3)

length(unique(couples$uniqueCode))

#117878 - 107070 > 10808 were removed

# remove withdrawal
w <- read.csv("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/w27449_20220222.csv")
colnames(w) <- "ID"
couples <- couples %>% filter(!(f.eid %in% w$ID))

couples <- couples %>% group_by(uniqueCode) %>%
  filter(f.31.0.0 %in% c(0,1)) %>% 
  filter(length(unique(f.31.0.0)) == 2) %>% ungroup()
couples <- couples1 %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

covid1 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_england_20210930.txt.gz") %>% dplyr::select(eid, specdate,result)
covid2 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_scotland_20210831.txt.gz") %>% dplyr::select(eid, specdate, result)
covid3 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_wales_20210831.txt.gz") %>% dplyr::select(eid, specdate, result)

covid <- bind_rows(covid1, covid2, covid3)
couples <- couples %>% arrange(uniqueCode)
couples <- couples %>% group_by(f.eid) %>%
  mutate(infection = ifelse(any(covid$result[covid$eid %in% f.eid]) == 1, 1, 0),
         positive_date = list(covid$specdate[covid$eid %in% f.eid & covid$result == 1]),
         negative_date = list(covid$specdate[covid$eid %in% f.eid & covid$result == 0]))

##remove deceased 

death <- fread("/scratch/richards/tomoko.nakanishi/DATA/UKB/death_20210930.txt.gz")
death <- death %>% filter(as.Date(date_of_death, format="%d/%m/%Y") < as.Date("2020-01-13"))
couples1 <- couples %>% group_by(uniqueCode) %>%
 filter(sum(infection) >= 1)
couples1 <- couples1 %>% group_by(uniqueCode) %>% 
  mutate(bothinfected = ifelse(sum(infection) == 2, TRUE, FALSE)) %>% ungroup()
couples1 <- couples1 %>% mutate(death = ifelse(f.eid %in% death$eid, 1, 0))

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
# couples1 <- couples %>% group_by(uniqueCode) %>%
#  filter(sum(infection) >= 1)

tmp <- couples1 %>% dplyr::select(uniqueCode, interval, bothinfected) %>% dplyr::filter(!is.na(interval) & bothinfected == TRUE) %>% unique()
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
saveRDS(couples1, file="couples_step3.rds")

couples1 <- couples1 %>% mutate(coinfection = case_when(bothinfected == TRUE & interval <= 7 ~ TRUE,
                                                        bothinfected == FALSE & interval <= 14 ~ FALSE,
                                                        bothinfected == TRUE & interval > 15 ~ FALSE))
table(couples1$coinfection)
saveRDS(couples1, file="couples_step3_incl_interval0.rds")

couples1 <- couples1 %>% mutate(coinfection = case_when(bothinfected == TRUE & interval == 0 ~ TRUE,
                                                        bothinfected == FALSE & interval <= 14 ~ FALSE))
table(couples1$coinfection)
saveRDS(couples1, file="couples_step3_interval0_only.rds")



