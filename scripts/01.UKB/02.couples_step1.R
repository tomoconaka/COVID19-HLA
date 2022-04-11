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

household1 <- household %>% drop_na(f.20074.0.0, f.20075.0.0, f.54.0.0, f.670.0.0, f.680.0.0, f.699.0.0, f.709.0.0, f.728.0.0) %>%
  mutate(uniqueCode = paste0("1_",f.20074.0.0, "_",f.20075.0.0, "_",f.54.0.0,"_",f.670.0.0, "_",f.680.0.0, "_",f.699.0.0, "_",f.709.0.0, "_",f.728.0.0))

couples1 <- household1 %>% group_by(uniqueCode) %>%
  filter(f.6141.0.0 == 1 | f.6141.0.1 == 1 | f.6141.0.2 == 1 | f.6141.0.3 == 1 | f.6141.0.4 == 1 ) %>% ungroup()

couples1 %>% group_by(uniqueCode) %>% summarise(length = length(f.eid)) %>% count(length)

couples1 <- couples1 %>% 
  group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% 
  ungroup()

dim(couples1)
# couples1 <- couples1 %>%
#   group_by(uniqueCode) %>%
#   filter(max(f.2734.0.0, na.rm=T) == max(f.2405.0.0, na.rm=T)) %>%
#   ungroup()

length(unique(couples1$uniqueCode)) #57219

couples1 <- couples1 %>% dplyr::select(f.eid, f.31.0.0, f.21022.0.0, uniqueCode)

household <- household %>% drop_na(f.20074.1.0, f.20075.1.0, f.54.1.0, f.670.1.0, f.680.1.0, f.699.1.0, f.709.1.0, f.728.1.0) %>%
  mutate(uniqueCode = paste0("2_",f.20074.1.0, "_",f.20075.1.0, "_",f.54.1.0,"_",f.670.1.0, "_",f.680.1.0, "_",f.699.1.0, "_",f.709.1.0, "_",f.728.1.0))

couples2 <- household %>% group_by(uniqueCode) %>%
  filter((f.6141.1.0 == 1 | f.6141.1.1 == 1 | f.6141.1.2 == 1 | f.6141.1.3 == 1 |
            f.6141.1.4 == 1 ) & !(f.eid %in% couples1$f.eid)) %>% ungroup()

couples2 %>% group_by(uniqueCode) %>% summarise(length = length(f.eid)) %>% count(length)

couples2 <- couples2 %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()
length(unique(couples2$uniqueCode)) #844
# couples2 <- couples2 %>%
#   group_by(uniqueCode) %>%
#   filter(max(f.2734.1.0, na.rm=T) == max(f.2405.1.0, na.rm=T)) %>%
#   ungroup()

couples2 <- couples2 %>% dplyr::select(f.eid, f.31.0.0, f.21022.0.0, uniqueCode)

household <- household %>% drop_na(f.20074.2.0, f.20075.2.0, f.54.2.0, f.670.2.0, f.680.2.0, f.699.2.0, f.709.2.0, f.728.2.0) %>%
  mutate(uniqueCode = paste0("3_",f.20074.2.0, "_",f.20075.2.0, "_",f.54.2.0,"_",f.670.2.0, "_",f.680.2.0, "_",f.699.2.0, "_",f.709.2.0, "_",f.728.2.0))

couples3 <- household %>% group_by(uniqueCode) %>%
  filter((f.6141.2.0 == 1 | f.6141.2.1 == 1 | f.6141.2.2 == 1 | f.6141.2.3 == 1 |
            f.6141.2.4 == 1 ) & !(f.eid %in% couples1$f.eid) & !(f.eid %in% couples2$f.eid)) %>% ungroup()

couples3 %>% group_by(uniqueCode) %>% summarise(length = length(f.eid)) %>% count(length)

couples3 <- couples3 %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()
length(unique(couples3$uniqueCode))#65
# couples3 <- couples3 %>%
#   group_by(uniqueCode) %>%
#   filter(max(f.2734.2.0, na.rm=T) == max(f.2405.2.0, na.rm=T)) %>%
#   ungroup()


couples3 <- couples3 %>% dplyr::select(f.eid, f.31.0.0, f.21022.0.0, uniqueCode)

couples <- bind_rows(couples1, couples2, couples3)

length(unique(couples$uniqueCode)) #58128

# remove withdrawal
w <- read.csv("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/w27449_20220222.csv")
colnames(w) <- "ID"
couples <- couples %>% filter(!(f.eid %in% w$ID))

couples1 <- couples %>% group_by(uniqueCode) %>%
  filter(f.31.0.0 %in% c(0,1)) %>% 
  filter(length(unique(f.31.0.0)) == 2) %>% ungroup()
couples1 <- couples1 %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

saveRDS(couples1, file="couples_step1.rds")
