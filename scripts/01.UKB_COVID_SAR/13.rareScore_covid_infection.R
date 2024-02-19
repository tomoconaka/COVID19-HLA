setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(data.table)
library(tidyr)
library(dplyr)

pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

batch10 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_10.tsv.gz")
batch11 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_11.tsv.gz")
batch12 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_12.tsv.gz")
batch13 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_13.tsv.gz")
batch14 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_14.tsv.gz")
batch15 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_15.tsv.gz")
batch16 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_16.tsv.gz")
batch17 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_17.tsv.gz")
batch18 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_18.tsv.gz")
batch19 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_19.tsv.gz")
batch20 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_20.tsv.gz")
batch21 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_21.tsv.gz")
batch22 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_22.tsv.gz")
batch23 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_23.tsv.gz")
batch24 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_24.tsv.gz")
batch25 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_25.tsv.gz")
batch26 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_26.tsv.gz")
batch27 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_27.tsv.gz")
batch28 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_28.tsv.gz")
batch29 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_29.tsv.gz")
batch30 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_30.tsv.gz")
batch31 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_31.tsv.gz")
batch32 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_32.tsv.gz")
batch33 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_33.tsv.gz")
batch34 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_34.tsv.gz")
batch35 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_35.tsv.gz")
batch36 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_36.tsv.gz")
batch37 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_37.tsv.gz")
batch38 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_38.tsv.gz")
batch39 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_39.tsv.gz")
batch40 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_40.tsv.gz")
batch41 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_41.tsv.gz")
batch42 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_42.tsv.gz")
batch43 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_43.tsv.gz")
batch44 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_44.tsv.gz")
batch45 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_45.tsv.gz")
batch46 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_46.tsv.gz")
batch47 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_47.tsv.gz")
batch48 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_48.tsv.gz")
batch49 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_49.tsv.gz")
batch50 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_50.tsv.gz")
batch51 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_51.tsv.gz")
batch52 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_52.tsv.gz")
batch53 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_53.tsv.gz")
batch54 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_54.tsv.gz")
batch55 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_55.tsv.gz")
batch56 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_56.tsv.gz")
batch57 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_57.tsv.gz")
batch58 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_58.tsv.gz")
batch59 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_59.tsv.gz")
batch60 <- fread("/scratch/richards/guillaume.butler-laporte/Tomoko/calls/hla_df_batch_60.tsv.gz")

hla_all <- bind_rows(batch10, batch11, batch12, batch13, batch14, batch15, batch16, batch17, batch18, batch19,
                     batch20, batch21, batch22, batch23, batch24, batch25, batch26, batch27, batch28, batch29,
                     batch30, batch31, batch32, batch33, batch34, batch35, batch36, batch37, batch38, batch39,
                     batch40, batch41, batch42, batch43, batch44, batch45, batch46, batch47, batch48, batch49,
                     batch50, batch51, batch52, batch53, batch54, batch55, batch56, batch57, batch58, batch59, batch60)

data <- pc %>% inner_join(hla_all, by=c("ID"))

covid1 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_england_20210930.txt.gz") %>% dplyr::select(eid, specdate,result)
covid2 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_scotland_20210831.txt.gz") %>% dplyr::select(eid, specdate, result)
covid3 <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/covid19_result_wales_20210831.txt.gz") %>% dplyr::select(eid, specdate, result)

covid <- bind_rows(covid1, covid2, covid3)
covid_pos <- covid %>% dplyr::filter(result == 1)
covid_neg <- covid %>% dplyr::filter(result == 0)
# covid2020 <- covid %>% filter(as.Date(specdate, format="%d/%m/%Y") < as.Date("2021-01-01"))
data <- data %>% mutate(covid19 = ifelse(ID %in% covid_pos$eid, 1, 0))

##remove deceased 

death <- fread("/scratch/richards/tomoko.nakanishi/DATA/UKB/death_20210930.txt.gz")
death <- death %>% filter(as.Date(date_of_death, format="%d/%m/%Y") < as.Date("2020-01-13"))
data <- data %>% mutate(death = ifelse(ID %in% death$eid, 1, 0))

data <- data %>% filter(death == 0)

##convert to dosage
data <- data %>% rowwise %>% 
  mutate(A_01 = length(names(.)[22:23][which(grepl("HLA-A\\*01", c_across(A_1:A_2)))]),
         A_02 = length(names(.)[22:23][which(grepl("HLA-A\\*02", c_across(A_1:A_2)))]),
         A_03 = length(names(.)[22:23][which(grepl("HLA-A\\*03", c_across(A_1:A_2)))]),
         A_11 = length(names(.)[22:23][which(grepl("HLA-A\\*11", c_across(A_1:A_2)))]),
         A_23 = length(names(.)[22:23][which(grepl("HLA-A\\*23", c_across(A_1:A_2)))]),
         A_24 = length(names(.)[22:23][which(grepl("HLA-A\\*24", c_across(A_1:A_2)))]),
         A_25 = length(names(.)[22:23][which(grepl("HLA-A\\*25", c_across(A_1:A_2)))]),
         A_26 = length(names(.)[22:23][which(grepl("HLA-A\\*26", c_across(A_1:A_2)))]),
         A_29 = length(names(.)[22:23][which(grepl("HLA-A\\*29", c_across(A_1:A_2)))]),
         A_30 = length(names(.)[22:23][which(grepl("HLA-A\\*30", c_across(A_1:A_2)))]),
         A_31 = length(names(.)[22:23][which(grepl("HLA-A\\*31", c_across(A_1:A_2)))]),
         A_32 = length(names(.)[22:23][which(grepl("HLA-A\\*32", c_across(A_1:A_2)))]),
         A_33 = length(names(.)[22:23][which(grepl("HLA-A\\*33", c_across(A_1:A_2)))]),
         A_34 = length(names(.)[22:23][which(grepl("HLA-A\\*34", c_across(A_1:A_2)))]),
         A_36 = length(names(.)[22:23][which(grepl("HLA-A\\*36", c_across(A_1:A_2)))]),
         A_40 = length(names(.)[22:23][which(grepl("HLA-A\\*40", c_across(A_1:A_2)))]),
         A_43 = length(names(.)[22:23][which(grepl("HLA-A\\*43", c_across(A_1:A_2)))]),
         A_66 = length(names(.)[22:23][which(grepl("HLA-A\\*66", c_across(A_1:A_2)))]),
         A_68 = length(names(.)[22:23][which(grepl("HLA-A\\*68", c_across(A_1:A_2)))]),
         A_69 = length(names(.)[22:23][which(grepl("HLA-A\\*69", c_across(A_1:A_2)))]),
         A_74 = length(names(.)[22:23][which(grepl("HLA-A\\*74", c_across(A_1:A_2)))]),
         A_80 = length(names(.)[22:23][which(grepl("HLA-A\\*80", c_across(A_1:A_2)))])
  )

data <- data %>% rowwise %>% 
  mutate(B_07 = length(names(.)[24:25][which(grepl("HLA-B\\*07", c_across(B_1:B_2)))]),
         B_08 = length(names(.)[24:25][which(grepl("HLA-B\\*08", c_across(B_1:B_2)))]),
         B_13 = length(names(.)[24:25][which(grepl("HLA-B\\*13", c_across(B_1:B_2)))]),
         B_14 = length(names(.)[24:25][which(grepl("HLA-B\\*14", c_across(B_1:B_2)))]),
         B_15 = length(names(.)[24:25][which(grepl("HLA-B\\*15", c_across(B_1:B_2)))]),
         B_18 = length(names(.)[24:25][which(grepl("HLA-B\\*18", c_across(B_1:B_2)))]),
         B_27 = length(names(.)[24:25][which(grepl("HLA-B\\*27", c_across(B_1:B_2)))]),
         B_35 = length(names(.)[24:25][which(grepl("HLA-B\\*35", c_across(B_1:B_2)))]),
         B_37 = length(names(.)[24:25][which(grepl("HLA-B\\*37", c_across(B_1:B_2)))]),
         B_38 = length(names(.)[24:25][which(grepl("HLA-B\\*38", c_across(B_1:B_2)))]),
         B_39 = length(names(.)[24:25][which(grepl("HLA-B\\*39", c_across(B_1:B_2)))]),
         B_40 = length(names(.)[24:25][which(grepl("HLA-B\\*40", c_across(B_1:B_2)))]),
         B_41 = length(names(.)[24:25][which(grepl("HLA-B\\*41", c_across(B_1:B_2)))]),
         B_42 = length(names(.)[24:25][which(grepl("HLA-B\\*42", c_across(B_1:B_2)))]),
         B_44 = length(names(.)[24:25][which(grepl("HLA-B\\*44", c_across(B_1:B_2)))]),
         B_45 = length(names(.)[24:25][which(grepl("HLA-B\\*45", c_across(B_1:B_2)))]),
         B_46 = length(names(.)[24:25][which(grepl("HLA-B\\*46", c_across(B_1:B_2)))]),
         B_47 = length(names(.)[24:25][which(grepl("HLA-B\\*47", c_across(B_1:B_2)))]),
         B_48 = length(names(.)[24:25][which(grepl("HLA-B\\*48", c_across(B_1:B_2)))]),
         B_49 = length(names(.)[24:25][which(grepl("HLA-B\\*49", c_across(B_1:B_2)))]),
         B_50 = length(names(.)[24:25][which(grepl("HLA-B\\*50", c_across(B_1:B_2)))]),
         B_51 = length(names(.)[24:25][which(grepl("HLA-B\\*51", c_across(B_1:B_2)))]),
         B_52 = length(names(.)[24:25][which(grepl("HLA-B\\*52", c_across(B_1:B_2)))]),
         B_53 = length(names(.)[24:25][which(grepl("HLA-B\\*53", c_across(B_1:B_2)))]),
         B_54 = length(names(.)[24:25][which(grepl("HLA-B\\*54", c_across(B_1:B_2)))]),
         B_55 = length(names(.)[24:25][which(grepl("HLA-B\\*55", c_across(B_1:B_2)))]),
         B_56 = length(names(.)[24:25][which(grepl("HLA-B\\*56", c_across(B_1:B_2)))]),
         B_57 = length(names(.)[24:25][which(grepl("HLA-B\\*57", c_across(B_1:B_2)))]),
         B_58 = length(names(.)[24:25][which(grepl("HLA-B\\*58", c_across(B_1:B_2)))]),
         B_59 = length(names(.)[24:25][which(grepl("HLA-B\\*59", c_across(B_1:B_2)))]),
         B_67 = length(names(.)[24:25][which(grepl("HLA-B\\*67", c_across(B_1:B_2)))]),
         B_70 = length(names(.)[24:25][which(grepl("HLA-B\\*70", c_across(B_1:B_2)))]),
         B_73 = length(names(.)[24:25][which(grepl("HLA-B\\*73", c_across(B_1:B_2)))]),
         B_78 = length(names(.)[24:25][which(grepl("HLA-B\\*78", c_across(B_1:B_2)))]),
         B_81 = length(names(.)[24:25][which(grepl("HLA-B\\*81", c_across(B_1:B_2)))]),
         B_82 = length(names(.)[24:25][which(grepl("HLA-B\\*82", c_across(B_1:B_2)))]))

data <- data %>% rowwise %>% 
  mutate(C_01 = length(names(.)[26:27][which(grepl("HLA-C\\*01", c_across(C_1:C_2)))]),
         C_02 = length(names(.)[26:27][which(grepl("HLA-C\\*02", c_across(C_1:C_2)))]),
         C_03 = length(names(.)[26:27][which(grepl("HLA-C\\*03", c_across(C_1:C_2)))]),
         C_04 = length(names(.)[26:27][which(grepl("HLA-C\\*04", c_across(C_1:C_2)))]),
         C_05 = length(names(.)[26:27][which(grepl("HLA-C\\*05", c_across(C_1:C_2)))]),
         C_06 = length(names(.)[26:27][which(grepl("HLA-C\\*06", c_across(C_1:C_2)))]),
         C_07 = length(names(.)[26:27][which(grepl("HLA-C\\*07", c_across(C_1:C_2)))]),
         C_08 = length(names(.)[26:27][which(grepl("HLA-C\\*08", c_across(C_1:C_2)))]),
         C_12 = length(names(.)[26:27][which(grepl("HLA-C\\*12", c_across(C_1:C_2)))]),
         C_14 = length(names(.)[26:27][which(grepl("HLA-C\\*14", c_across(C_1:C_2)))]),
         C_15 = length(names(.)[26:27][which(grepl("HLA-C\\*15", c_across(C_1:C_2)))]),
         C_16 = length(names(.)[26:27][which(grepl("HLA-C\\*16", c_across(C_1:C_2)))]),
         C_17 = length(names(.)[26:27][which(grepl("HLA-C\\*17", c_across(C_1:C_2)))]),
         C_18 = length(names(.)[26:27][which(grepl("HLA-C\\*18", c_across(C_1:C_2)))]))

data <- data %>% rowwise %>% 
  mutate(DRB1_01 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*01", c_across(DRB1_1:DRB1_2)))]),
         DRB1_03 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*03", c_across(DRB1_1:DRB1_2)))]),
         DRB1_04 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*04", c_across(DRB1_1:DRB1_2)))]),
         DRB1_07 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*07", c_across(DRB1_1:DRB1_2)))]),
         DRB1_08 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*08", c_across(DRB1_1:DRB1_2)))]),
         DRB1_09 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*09", c_across(DRB1_1:DRB1_2)))]),
         DRB1_10 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*10", c_across(DRB1_1:DRB1_2)))]),
         DRB1_11 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*11", c_across(DRB1_1:DRB1_2)))]),
         DRB1_12 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*12", c_across(DRB1_1:DRB1_2)))]),
         DRB1_13 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*13", c_across(DRB1_1:DRB1_2)))]),
         DRB1_14 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*14", c_across(DRB1_1:DRB1_2)))]),
         DRB1_15 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*15", c_across(DRB1_1:DRB1_2)))]),
         DRB1_16 = length(names(.)[28:29][which(grepl("HLA-DRB1\\*16", c_across(DRB1_1:DRB1_2)))]))

data <- data %>% rowwise %>% 
  mutate(DQA1_01 = length(names(.)[30:31][which(grepl("HLA-DQA1\\*01", c_across(DQA1_1:DQA1_2)))]),
         DQA1_02 = length(names(.)[30:31][which(grepl("HLA-DQA1\\*02", c_across(DQA1_1:DQA1_2)))]),
         DQA1_03 = length(names(.)[30:31][which(grepl("HLA-DQA1\\*03", c_across(DQA1_1:DQA1_2)))]),
         DQA1_04 = length(names(.)[30:31][which(grepl("HLA-DQA1\\*04", c_across(DQA1_1:DQA1_2)))]),
         DQA1_05 = length(names(.)[30:31][which(grepl("HLA-DQA1\\*05", c_across(DQA1_1:DQA1_2)))]),
         DQA1_06 = length(names(.)[30:31][which(grepl("HLA-DQA1\\*06", c_across(DQA1_1:DQA1_2)))]))

data <- data %>% rowwise %>% 
  mutate(DQB1_02 = length(names(.)[32:33][which(grepl("HLA-DQB1\\*02", c_across(DQB1_1:DQB1_2)))]),
         DQB1_03 = length(names(.)[32:33][which(grepl("HLA-DQB1\\*03", c_across(DQB1_1:DQB1_2)))]),
         DQB1_04 = length(names(.)[32:33][which(grepl("HLA-DQB1\\*04", c_across(DQB1_1:DQB1_2)))]),
         DQB1_05 = length(names(.)[32:33][which(grepl("HLA-DQB1\\*05", c_across(DQB1_1:DQB1_2)))]),
         DQB1_06 = length(names(.)[32:33][which(grepl("HLA-DQB1\\*06", c_across(DQB1_1:DQB1_2)))]))

saveRDS(data, file="All_wb_exome.rds")


## HLA score calculations
data <- readRDS("All_wb_exome.rds")
data <- data %>% dplyr::select(ID:PC20,covid19,A_01:DQB1_06)

demo <- fread("ukb27449_20688_household.tab.gz")
demo <- demo %>% dplyr::select(f.eid, f.31.0.0,f.21022.0.0)

data <- data %>% inner_join(demo, by=c("ID"="f.eid"))

## calculate frequencies
Alleles <- colnames(data)[23:118]
out <- data.frame(matrix(0, length(Alleles), 1))
colnames(out) <- "Freq"
for(i in seq(1, length(Alleles))){
  out$Freq[i] <- sum(data[,Alleles[i]])/(dim(data)[1]*2)
}

## impose dosage 2 to 1
data <- data %>% mutate_at(.vars=vars(A_01:DQB1_06), .funs=funs(ifelse(.>=2, 1, .)))
data <- data %>% mutate(covid19 = ifelse(ID %in% covid_pos$eid, 1, 0))

## matrix of N x M (sample x HLA allele {0 or 1})
allele_mat <- as.matrix(data[,23:118])

## score calculations
score <- allele_mat %*% out$Freq

data1 <- data %>% dplyr::select(ID:covid19,f.31.0.0,f.21022.0.0)
data1$score <- score

quantile(data1$score, probs = seq(.1, .9, by = .1))
#10th centile 1.34, 90th centile 2.38

#defining Rare and high scores
data1 <- data1 %>% mutate(RareScore = ifelse(score <= 1.342034, 1, 0))
data1 <- data1 %>% mutate(HighScore = ifelse(score > 2.388517, 1, 0))

LM <- glm(covid19 ~ RareScore + f.31.0.0 + f.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=data1, family="binomial")
summary(LM)

LM <- glm(covid19 ~ HighScore + f.31.0.0 + f.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=data1, family="binomial")
summary(LM)

LM <- glm(covid19 ~ score + f.31.0.0 + f.21022.0.0 + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=data1, family="binomial")
summary(LM)
