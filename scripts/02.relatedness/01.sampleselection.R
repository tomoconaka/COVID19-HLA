setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

couples1 <- readRDS("couples_step1.rds")
wb <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")

couples1 <- couples1 %>% dplyr::filter(f.eid %in% wb$V1)
hla_all <- fread("/project/richards/restricted/ukb-general/scratch/genetic-data/genome/imputed.v3/hla/ukb_hla_v2.txt.gz")
sample <- fread("/project/richards/restricted/ukb-27449/scratch/old-storage/full_release/v3/genotyped/w27449_20200204/ukb27449_cal_chr1_v2_20200204.fixCol6.fam", header=F)

hla_all <- bind_cols(sample$V1, hla_all)
colnames(hla_all)[1] <- "ID"

couples1 <- couples1 %>% inner_join(hla_all, by=c("f.eid"="ID"))

couples1 <- couples1 %>% group_by(uniqueCode) %>%
  dplyr::filter(length(f.eid) == 2) %>% ungroup()

couples1 <- couples1 %>% arrange(uniqueCode)

write.table(cbind(couples1$f.eid, couples1$f.eid, couples1$uniqueCode), file = "realcouples/WBcouples.sample", sep="\t", quote=F, col.names = F, row.names = F)

dat_female_index <- sample(seq(1,dim(couples1)[1]/2), size = dim(couples1)[1]/2,  replace = FALSE, prob = NULL)
dat_male_index <- sample(seq(1,dim(couples1)[1]/2), size = dim(couples1)[1]/2,  replace = FALSE, prob = NULL)
dat_female <- couples1 %>% dplyr::filter(f.31.0.0 == 0) %>% mutate(group = paste0("pair",dat_female_index))
dat_male <- couples1 %>% dplyr::filter(f.31.0.0 == 1) %>% mutate(group = paste0("pair",dat_male_index))
dat <- bind_rows(dat_female, dat_male)

dat <- dat %>% arrange(group)

dat <- dat %>% group_by(group) %>%
  dplyr::filter(length(f.eid) == 2) %>% ungroup()

write.table(cbind(dat$f.eid, dat$f.eid, dat$group), file = "randompairs/WBcouples.sample", sep="\t", quote=F, col.names = F, row.names = F)






