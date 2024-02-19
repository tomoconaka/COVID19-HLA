setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA/")

#HLA Nsnps = 4755
realcouples_HLA <- readRDS("realcouples/SNP_diversity_WBcouples_HLA.rds")
randompairs_HLA <- readRDS("randompairs/SNP_diversity_WBcouples_HLA.rds")

realcouples_HLA <- realcouples %>% mutate(R = (Qc - mean(randompairs$Qc))/(1 - mean(randompairs$Qc)))

#Other genome
realcouples_chr1 <- readRDS("realcouples/SNP_diversity_WBcouples_1.rds")# %>% mutate(Qc_1 = Qc*5523) %>% select(uniqueCode, Qc_1)
realcouples_chr2 <- readRDS("realcouples/SNP_diversity_WBcouples_2.rds")# %>% mutate(Qc_2 = Qc*5316) %>% select(uniqueCode, Qc_2)
realcouples_chr3 <- readRDS("realcouples/SNP_diversity_WBcouples_3.rds")# %>% mutate(Qc_3 = Qc*4616) %>% select(uniqueCode, Qc_3)
realcouples_chr4 <- readRDS("realcouples/SNP_diversity_WBcouples_4.rds")# %>% mutate(Qc_4 = Qc*4339) %>% select(uniqueCode, Qc_4)
realcouples_chr5 <- readRDS("realcouples/SNP_diversity_WBcouples_5.rds")# %>% mutate(Qc_5 = Qc*4103) %>% select(uniqueCode, Qc_5)
realcouples_chr6 <- readRDS("realcouples/SNP_diversity_WBcouples_6.rds")# %>% mutate(Qc_6 = Qc*3871) %>% select(uniqueCode, Qc_6)
realcouples_chr7 <- readRDS("realcouples/SNP_diversity_WBcouples_7.rds")# %>% mutate(Qc_7 = Qc*3763) %>% select(uniqueCode, Qc_7)
realcouples_chr8 <- readRDS("realcouples/SNP_diversity_WBcouples_8.rds")# %>% mutate(Qc_8 = Qc*3480) %>% select(uniqueCode, Qc_8)
realcouples_chr9 <- readRDS("realcouples/SNP_diversity_WBcouples_9.rds")# %>% mutate(Qc_9 = Qc*3250) %>% select(uniqueCode, Qc_9)
realcouples_chr10 <- readRDS("realcouples/SNP_diversity_WBcouples_10.rds")# %>% mutate(Qc_10 = Qc*3566) %>% select(uniqueCode, Qc_10)
realcouples_chr11 <- readRDS("realcouples/SNP_diversity_WBcouples_11.rds")# %>% mutate(Qc_11 = Qc*3276) %>% select(uniqueCode, Qc_11)
realcouples_chr12 <- readRDS("realcouples/SNP_diversity_WBcouples_12.rds")# %>% mutate(Qc_12 = Qc*3397) %>% select(uniqueCode, Qc_12)
realcouples_chr13 <- readRDS("realcouples/SNP_diversity_WBcouples_13.rds")# %>% mutate(Qc_13 = Qc*2595) %>% select(uniqueCode, Qc_13)
realcouples_chr14 <- readRDS("realcouples/SNP_diversity_WBcouples_14.rds")# %>% mutate(Qc_14 = Qc*2372) %>% select(uniqueCode, Qc_14)
realcouples_chr15 <- readRDS("realcouples/SNP_diversity_WBcouples_15.rds")# %>% mutate(Qc_15 = Qc*2375) %>% select(uniqueCode, Qc_15)
realcouples_chr16 <- readRDS("realcouples/SNP_diversity_WBcouples_16.rds")# %>% mutate(Qc_16 = Qc*2587) %>% select(uniqueCode, Qc_16)
realcouples_chr17 <- readRDS("realcouples/SNP_diversity_WBcouples_17.rds")# %>% mutate(Qc_17 = Qc*2436) %>% select(uniqueCode, Qc_17)
realcouples_chr18 <- readRDS("realcouples/SNP_diversity_WBcouples_18.rds")# %>% mutate(Qc_18 = Qc*2423) %>% select(uniqueCode, Qc_18)
realcouples_chr19 <- readRDS("realcouples/SNP_diversity_WBcouples_19.rds")# %>% mutate(Qc_19 = Qc*2021) %>% select(uniqueCode, Qc_19)
realcouples_chr20 <- readRDS("realcouples/SNP_diversity_WBcouples_20.rds")# %>% mutate(Qc_20 = Qc*2044) %>% select(uniqueCode, Qc_20)
realcouples_chr21 <- readRDS("realcouples/SNP_diversity_WBcouples_21.rds")# %>% mutate(Qc_21 = Qc*1179) %>% select(uniqueCode, Qc_21)
realcouples_chr22 <- readRDS("realcouples/SNP_diversity_WBcouples_22.rds")# %>% mutate(Qc_22 = Qc*1252) %>% select(uniqueCode, Qc_22)

realcouples_all_genome <- chr1 %>% inner_join(chr2, by="uniqueCode") %>% inner_join(chr3, by="uniqueCode") %>% inner_join(chr4, by="uniqueCode") %>%
  inner_join(chr5, by="uniqueCode") %>% inner_join(chr6, by="uniqueCode") %>% inner_join(chr7, by="uniqueCode") %>% inner_join(chr8, by="uniqueCode") %>% 
  inner_join(chr9, by="uniqueCode") %>% inner_join(chr10, by="uniqueCode") %>% inner_join(chr11, by="uniqueCode") %>% inner_join(chr12, by="uniqueCode") %>% 
  inner_join(chr13, by="uniqueCode") %>% inner_join(chr14, by="uniqueCode") %>% inner_join(chr15, by="uniqueCode") %>% inner_join(chr16, by="uniqueCode") %>% 
  inner_join(chr17, by="uniqueCode") %>% inner_join(chr18, by="uniqueCode") %>% inner_join(chr19, by="uniqueCode") %>% inner_join(chr20, by="uniqueCode") %>% 
  inner_join(chr21, by="uniqueCode") %>% inner_join(chr22, by="uniqueCode")

realcouples_all_genome <- realcouples_all_genome %>% mutate(Qc = realcouples_all_genome %>% select(Qc_1:Qc_22) %>% rowSums)

realcouples_all_genome <- realcouples_all_genome %>% mutate(Qc = Qc/69784)


#Random
randompairs_chr1 <- readRDS("randompairs/SNP_diversity_WBcouples_1.rds")# %>% mutate(Qc_1 = Qc*5523) %>% select(uniqueCode, Qc_1)
randompairs_chr2 <- readRDS("randompairs/SNP_diversity_WBcouples_2.rds")# %>% mutate(Qc_2 = Qc*5316) %>% select(uniqueCode, Qc_2)
randompairs_chr3 <- readRDS("randompairs/SNP_diversity_WBcouples_3.rds")# %>% mutate(Qc_3 = Qc*4616) %>% select(uniqueCode, Qc_3)
randompairs_chr4 <- readRDS("randompairs/SNP_diversity_WBcouples_4.rds")# %>% mutate(Qc_4 = Qc*4339) %>% select(uniqueCode, Qc_4)
randompairs_chr5 <- readRDS("randompairs/SNP_diversity_WBcouples_5.rds")# %>% mutate(Qc_5 = Qc*4103) %>% select(uniqueCode, Qc_5)
randompairs_chr6 <- readRDS("randompairs/SNP_diversity_WBcouples_6.rds")# %>% mutate(Qc_6 = Qc*3871) %>% select(uniqueCode, Qc_6)
randompairs_chr7 <- readRDS("randompairs/SNP_diversity_WBcouples_7.rds")# %>% mutate(Qc_7 = Qc*3763) %>% select(uniqueCode, Qc_7)
randompairs_chr8 <- readRDS("randompairs/SNP_diversity_WBcouples_8.rds")# %>% mutate(Qc_8 = Qc*3480) %>% select(uniqueCode, Qc_8)
randompairs_chr9 <- readRDS("randompairs/SNP_diversity_WBcouples_9.rds")# %>% mutate(Qc_9 = Qc*3250) %>% select(uniqueCode, Qc_9)
randompairs_chr10 <- readRDS("randompairs/SNP_diversity_WBcouples_10.rds")# %>% mutate(Qc_10 = Qc*3566) %>% select(uniqueCode, Qc_10)
randompairs_chr11 <- readRDS("randompairs/SNP_diversity_WBcouples_11.rds")# %>% mutate(Qc_11 = Qc*3276) %>% select(uniqueCode, Qc_11)
randompairs_chr12 <- readRDS("randompairs/SNP_diversity_WBcouples_12.rds")# %>% mutate(Qc_12 = Qc*3397) %>% select(uniqueCode, Qc_12)
randompairs_chr13 <- readRDS("randompairs/SNP_diversity_WBcouples_13.rds")# %>% mutate(Qc_13 = Qc*2595) %>% select(uniqueCode, Qc_13)
randompairs_chr14 <- readRDS("randompairs/SNP_diversity_WBcouples_14.rds")# %>% mutate(Qc_14 = Qc*2372) %>% select(uniqueCode, Qc_14)
randompairs_chr15 <- readRDS("randompairs/SNP_diversity_WBcouples_15.rds")# %>% mutate(Qc_15 = Qc*2375) %>% select(uniqueCode, Qc_15)
randompairs_chr16 <- readRDS("randompairs/SNP_diversity_WBcouples_16.rds")# %>% mutate(Qc_16 = Qc*2587) %>% select(uniqueCode, Qc_16)
randompairs_chr17 <- readRDS("randompairs/SNP_diversity_WBcouples_17.rds")# %>% mutate(Qc_17 = Qc*2436) %>% select(uniqueCode, Qc_17)
randompairs_chr18 <- readRDS("randompairs/SNP_diversity_WBcouples_18.rds")# %>% mutate(Qc_18 = Qc*2423) %>% select(uniqueCode, Qc_18)
randompairs_chr19 <- readRDS("randompairs/SNP_diversity_WBcouples_19.rds")# %>% mutate(Qc_19 = Qc*2021) %>% select(uniqueCode, Qc_19)
randompairs_chr20 <- readRDS("randompairs/SNP_diversity_WBcouples_20.rds")# %>% mutate(Qc_20 = Qc*2044) %>% select(uniqueCode, Qc_20)
randompairs_chr21 <- readRDS("randompairs/SNP_diversity_WBcouples_21.rds")# %>% mutate(Qc_21 = Qc*1179) %>% select(uniqueCode, Qc_21)
randompairs_chr22 <- readRDS("randompairs/SNP_diversity_WBcouples_22.rds")# %>% mutate(Qc_22 = Qc*1252) %>% select(uniqueCode, Qc_22)

randompairs_all_genome <- chr1 %>% inner_join(chr2, by="uniqueCode") %>% inner_join(chr3, by="uniqueCode") %>% inner_join(chr4, by="uniqueCode") %>%
  inner_join(chr5, by="uniqueCode") %>% inner_join(chr6, by="uniqueCode") %>% inner_join(chr7, by="uniqueCode") %>% inner_join(chr8, by="uniqueCode") %>% 
  inner_join(chr9, by="uniqueCode") %>% inner_join(chr10, by="uniqueCode") %>% inner_join(chr11, by="uniqueCode") %>% inner_join(chr12, by="uniqueCode") %>% 
  inner_join(chr13, by="uniqueCode") %>% inner_join(chr14, by="uniqueCode") %>% inner_join(chr15, by="uniqueCode") %>% inner_join(chr16, by="uniqueCode") %>% 
  inner_join(chr17, by="uniqueCode") %>% inner_join(chr18, by="uniqueCode") %>% inner_join(chr19, by="uniqueCode") %>% inner_join(chr20, by="uniqueCode") %>% 
  inner_join(chr21, by="uniqueCode") %>% inner_join(chr22, by="uniqueCode")

randompairs_all_genome <- randompairs_all_genome %>% mutate(Qc = randompairs_all_genome %>% select(Qc_1:Qc_22) %>% rowSums)

randompairs_all_genome <- randompairs_all_genome %>% mutate(Qc = Qc/69784)

realcouples_all_genome <- realcouples_all_genome %>% select(uniqueCode, Qc) %>% 
  mutate(R = (Qc - mean(randompairs_all_genome$Qc))/(1 - mean(randompairs_all_genome$Qc)))

t.test(realcouples_all_genome$Qc, randompairs_all_genome$Qc)


##
out <- data.frame(matrix(0, 23, 5))
colnames(out) <- c("genome", "Nsnps", "mean_real", "mean_couples", "p.value")

res <- t.test(realcouples_HLA$Qc, randompairs_HLA$Qc)
out$genome[1] <- "HLA"
out$Nsnps[1] <- unique(realcouples_HLA$Nsnp)
out[1,3:5] <- c(res$estimate[1], res$estimate[2], res$p.value)

res <- t.test(realcouples_chr1$Qc, randompairs_chr1$Qc)
out$genome[2] <- "chr1"
out$Nsnps[2] <- 4755
out[2,3:5] <- c(res$estimate[1], res$estimate[2], res$p.value)

res <- t.test(realcouples_chr2$Qc, randompairs_chr2$Qc)
out$genome[3] <- "chr2"
out$Nsnps[3] <- 4755
out[3,3:5] <- c(res$estimate[1], res$estimate[2], res$p.value)

res <- t.test(realcouples_chr3$Qc, randompairs_chr3$Qc)
out$genome[4] <- "chr3"
out$Nsnps[4] <- 4755
out[4,3:5] <- c(res$estimate[1], res$estimate[2], res$p.value)

