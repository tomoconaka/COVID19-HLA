setwd("/home/richards/tomoko.nakanishi/scratch/09.COVID19/01.UKBB/03.HLA")

library(stringr)
couples <- readRDS(file="couples_step4_whichfirst.rds")

exome <- readRDS("couples_step4_WB_withHLA.rds")
# wb
pc <- fread("/home/richards/tomoko.nakanishi/scratch/DATA/UKB/EMC-White-British.20pc.txt")
colnames(pc) <- c("ID", paste0("PC",1:20))

couples <- couples %>% inner_join(pc, by=c("f.eid" = "ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

hla_all <- fread("/project/richards/restricted/ukb-general/scratch/genetic-data/genome/imputed.v3/hla/ukb_hla_v2.txt.gz")
sample <- fread("/project/richards/restricted/ukb-27449/scratch/old-storage/full_release/v3/genotyped/w27449_20200204/ukb27449_cal_chr1_v2_20200204.fixCol6.fam", header=F)

hla_all <- bind_cols(sample$V1, hla_all)
colnames(hla_all)[1] <- "ID"

couples <- couples %>% inner_join(hla_all, by=c("f.eid"="ID"))
couples <- couples %>% group_by(uniqueCode) %>%
  filter(length(f.eid) == 2) %>% ungroup()

couples <- data.frame(couples)
couples <- couples %>% drop_na(A_101:DQA1_601)
couples <- couples %>% rowwise() %>% mutate(A_01 = round(sum(c_across(A_101:A_103))),
                                            A_02 = round(sum(c_across(A_201:A_264))),
                                            A_03 = round(sum(c_across(A_301:A_302))),
                                            A_11 = round(sum(c_across(A_1101:A_1103))),
                                            A_23 = round(sum(c_across(A_2301))),
                                            A_24 = round(sum(c_across(A_2402:A_2410))),
                                            A_25 = round(sum(c_across(A_2501))),
                                            A_26 = round(sum(c_across(A_2601:A_2608))),
                                            A_29 = round(sum(c_across(A_2901:A_2902))),
                                            A_30 = round(sum(c_across(A_3001:A_3004))),
                                            A_31 = round(sum(c_across(A_3101:A_3104))),
                                            A_32 = round(sum(c_across(A_3201))),
                                            A_33 = round(sum(c_across(A_3301:A_3303))),
                                            A_34 = round(sum(c_across(A_3401:A_3402))),
                                            A_36 = round(sum(c_across(A_3601))),
                                            A_40 = round(sum(c_across(A_4027))),
                                            A_43 = round(sum(c_across(A_4386))),
                                            A_66 = round(sum(c_across(A_6601:A_6602))),
                                            A_68 = round(sum(c_across(A_6801:A_6817))),
                                            A_69 = round(sum(c_across(A_6901))),
                                            A_74 = round(sum(c_across(A_7401:A_7403))),
                                            A_80 = round(sum(c_across(A_8001))))

couples <- couples %>% rowwise() %>% mutate(B_07 = round(sum(c_across(B_702:B_705))),
                                            B_08 = round(sum(c_across(B_801))),
                                            B_13 = round(sum(c_across(B_1301:B_1302))),
                                            B_14 = round(sum(c_across(B_1401:B_1402))),
                                            B_15 = round(sum(c_across(B_1501:B_1537))),
                                            B_18 = round(sum(c_across(B_1801:B_1804))),
                                            B_27 = round(sum(c_across(B_2702:B_2707))),
                                            B_35 = round(sum(c_across(B_3501:B_3543))),
                                            B_37 = round(sum(c_across(B_3701))),
                                            B_38 = round(sum(c_across(B_3801:B_3802))),
                                            B_39 = round(sum(c_across(B_3901:B_3924))),
                                            B_40 = round(sum(c_across(B_4001:B_4012))),
                                            B_41 = round(sum(c_across(B_4101:B_4104))),
                                            B_42 = round(sum(c_across(B_4201))),
                                            B_44 = round(sum(c_across(B_4402:B_4407))),
                                            B_45 = round(sum(c_across(B_4501))),
                                            B_46 = round(sum(c_across(B_4601))),
                                            B_47 = round(sum(c_across(B_4701))),
                                            B_48 = round(sum(c_across(B_4801:B_4803))),
                                            B_49 = round(sum(c_across(B_4901))),
                                            B_50 = round(sum(c_across(B_5001:B_5002))),
                                            B_51 = round(sum(c_across(B_5101:B_5123))),
                                            B_52 = round(sum(c_across(B_5201))),
                                            B_53 = round(sum(c_across(B_5301))),
                                            B_54 = round(sum(c_across(B_5401))),
                                            B_55 = round(sum(c_across(B_5501:B_5502))),
                                            B_56 = round(sum(c_across(B_5601:B_5604))),
                                            B_57 = round(sum(c_across(B_5701:B_5703))),
                                            B_58 = round(sum(c_across(B_5801:B_5802))),
                                            B_59 = round(sum(c_across(B_5901))),
                                            B_67 = round(sum(c_across(B_6701))),
                                            B_70 = round(sum(c_across(B_7020))),
                                            B_73 = round(sum(c_across(B_7301))),
                                            B_78 = round(sum(c_across(B_7801))),
                                            B_81 = round(sum(c_across(B_8101))),
                                            B_82 = round(sum(c_across(B_8201:B_8202))))
                                            

couples <- couples %>% rowwise() %>% mutate(C_01 = round(sum(c_across(C_102))),
                                            C_02 = round(sum(c_across(C_202:C_210))),
                                            C_03 = round(sum(c_across(C_302:C_306))),
                                            C_04 = round(sum(c_across(C_401:C_407))),
                                            C_05 = round(sum(c_across(C_501))),
                                            C_06 = round(sum(c_across(C_602))),
                                            C_07 = round(sum(c_across(C_701:C_726))),
                                            C_08 = round(sum(c_across(C_801:C_804))),
                                            C_12 = round(sum(c_across(C_1202:C_1203))),
                                            C_14 = round(sum(c_across(C_1402:C_1403))),
                                            C_15 = round(sum(c_across(C_1502:C_1505))),
                                            C_16 = round(sum(c_across(C_1601:C_1604))),
                                            C_17 = round(sum(c_across(C_1701))),
                                            C_18 = round(sum(c_across(C_1801))))

couples <- couples %>% rowwise() %>% mutate(DRB1_01 = round(sum(c_across(DRB1_101:DRB1_103))),
                                            DRB1_03 = round(sum(c_across(DRB1_301:DRB1_302))),
                                            DRB1_04 = round(sum(c_across(DRB1_401:DRB1_410))),
                                            DRB1_07 = round(sum(c_across(DRB1_701))),
                                            DRB1_08 = round(sum(c_across(DRB1_801:DRB1_811))),
                                            DRB1_09 = round(sum(c_across(DRB1_901))),
                                            DRB1_10 = round(sum(c_across(DRB1_1001))),
                                            DRB1_11 = round(sum(c_across(DRB1_1101:DRB1_1143))),
                                            DRB1_12 = round(sum(c_across(DRB1_1201:DRB1_1202))),
                                            DRB1_13 = round(sum(c_across(DRB1_1301:DRB1_1321))),
                                            DRB1_14 = round(sum(c_across(DRB1_1401:DRB1_1444))),
                                            DRB1_15 = round(sum(c_across(DRB1_1501:DRB1_1503))),
                                            DRB1_16 = round(sum(c_across(DRB1_1601:DRB1_1602))))

couples <- couples %>% rowwise() %>% mutate(DQA1_01 = round(sum(c_across(DQA1_101:DQA1_104))),
                                            DQA1_02 = round(sum(c_across(DQA1_201))),
                                            DQA1_03 = round(sum(c_across(DQA1_301:DQA1_303))),
                                            DQA1_04 = round(sum(c_across(DQA1_401))),
                                            DQA1_05 = round(sum(c_across(DQA1_501:DQA1_509))),
                                            DQA1_06 = round(sum(c_across(DQA1_601))))

couples <- couples %>% rowwise() %>% mutate(DQB1_02 = round(sum(c_across(DQB1_201:DQB1_202))),
                                            DQB1_03 = round(sum(c_across(DQB1_301:DQB1_304))),
                                            DQB1_04 = round(sum(c_across(DQB1_401:DQB1_402))),
                                            DQB1_05 = round(sum(c_across(DQB1_501:DQB1_504))),
                                            DQB1_06 = round(sum(c_across(DQB1_601:DQB1_609))))

impute <- couples %>% dplyr::select(f.eid,A_01:DQB1_06)

exome <- exome %>% rowwise %>% 
  mutate(A_01 = length(names(.)[34:35][which(grepl("HLA-A\\*01", c_across(A_1:A_2)))]),
         A_02 = length(names(.)[34:35][which(grepl("HLA-A\\*02", c_across(A_1:A_2)))]),
         A_03 = length(names(.)[34:35][which(grepl("HLA-A\\*03", c_across(A_1:A_2)))]),
         A_11 = length(names(.)[34:35][which(grepl("HLA-A\\*11", c_across(A_1:A_2)))]),
         A_23 = length(names(.)[34:35][which(grepl("HLA-A\\*23", c_across(A_1:A_2)))]),
         A_24 = length(names(.)[34:35][which(grepl("HLA-A\\*24", c_across(A_1:A_2)))]),
         A_25 = length(names(.)[34:35][which(grepl("HLA-A\\*25", c_across(A_1:A_2)))]),
         A_26 = length(names(.)[34:35][which(grepl("HLA-A\\*26", c_across(A_1:A_2)))]),
         A_29 = length(names(.)[34:35][which(grepl("HLA-A\\*29", c_across(A_1:A_2)))]),
         A_30 = length(names(.)[34:35][which(grepl("HLA-A\\*30", c_across(A_1:A_2)))]),
         A_31 = length(names(.)[34:35][which(grepl("HLA-A\\*31", c_across(A_1:A_2)))]),
         A_32 = length(names(.)[34:35][which(grepl("HLA-A\\*32", c_across(A_1:A_2)))]),
         A_33 = length(names(.)[34:35][which(grepl("HLA-A\\*33", c_across(A_1:A_2)))]),
         A_34 = length(names(.)[34:35][which(grepl("HLA-A\\*34", c_across(A_1:A_2)))]),
         A_36 = length(names(.)[34:35][which(grepl("HLA-A\\*36", c_across(A_1:A_2)))]),
         A_40 = length(names(.)[34:35][which(grepl("HLA-A\\*40", c_across(A_1:A_2)))]),
         A_43 = length(names(.)[34:35][which(grepl("HLA-A\\*43", c_across(A_1:A_2)))]),
         A_66 = length(names(.)[34:35][which(grepl("HLA-A\\*66", c_across(A_1:A_2)))]),
         A_68 = length(names(.)[34:35][which(grepl("HLA-A\\*68", c_across(A_1:A_2)))]),
         A_69 = length(names(.)[34:35][which(grepl("HLA-A\\*69", c_across(A_1:A_2)))]),
         A_74 = length(names(.)[34:35][which(grepl("HLA-A\\*74", c_across(A_1:A_2)))]),
         A_80 = length(names(.)[34:35][which(grepl("HLA-A\\*80", c_across(A_1:A_2)))])
         )

exome <- exome %>% rowwise %>% 
  mutate(B_07 = length(names(.)[36:37][which(grepl("HLA-B\\*07", c_across(B_1:B_2)))]),
         B_08 = length(names(.)[36:37][which(grepl("HLA-B\\*08", c_across(B_1:B_2)))]),
         B_13 = length(names(.)[36:37][which(grepl("HLA-B\\*13", c_across(B_1:B_2)))]),
         B_14 = length(names(.)[36:37][which(grepl("HLA-B\\*14", c_across(B_1:B_2)))]),
         B_15 = length(names(.)[36:37][which(grepl("HLA-B\\*15", c_across(B_1:B_2)))]),
         B_18 = length(names(.)[36:37][which(grepl("HLA-B\\*18", c_across(B_1:B_2)))]),
         B_27 = length(names(.)[36:37][which(grepl("HLA-B\\*27", c_across(B_1:B_2)))]),
         B_35 = length(names(.)[36:37][which(grepl("HLA-B\\*35", c_across(B_1:B_2)))]),
         B_37 = length(names(.)[36:37][which(grepl("HLA-B\\*37", c_across(B_1:B_2)))]),
         B_38 = length(names(.)[36:37][which(grepl("HLA-B\\*38", c_across(B_1:B_2)))]),
         B_39 = length(names(.)[36:37][which(grepl("HLA-B\\*39", c_across(B_1:B_2)))]),
         B_40 = length(names(.)[36:37][which(grepl("HLA-B\\*40", c_across(B_1:B_2)))]),
         B_41 = length(names(.)[36:37][which(grepl("HLA-B\\*41", c_across(B_1:B_2)))]),
         B_42 = length(names(.)[36:37][which(grepl("HLA-B\\*42", c_across(B_1:B_2)))]),
         B_44 = length(names(.)[36:37][which(grepl("HLA-B\\*44", c_across(B_1:B_2)))]),
         B_45 = length(names(.)[36:37][which(grepl("HLA-B\\*45", c_across(B_1:B_2)))]),
         B_46 = length(names(.)[36:37][which(grepl("HLA-B\\*46", c_across(B_1:B_2)))]),
         B_47 = length(names(.)[36:37][which(grepl("HLA-B\\*47", c_across(B_1:B_2)))]),
         B_48 = length(names(.)[36:37][which(grepl("HLA-B\\*48", c_across(B_1:B_2)))]),
         B_49 = length(names(.)[36:37][which(grepl("HLA-B\\*49", c_across(B_1:B_2)))]),
         B_50 = length(names(.)[36:37][which(grepl("HLA-B\\*50", c_across(B_1:B_2)))]),
         B_51 = length(names(.)[36:37][which(grepl("HLA-B\\*51", c_across(B_1:B_2)))]),
         B_52 = length(names(.)[36:37][which(grepl("HLA-B\\*52", c_across(B_1:B_2)))]),
         B_53 = length(names(.)[36:37][which(grepl("HLA-B\\*53", c_across(B_1:B_2)))]),
         B_54 = length(names(.)[36:37][which(grepl("HLA-B\\*54", c_across(B_1:B_2)))]),
         B_55 = length(names(.)[36:37][which(grepl("HLA-B\\*55", c_across(B_1:B_2)))]),
         B_56 = length(names(.)[36:37][which(grepl("HLA-B\\*56", c_across(B_1:B_2)))]),
         B_57 = length(names(.)[36:37][which(grepl("HLA-B\\*57", c_across(B_1:B_2)))]),
         B_58 = length(names(.)[36:37][which(grepl("HLA-B\\*58", c_across(B_1:B_2)))]),
         B_59 = length(names(.)[36:37][which(grepl("HLA-B\\*59", c_across(B_1:B_2)))]),
         B_67 = length(names(.)[36:37][which(grepl("HLA-B\\*67", c_across(B_1:B_2)))]),
         B_70 = length(names(.)[36:37][which(grepl("HLA-B\\*70", c_across(B_1:B_2)))]),
         B_73 = length(names(.)[36:37][which(grepl("HLA-B\\*73", c_across(B_1:B_2)))]),
         B_78 = length(names(.)[36:37][which(grepl("HLA-B\\*78", c_across(B_1:B_2)))]),
         B_81 = length(names(.)[36:37][which(grepl("HLA-B\\*81", c_across(B_1:B_2)))]),
         B_82 = length(names(.)[36:37][which(grepl("HLA-B\\*82", c_across(B_1:B_2)))]))

exome <- exome %>% rowwise %>% 
  mutate(C_01 = length(names(.)[38:39][which(grepl("HLA-C\\*01", c_across(C_1:C_2)))]),
         C_02 = length(names(.)[38:39][which(grepl("HLA-C\\*02", c_across(C_1:C_2)))]),
         C_03 = length(names(.)[38:39][which(grepl("HLA-C\\*03", c_across(C_1:C_2)))]),
         C_04 = length(names(.)[38:39][which(grepl("HLA-C\\*04", c_across(C_1:C_2)))]),
         C_05 = length(names(.)[38:39][which(grepl("HLA-C\\*05", c_across(C_1:C_2)))]),
         C_06 = length(names(.)[38:39][which(grepl("HLA-C\\*06", c_across(C_1:C_2)))]),
         C_07 = length(names(.)[38:39][which(grepl("HLA-C\\*07", c_across(C_1:C_2)))]),
         C_08 = length(names(.)[38:39][which(grepl("HLA-C\\*08", c_across(C_1:C_2)))]),
         C_12 = length(names(.)[38:39][which(grepl("HLA-C\\*12", c_across(C_1:C_2)))]),
         C_14 = length(names(.)[38:39][which(grepl("HLA-C\\*14", c_across(C_1:C_2)))]),
         C_15 = length(names(.)[38:39][which(grepl("HLA-C\\*15", c_across(C_1:C_2)))]),
         C_16 = length(names(.)[38:39][which(grepl("HLA-C\\*16", c_across(C_1:C_2)))]),
         C_17 = length(names(.)[38:39][which(grepl("HLA-C\\*17", c_across(C_1:C_2)))]),
         C_18 = length(names(.)[38:39][which(grepl("HLA-C\\*18", c_across(C_1:C_2)))]))
         
exome <- exome %>% rowwise %>% 
  mutate(DRB1_01 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*01", c_across(DRB1_1:DRB1_2)))]),
         DRB1_03 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*03", c_across(DRB1_1:DRB1_2)))]),
         DRB1_04 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*04", c_across(DRB1_1:DRB1_2)))]),
         DRB1_07 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*07", c_across(DRB1_1:DRB1_2)))]),
         DRB1_08 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*08", c_across(DRB1_1:DRB1_2)))]),
         DRB1_09 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*09", c_across(DRB1_1:DRB1_2)))]),
         DRB1_10 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*10", c_across(DRB1_1:DRB1_2)))]),
         DRB1_11 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*11", c_across(DRB1_1:DRB1_2)))]),
         DRB1_12 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*12", c_across(DRB1_1:DRB1_2)))]),
         DRB1_13 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*13", c_across(DRB1_1:DRB1_2)))]),
         DRB1_14 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*14", c_across(DRB1_1:DRB1_2)))]),
         DRB1_15 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*15", c_across(DRB1_1:DRB1_2)))]),
         DRB1_16 = length(names(.)[40:41][which(grepl("HLA-DRB1\\*16", c_across(DRB1_1:DRB1_2)))]))

exome <- exome %>% rowwise %>% 
  mutate(DQA1_01 = length(names(.)[42:43][which(grepl("HLA-DQA1\\*01", c_across(DQA1_1:DQA1_2)))]),
         DQA1_02 = length(names(.)[42:43][which(grepl("HLA-DQA1\\*02", c_across(DQA1_1:DQA1_2)))]),
         DQA1_03 = length(names(.)[42:43][which(grepl("HLA-DQA1\\*03", c_across(DQA1_1:DQA1_2)))]),
         DQA1_04 = length(names(.)[42:43][which(grepl("HLA-DQA1\\*04", c_across(DQA1_1:DQA1_2)))]),
         DQA1_05 = length(names(.)[42:43][which(grepl("HLA-DQA1\\*05", c_across(DQA1_1:DQA1_2)))]),
         DQA1_06 = length(names(.)[42:43][which(grepl("HLA-DQA1\\*06", c_across(DQA1_1:DQA1_2)))]))

exome <- exome %>% rowwise %>% 
  mutate(DQB1_02 = length(names(.)[44:45][which(grepl("HLA-DQB1\\*02", c_across(DQB1_1:DQB1_2)))]),
         DQB1_03 = length(names(.)[44:45][which(grepl("HLA-DQB1\\*03", c_across(DQB1_1:DQB1_2)))]),
         DQB1_04 = length(names(.)[44:45][which(grepl("HLA-DQB1\\*04", c_across(DQB1_1:DQB1_2)))]),
         DQB1_05 = length(names(.)[44:45][which(grepl("HLA-DQB1\\*05", c_across(DQB1_1:DQB1_2)))]),
         DQB1_06 = length(names(.)[44:45][which(grepl("HLA-DQB1\\*06", c_across(DQB1_1:DQB1_2)))]))


data_exome <- exome %>% dplyr::select(f.eid,A_01:DQB1_06)

df <- inner_join(impute, data_exome, by="f.eid")
saveRDS(df, file="imputation_exome_data.rds")

A_alleles <- colnames(data_exome)[-1]
out <- data.frame(matrix(0, length(A_alleles), 7))
colnames(out) <- c("Allele","ImputeFreq", "ExomeFreq", "R", "R.LL","R.UL","R.pvalue")
for(i in seq(1, length(A_alleles))){
  out$Allele[i] <- A_alleles[i]
  out$ImputeFreq[i] <- sum(df[,paste0(A_alleles[i], ".x")])/(dim(df)[1]*2)
  out$ExomeFreq[i] <- sum(df[,paste0(A_alleles[i], ".y")])/(dim(df)[1]*2)
  out$R[i] <- cor.test(unlist(df[,paste0(A_alleles[i], ".x")]), unlist(df[,paste0(A_alleles[i], ".y")]))$estimate
  out[i,5:6] <- cor.test(unlist(df[,paste0(A_alleles[i], ".x")]), unlist(df[,paste0(A_alleles[i], ".y")]))$conf.int
  out$R.pvalue[i] <- cor.test(unlist(df[,paste0(A_alleles[i], ".x")]), unlist(df[,paste0(A_alleles[i], ".y")]))$p.value
}

write.table(out, file="imputation_exome_r2.tsv", quote=F, col.names = T, row.names = F, sep="\t")