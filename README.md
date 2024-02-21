# COVID19-HLA

## Analysis steps
### 1. define pairs according to the following criteria. `02.pairs_step1.R`.

In UK Biobank, participants were assigned to couple pairs on the basis of a shared household identifier. 
* Make a unique “couple code” representing “a couple” using following key variables.
1) Home location at assessment (rounded east/north co-ordinate)  i.e 367000 and 437000
2) Assessment centre
3) Type of accommodation lived in
4) Own or rent accommodation lived in
5) Length of time at current address
6) Number in household
7) Number of vehicles in household
only two

* who reported living with a husband, wife, or partner were selected.
* If more than two individuals shared identical information across all variables, these individuals were excluded from analysis. 

These criteria is based on [ref](https://www.nature.com/articles/s41467-019-12424-x#Sec12) 

### 2. Identifying pairs with at least one infected with COVID-19 before Jan 2021. `03.pairs_step23.R`

* COVID-19 PCR tests date before Dec 31 2020 were considered in the analysis. (To include only those who were unvaccinated)
* If both individuals in the pairs were infected, `interval` was calculated based on the date of COVID-19 PCR tests.
* Any pairs with `interval` ≥1 and ≤7 were considered as **secondary attack**.
![662e0dd7-6604-4ca9-b5aa-1218a0921818](https://github.com/tomoconaka/COVID19-HLA/assets/48235580/fd4d5062-7994-4916-b1b7-72cfca52b773)
* pairs with at least one COVID-19 infection and both of individuals are alive before COVID-19 pandemic.


### 3. Identifying the index person in the pairs. `04.pairs_step4_whichfirst.R`

### 4. Ancestry assignment and HLA typing based on WES data. `05.HLA_ancestry.R`

* details of ancestry assignment strategy are described in [Morris et al, Nat Genet, 2018](https://www.nature.com/articles/s41588-018-0302-x)
* details of HLA typing based on WES are described in [Butler-Laporte, et al, Commun Biol 2023](https://www.nature.com/articles/s42003-023-05496-5)
* Use only pairs with the same continental ancestry. (i.e. not include EUR and AFR couples.) 

### 5. association with COVID-19 infection and HLA incompatibility `06.HLA_onefield.R`

HLA incompatibility is defined for each HLA genotype in class 1 (HLA-A, HLA-B, HLA-C) and class 2 (HLA-DRB1, HLA-DQA1, HLA-DQB1) using one-filed.
For each HLA genotypes, 
HLA compatibility (`HLAmatch`) is defined if both of two HLA alleles in the index person are included in the second person.

i.e) `HLAAmatch == TRUE` if the index person has `A*02` and `A*02` and the second person has `A*02` and `A*25`.
     `HLAAmatch == FALSE` if the index person has `A*02` and `A*25` and the second person has `A*02` and `A*02`.

Logistic regression was performed with the following formula.

Secondary attack ~ HLAmatch + age_diff + age_sum + sex_index + sex_second * PC_diff{1..10}

`age_diff`: age difference of the pairs
`age_sum`: age sum of the pairs
`sex_index`: sex of the index person
`sex_second`: sex of the second person
`PC_diff{1..10}`: absolute value of the difference of genetic PCs

Send us the summary statistics of `beta` , `se` , `pval` , `Ncase` , `Ncontrol`. `Ncase` is the number of pairs with co-infection, `Ncontrols` is the number of pairs without co-infection.

### 6. correlation analyses between each HLA incompatibility

To run meta-analyses, we need correlation information between each HLA incompatibility.

```
tmp_cor <- tmp %>% dplyr::select(colnames(couples)[grepl("match", colnames(couples))])
library(corrplot)
M <- rcorr(as.matrix(tmp_cor))
```
`tmp_cor` : N(number of spausal pairs) x M (number of HLA genes, 6 or 5) matrix, representing whether their HLAs are matched or not (1 or 0)
send us `M` data.



