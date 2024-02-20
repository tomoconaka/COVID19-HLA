# COVID19-HLA

## Analysis steps
### 1. define pairs according to the following criteria. `02.pairs_step1.R`.


Participants were assigned to couple pairs on the basis of a shared household identifier. 
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

