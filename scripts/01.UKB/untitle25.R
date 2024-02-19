
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2519788/

## theoretical 
# p, allele frequency
p <- 0.2
calculateQc <- function(p){
  Qc <- 1*p^4 + 0.5*p^2*2*p*(1-p) + 0*(1-p)^2*p^2 +  
    0.5*p^2*2*p*(1-p) + 0.5*4*p^2*(1-p)^2 + 0.5*2*p*(1-p)^3 + 
    0*(1-p)^2*p^2 + 0.5*2*p*(1-p)^3 + 1*(1-p)^4
}



