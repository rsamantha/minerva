library(minerva)
## library(parallel)
myalpha <- 0.6
myc <- 15
data(Spellman)
Spellman <- as.matrix(Spellman)

Spellman

res <- mine(Spellman,master=2,n.cores=2,alpha=myalpha,C=myc)


res <- mine(mydata,n.cores=6,alpha=myalpha,C=myc)



aval.cores <- detectCores()
if (aval.cores > 1){
  cat("Multicore: On this machinhe you have ",aval.cores," computational cores.\n")
}
if (aval.cores > 2){
  cat("We suggest to exploit minerva parallel computing possibilities with ",aval.cores-1," cores.\n(where possible set 'n.cores = 3')\n")
cat("Test ok!!!\n")
}
