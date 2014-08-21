library(minerva)
## library(parallel)
myalpha <- 0.6
myc <- 15
data(Spellman)
Spellman <- as.matrix(Spellman)

## Spellman[,2] <- 0
mydata <- Spellman[,1:10]

res <- mine(Spellman,master=1,n.cores=5,alpha=myalpha,C=myc)

mydata[,2] <- 0
res <- mine(mydata,master=2,n.cores=6,alpha=myalpha,C=myc)



aval.cores <- detectCores()
if (aval.cores > 1){
  cat("Multicore: On this machinhe you have ",aval.cores," computational cores.\n")
}
if (aval.cores > 2){
  cat("We suggest to exploit minerva parallel computing possibilities with ",aval.cores-1," cores.\n(where possible set 'n.cores = 3')\n")
cat("Test ok!!!\n")
}


library(nettools)
for (i in 1:10) {
  idx <- sample(ncol(Spellman),1000)
  aa <- mat2adj(Spellman[,idx],method="ARACNE")
  write.table(aa,paste("test_data_",i,".csv",sep=""),
              col.names=NA,row.nam=TRUE, sep="\t",
              quote=FALSE)
}
