library(minerva)
mydata <- read.table("../mictools/examples/datasaurus.txt", header=TRUE, row.names = 1, sep="\t", as.is=TRUE, stringsAsFactors = FALSE)
mydata <- t(mydata)
mmres <- mictools(mydata, nperm=1000)

obstic.df <- as.data.frame(mmres$obstic)
obstic.df$Var1 <- colnames(mydata)[obstic.df[,2]]
obstic.df$Var2 <- colnames(mydata)[obstic.df[,3]]

nullpython <- read.table("../mictools/examples/datasaurus_results/null_dist.txt", header=TRUE, stringsAsFactors = FALSE)
obspython <- read.table("../mictools/examples/datasaurus_results/obs.txt", header=TRUE, stringsAsFactors = FALSE)
obsdistpython <- read.table("../mictools/examples/datasaurus_results/obs_dist.txt", header=TRUE, stringsAsFactors = FALSE)
obspvalpython <- read.table("../mictools/examples/datasaurus_results/pval.txt", header=TRUE, stringsAsFactors = FALSE)

## Check observed distribution with original mictools version
any(sapply(obspython[,3] - obstic.df[,1], all.equal, current=0, tolerance = 1e-6)==FALSE)
any(sapply(obsdistpython$ObsCount - mmres$obsdist$Count, all.equal, current=0, tolerance = 1e-6)==FALSE)
any(sapply(obsdistpython$ObsCountCum - mmres$obsdist$CountCum, all.equal, current=0, tolerance = 1e-6)==FALSE)
any(sapply(obspvalpython$None - mmres$pval$pval, all.equal, current=0, tolerance = 1e-6)==FALSE)

p.adjust(mmres$pval$pval) <0.05
p.adjust(obspvalpython$None) < 0.05


obspvalpython$None[1:10]
mmres$pval$pval[1:10]
obspvalpython$None - mmres$pval$pval
