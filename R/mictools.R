# mictools_null_r <- function(x){
#   ix <- sample(ncol(x), 2)
#   mine_stat(x[,ix[1]], sample(x[,ix[2]]), measure = 'tic', norm=TRUE, est="mic_e", alpha = 9, C=5)
# }
# 
# 
# aa <- matrix(rep(1:10, 4), ncol=4, nrow=10) + matrix(rep((0:3), each=10), ncol=4, nrow=10)
# 
# 
# nperm <- 2500
# set.seed(0)
# microbenchmark(res <- replicate(nperm, mictools_null_r(aa)),
#                re2 <- mictools(aa, nperm=nperm, seed=0))



#' Function for mictools pipeline
#' @aliases mictools
#' @param nperm integer, number of permutation to perform
#' @export
mictools <- function(x, alpha=9, C=5, seed=0, nperm=200000){
  
  bins <- seq(0,1, length.out = 10000 + 1)
  
  ## Compute null distribution for TIC
  ticnull <- mictools_null(x, B=alpha, C=C, nperm=nperm, seed=seed)
  
  ## Compute histogram of the null distribution and right tailed area
  histtic <- hist(ticnull, breaks = bins, plot = FALSE, right=FALSE)
  histcumsum <- rev(cumsum(rev(histtic$counts)))
  
  ## Compute observed values
  ticmicmat <- mictools_pval(x, alpha=alpha, C=C, est="mic_e")
  ix <- t(combn(1:ncol(x),2))
  ticmicmat <- cbind(ticmicmat, ix)
  
  ## Histogram of observed values
  histtictrue <- hist(ticmicmat[,1], breaks = bins, plot=FALSE, right=FALSE)
  obscumtic <- rev(cumsum(rev(histtictrue$counts)))
  
  obsdist <- data.frame(BinStart=bins[1:(length(bins)-1)], BinEnd=bins[2:(length(bins))],
                        Count=histtictrue$counts, CountCum=obscumtic)
  
  ## One-dimensional linear interpolation
  pval <- approx(bins[1:(length(bins)-1)], y = histcumsum, xout=ticmicmat[,1], method="linear", rule=2)$y / (histcumsum[1] + 1)
  pval.df <- data.frame(pval=pval, Var1=ix[,1], Var2=ix[,2])
  
  ## Return values as list
  retlist <- list(tic=ticnull, hist=histtic, histcum=histcumsum, 
                  obstic=ticmicmat[,c(1, 3, 4)], obsdist=obsdist, pval=pval.df)
  return(retlist)
}