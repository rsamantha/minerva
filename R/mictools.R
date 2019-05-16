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
#' @inheritParams mine
#' @param nperm integer, number of permutation to perform
#' @param seed set the seed for random number generation repreoducibility
#' @param p.adjust.method method for pvalue adjustment, see [stats::p.adjust.methods] for further details.
#' @export
mictools <- function(x, alpha=9, C=5, seed=0, nperm=200000, p.adjust.method="BH"){
  
  ## Setup
  bins <- seq(0,1, length.out = 10000 + 1)
  varnames <- colnames(x)
  if (is.null(varnames))
    varnames <- paste("Var", 1:ncol(x))
  
  ## Compute null distribution for TIC
  ticnull <- mictools_null(x, alpha=alpha, C=C, nperm=nperm, seed=seed)
  
  ## Compute histogram of the null distribution and right tailed area
  histtic <- graphics::hist(ticnull, breaks = bins, plot = FALSE, right=FALSE)
  histcumsum <- rev(cumsum(rev(histtic$counts)))
  
  ## Compute observed values
  ticmicmat <- mictools_pstats(x, alpha=alpha, C=C, est="mic_e")
  ix <- t(utils::combn(1:ncol(x),2))
  ticmicmat <- data.frame(TIC=ticmicmat[,2], I1=ix[,1], I2=ix[,2], Var1=varnames[ix[,1]], Var2=varnames[ix[,2]], stringsAsFactors = FALSE)

  ## Histogram of observed values
  histtictrue <- graphics::hist(ticmicmat[,1], breaks = bins, plot=FALSE, right=FALSE)
  obscumtic <- rev(cumsum(rev(histtictrue$counts)))
  obsdist <- data.frame(BinStart=bins[1:(length(bins)-1)], BinEnd=bins[2:(length(bins))],
                        Count=histtictrue$counts, CountCum=obscumtic, stringsAsFactors = FALSE)
  
  ## One-dimensional linear interpolation
  pval <- stats::approx(bins[1:(length(bins)-1)], y = histcumsum, xout=ticmicmat[,1], method="linear", rule=2)$y / (histcumsum[1] + 1)
  pval.df <- data.frame(pval=pval, I1=ix[,1], I2=ix[,2], Var1=varnames[ix[,1]], Var2=varnames[ix[,2]], stringsAsFactors = FALSE)
  pval.df$adj.P.Val <- stats::p.adjust(pval.df$pval, method=p.adjust.method)
  
  ## Return values as list
  retlist <- list(tic=ticnull, hist=histtic, histcum=histcumsum, 
                  obstic=ticmicmat, obsdist=obsdist, pval=pval.df)
  return(retlist)
}


#' This function compute the mic strenght based on statistical significance of the pvalue computed on null distribution of TIC.
#' @param pval a data.frame with pvalues for each pair of association of the x input matrix. It should contain two colums with 
#' the indices of the computed association according to the x input matrix
#' @param pthr threshold on pvalue for measure to consider for computing mic_e
#' @param pval.col a vector of integer or character of length 3 indicating the columns of the pval data.frame with the index of the 
#' matrix x of the association computed. The first index refers to the column containing the p.values, while the remaining 2 indices
#' refers to x variables.
#' @inheritParams mictools
#' @export
mic_strength <- function(x, pval, alpha=NULL, C=5, pthr=0.05, pval.col=1:3)
{
  nbins <- c(1,    25,   50,   250,   500, 1000, 2500, 5000, 10000, 40000)
  alphas <- c(0.85, 0.80, 0.75, 0.70, 0.65, 0.6,  0.55, 0.5,  0.45,  0.4)
  
  if (is.null(alpha))
    alpha <- alphas[.bincode(nrow(x), nbins, include.lowest=TRUE)]
  
  ##
  ##  Check parameters.....
  ##
  
  ## Select association by pvalue
  pix <- pval[, pval.col[1]] < 0.05
  
  ## Compute mic in res[,2]
  res <- sapply(which(pix), function(y, x.sub, vartouse, alpha, C){
    xtmp <- x.sub[,vartouse[y, 1]]
    ytmp <- x.sub[,vartouse[y, 2]]
    
    mine_compute(xtmp, ytmp, alpha=alpha, C=C, est="mic_e", measure = 1)
    
  }, x.sub=x, alpha=alpha, C=C, vartouse=pval[, pval.col[2:3]])
  
  ## Add variable index of the MIC computed
  res.df <- data.frame(TicePval=pval[pix, pval.col[1]], MIC=res, 
                       I1=pval[pix, pval.col[2]], I2=pval[pix, pval.col[3]], 
                       stringsAsFactors = FALSE)
  
  return(res.df)
}