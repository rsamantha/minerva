#' Function that implements the \code{mictools} pipeline.
#' In particular it computes the null and observed distribution of the \code{tic_e} measure
#' 
#' @aliases mictools
#' @param x a numeric matrix with N samples on the rows and M variables on the columns (NxM). 
#' @param C a positive integer number, the \code{C} parameter of the \code{mine} statistic. 
#' See \code{\link[minerva]{mine}} function for further details.
#' @param nperm integer, number of permutation to perform
#' @param seed seed for random number generation reproducibility
#' @param p.adjust.method method for pvalue adjustment, see \code{\link[stats]{p.adjust}} for available methods.
#' @inheritParams mine
#' @details This is a function to implement the `mictools` pipeline.
#' Differently from the python pipeline available on github we consider a data matrix of NxM with N samples by rows 
#' and M variables by columns as standard for R.
#' @return A list of 5 named elements containing the following information of the computed statistic:
#' \describe{
#' \item{tic}{This is a vector with the null distribution of tic_e values based on the permutation.}
#' \item{nulldist}{Null distribution of the \code{tic_e} measure. It is a \code{data.frame} of 4 columns
#'  containing the histogram of the distribution of \code{tic_e} for each bin delimited by \code{BinStart} 
#'  and \code{BinEnd}, the count for each bin \code{NullCount} and the cumulative distribution 
#'  of the right tail area \code{NullCumSum}}
#' \item{obstic}{\code{data.frame} with the observed \code{tic_e} values, the indexes of the variables between the tic is computed. 
#' If the input matrix has column names then the names are reported in the dataframe, otherwise "Var<i>" is added for each variable.}
#' \item{obsdists}{\code{data.frame} similar to \code{nulldist} but with observed values of \code{tic_e}}
#' \item{pval}{data.frame with the pvalue computed for each comparison. The adjusted pvalue is also reported based 
#' on the method chosen with the parameter \code{p.adjust.method}}
#' }
#' @seealso \code{\link[stats]{p.adjust}}, \code{\link[graphics]{hist}}, \code{\link[minerva]{mine}}
#' @references 
#' D. Albanese, S. Riccadonna, C. Donati, P. Franceschi (2018)
#' _A practical tool for Maximal Information Coefficient Analysis_
#' GigaScience, 7, 4, \url{https://doi.org/10.1093/gigascience/giy032}
#' @examples 
#' data(Spellman)
#' Spellman <- as.matrix(Spellman)
#' spellress <- mictools(Spellman[, 10:20], nperm=1000)
#' 
#' ## Use a different pvalue correction method
#' spellressb <- mictools(Spellman[,10:20], nperm=1000, seed=1234, p.adjust.method="bonferroni")
#' 
#' ## Distribution of tic_e null
#' hist(spellress$tic, breaks=100, main="Tic_e null distribution")
#' barplot(spellress$nulldist$NullCount)
#' 
#' ## Distribution of the observed tic
#' hist(spellress$obstic$TIC)
#' barplot(spellress$obsdist$Count)
#' 
#' ## Distribution of empirical pvalues
#' hist(spellress$pval$pval, breaks=50)
#' 
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
  ## ticmicmat <- mictools_pstats(x, alpha=alpha, C=C, est="mic_e")
  ticmicmat <- pstats(x, alpha=alpha, C=C, est="mic_e")
  
  # ix <- t(utils::combn(1:ncol(x),2))
  ticmicmat <- data.frame(TIC=ticmicmat[, 4], I1=ticmicmat[, 1], I2=ticmicmat[,2], 
                          VarNames1=varnames[ticmicmat[, 1]], VarNames2=varnames[ticmicmat[, 2]], 
                          stringsAsFactors = FALSE)
  
  nulldist <- data.frame(BinStart=bins[1:(length(bins)-1)], BinEnd=bins[2:(length(bins))], 
                         NullCount=histtic$counts, NullCumSum=histcumsum, stringsAsFactors = FALSE)
  
  ## Histogram of observed values
  histtictrue <- graphics::hist(ticmicmat[,1], breaks = bins, plot=FALSE, right=FALSE)
  obscumtic <- rev(cumsum(rev(histtictrue$counts)))
  obsdist <- data.frame(BinStart=bins[1:(length(bins)-1)], BinEnd=bins[2:(length(bins))],
                        Count=histtictrue$counts, CountCum=obscumtic, stringsAsFactors = FALSE)
  
  ## One-dimensional linear interpolation
  pval <- stats::approx(bins[1:(length(bins)-1)], y = histcumsum, xout=ticmicmat[,1], method="linear", rule=2)$y / (histcumsum[1] + 1)
  pval.df <- data.frame(pval=pval, I1=ticmicmat[,2], I2=ticmicmat[,3], 
                        Var1=varnames[ticmicmat[,2]], Var2=varnames[ticmicmat[,3]], 
                        stringsAsFactors = FALSE)
  pval.df$adj.P.Val <- stats::p.adjust(pval.df$pval, method=p.adjust.method)
  
  ## Return values as list
  retlist <- list(tic=ticnull, nulldist=nulldist, 
                  obstic=ticmicmat, obsdist=obsdist, pval=pval.df)
  return(retlist)
}


#' Compute the association strengh 
#' 
#' This function uses the null distribution of the \code{tic_e} computed with the function \code{\link[minerva]{mictools}}. 
#' Based on the available pvalue and the permutation null distribution it identifies reliable association between variables.
#'
#' 
#' @param pval a data.frame with pvalues for each pair of association of the \code{x} input matrix. It should contain two colums with 
#' the indices of the computed association according to the x input matrix
#' @param pthr threshold on pvalue for measure to consider for computing mic_e
#' @param pval.col an integer or character or vector relative to the columns of \code{pval} dataframe respectively for \code{pvalue}, 
#' association between variable 1, variable 2 in the \code{x} input matrix. See Details for further information.
#' @inheritParams mictools
#' @details The method implemented here is a wrapper for the original method published by Albaese et al. (2018). The python version
#' is available at \url{https://github.com/minepy/mictools}.
#' 
#' This function should be called after the estimation of the null distribution of \code{tic_e} scores based on permutations of the input data.
#' 
#' The \code{mic} association is computed only for the variables for which the pvalue in the \code{pval} \code{data.frame} is less then 
#' the threshold set with the \code{pthr} input parameter. 
#' We assume the first column of the \code{pval} \code{data.frame} contains the pvalue, this value can be changed using 
#' the \code{pval.col}[1] parameter. 
#' 
#' The \code{pval.col} parameter, by default takes the first three columns in the \code{pval} \code{data.frame}, in particular the first column containing the \code{pvalues} 
#' of the association between variable in column \code{pval.col[2]} and \code{pval.col[3]}.
#' If a character vector is provided names in \code{pval.col} are matched with the names in \code{pval} \code{data.frame}.
#' If \code{NULL} is passed it is assumed the first column contains pvalue, while the 2 and 3 the index or name of the variable in \code{x}.
#' If one value is passed it refers to the \code{pvalue} column and the consecutive two columns are assume to contain variable indexes.
#' 
#' @return A dataframe with the \code{tic_e} Pvalue, the \code{mic} value and the column identifier regarding the input matrix
#' \code{x} of the variables of which the association is computed.
#' @seealso \code{\link[minerva]{mine}}, \code{\link[minerva]{mictools}}, \code{\link[stats]{p.adjust}}
#' @examples
#' data(Spellman)
#' mydata <- as.matrix(Spellman[, 10:20])
#' ticenull <- mictools(mydata, nperm=1000)
#' 
#' ## Use the nominal pvalue:
#' ms <- mic_strength(mydata, pval=ticenull$pval, alpha=NULL, pval.col = c(6, 4,5))
#' 
#' ## Use the adjusted pvalue:
#' ms <- mic_strength(mydata, pval=ticenull$pval, alpha=NULL, pval.col = c(6, 4,5))
#' 
#' ms 
#' 
#' \dontrun{
#' ## Use qvalue
#' require(qvalue)
#' qobj <- qvalue(ticenull$pval$pval)
#' ticenull$pval$qvalue <- qobj$qvalue
#' ms <- mic_strength(mydata, pval=ticenull$pval, alpha=NULL, pval.col = c("qvalue", "Var1", "Var2"))
#' 
#' ## Get the data from mictools repository
#'
#' lnf <- "https://raw.githubusercontent.com/minepy/mictools/master/examples/datasaurus.txt"
#' datasaurus <- read.table(lnf, header=TRUE, row.names = 1, stringsAsFactors = FALSE)
#' datasaurus <- t(datasaurus)
#' ticenull <- mictools(datasaurus, nperm=200000)
#' micres <- mic_strength(mydata, ticenull$pval, pval.col=c(6, 4, 5))
#' 
#' ## Plot distribution of pvalues
#' hist(ticenull$pval, breaks=50, freq=FALSE)
#' 
#' ## Plot distribution of tic_e values
#' hist(ticenull$tic)
#' 
#' ## Correct pvalues using qvalue package
#' require(qvalue)
#' require(ggplot2)
#' qobj <- qvalue(ticenull$pval$pval)
#' ticenull$pval$qvalue <- qobj$qvalue
#' micres <- mic_strength(datasaurus, ticenull$pval, pval.col=c("qvalue", "Var1", "Var2"))
#' 
#' hist(qobj$qvalue)
#' 
#' df <- data.frame(pi0.labmda=qobj$pi0.lambda, lambda=qobj$lambda, pi0.smooth=qobj$pi0.smooth)
#' gp0 <- ggplot(df, aes(lambda, pi0.labmda)) + geom_point() 
#' gp0 <- gp0 + geom_line(aes(lambda, pi0.smooth))
#' gp0 <- gp0 + geom_hline(yintercept = qobj$pi0, linetype="dashed", col="red")
#' }
#' @export
mic_strength <- function(x, pval, alpha=NULL, C=5, pthr=0.05, pval.col=NULL)
{
  nbins <- c(1,    25,   50,   250,   500, 1000, 2500, 5000, 10000, 40000)
  alphas <- c(0.85, 0.80, 0.75, 0.70, 0.65, 0.6,  0.55, 0.5,  0.45,  0.4)
  
  if (is.null(alpha))
    alpha <- alphas[.bincode(nrow(x), nbins, include.lowest=TRUE)]
  
  ##  Check parameters.....
  pval.col <- check_pvalcol(pval, pval.col)
  
  ## Select association by pvalue
  pix <- pval[, pval.col[1]] < 0.05
  
  ## Compute mic in res[,2]
  res <- sapply(which(pix), function(y, x.sub, vartouse, alpha, C){
    xtmp <- x.sub[,vartouse[y, 1]]
    ytmp <- x.sub[,vartouse[y, 2]]
    
    mine_stat(xtmp, ytmp, alpha=alpha, C=C, est="mic_e", measure = "mic")
    
  }, x.sub=x, alpha=alpha, C=C, vartouse=pval[, pval.col[2:3]])
  
  ## Add variable index of the MIC computed
  res.df <- data.frame(TicePval=pval[pix, pval.col[1]], MIC=res, 
                       I1=pval[pix, pval.col[2]], I2=pval[pix, pval.col[3]], 
                       stringsAsFactors = FALSE)
  
  return(res.df)
}

## Helper function to check for correct pval.col parameter
## Regarding pval dataframe
check_pvalcol <- function(pval, pval.col)
{
  if (length(pval.col) > 3)
    pval.col <- pval.col[1:3]
  
  if (is.character(pval.col)){
    mx <- match(pval.col, names(pval))
    pval.col <- mx[!is.na(mx)]
  }
  
  if ((length(pval.col) == 0) || (is.null(pval.col)))
    pval.col <- 1:3
  
  if (length(pval.col) == 1)
    pval.col <- seq(pval.col, pval.col+2)
  
  if (length(pval.col)==2)
    pval.col <- c(pval.col, pval.col[length(pval.col)] + 1)
  
  if (max(pval.col) > ncol(pval))
    stop(paste0("Not available column "), max(pval.col), call. = FALSE)
  
  return(pval.col)
}