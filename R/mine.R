# This code is written by Michele Filosi  <michele.filosi@gmail.com>
## Roberto Visintainer <r.visintainer@gmail.com>.
## 2012 

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' MINE family statistics
#' Maximal Information-Based Nonparametric Exploration (MINE) 
#' statistics. \code{mine} computes the MINE family measures between two variables. 
#'
#' @aliases mine MINE  MIC-R2  mic-r2
#' 
#' @param x a numeric vector (of size \emph{n}), matrix or data frame (which is coerced to matrix). 
#' @param y  NULL (default) or a numeric vector of size \emph{n} (\emph{i.e.}, with compatible dimensions to x).
#' @param master an optional vector of indices (numeric or character) to
#' be given when \code{y} is not set, otherwise master is ignored. It can be 
#' either one column index to be used as reference for the comparison 
#' (versus all other columns) or a vector of column indices to be used for 
#' computing all mutual statistics.
# If not specified it is set to \code{1:ncol(x)}. 
#' @param alpha an optional number of cells allowed in the \emph{X-by-Y} search-grid. Default value is 0.6 (see Details).
#' @param C an optional number determining the starting point of the 
#' \emph{X-by-Y} search-grid. When trying to partition the \emph{x}-axis into 
#' \emph{X} columns, the algorithm will start with at most \code{C}\emph{X} 
#' \emph{clumps}. Default value is 15 (see Details). 
#' @param n.cores ooptional number of cores to be used in the
#' computations, when master is specified. It requires the
#' \pkg{parallel} package, which provides support for parallel 
#' computing, released with \R >= 2.14.0. Defaults is 1 (\emph{i.e.}, not performing parallel computing). 
#' @param var.thr  minimum value allowed for the variance of the input
#' variables, since \code{mine} can not be computed in case of variance
#' close to 0. Default value is 1e-5. Information about failed check
#' are reported in \emph{var_thr.log} file. 
#' @param eps integer in [0,1].  If 'NULL' (default) it is set to
#' 1-MIC. It can be set to zero for noiseless functions,
#' but the default choice is the most appropriate parametrization
#' for general cases (as stated in Reshef et al. SOM).
#' It provides robustness. 
#' @param est Default value is "mic_approx". With est="mic_approx" the original MINE statistics will
#' be computed, with est="mic_e" the equicharacteristic matrix is
#' is evaluated and the mic() and tic() methods will return MIC_e and
#' TIC_e values respectively.
#' @param na.rm boolean. This variable is passed directly to the
#' \code{cor}-based functions. See \code{cor} for further details.
#' @param use Default value is "all.obs". This variable is passed directly to the 
#' \code{cor}-based functions. See \code{cor} for further details.
#' @param \dots currently ignored
#' @details \code{mine} is an R wrapper for the C engine \emph{cmine} 
#' (\url{http://minepy.readthedocs.io/en/latest/}), 
#' an implementation of Maximal Information-Based Nonparametric Exploration (MINE) 
#' statistics. The MINE statistics were firstly detailed in 
#' D. Reshef et al. (2011) \emph{Detecting novel associations in large datasets}. 
#' Science 334, 6062 (\url{http://www.exploredata.net}).
#' 
#' Here we recall the main concepts of the MINE family statistics.
#' Let \eqn{D={(x,y)}} be the set of \emph{n} ordered pairs of elements of \code{x}
#' and \code{y}. The data space is partitioned in 
#' an \emph{X-by-Y} grid, grouping the \emph{x} and \emph{y} values 
#' in \emph{X} and \emph{Y} bins respectively.\cr
#' 
#' The \strong{Maximal Information Coefficient (MIC)} is defined as 
#' \deqn{\textrm{MIC}(D)=\max_{XY<B(n)} M(D)_{X,Y} = \max_{XY<B(n)} \frac{I^*(D,X,Y)}{log(\min{X,Y})},}{MIC(D)=max_{XY<B(n)} M(D)_{X,Y}=max_{XY<B(n)} I*(D,X,Y)/log(min(X,Y)),} where
#' \eqn{B(n)=n^{\alpha}} is the search-grid size,
#' \eqn{I^*(D,X,Y)}{I*(D,X,Y)}
#' is the maximum mutual information over all grids \emph{X-by-Y}, of the distribution induced by D on 
#' a grid having \emph{X} and \emph{Y} bins (where the probability mass on a cell 
#'                                           of the grid is the fraction of points of D falling in that cell).
#' The other statistics of the MINE family are derived from the mutual information 
#' matrix achieved by an \emph{X-by-Y} grid on D.\cr
#' 
#' The \strong{Maximum Asymmetry Score (MAS)} is defined as
#' \deqn{\textrm{MAS}(D) = \max_{XY<B(n)} |M(D)_{X,Y} - M(D)_{Y,X}|.}{MAS(D) = max_{XY<B(n)} |M(D)_{X,Y} - M(D)_{Y,X}|.}
#' 
#' The \strong{Maximum Edge Value (MEV)} is defined as
#' \deqn{\textrm{MEV}(D) = \max_{XY<B(n)} \{M(D)_{X,Y}: X=2~or~Y=2\}.}{MEV(D) = max_{XY<B(n)} {M(D)_{X,Y}: X=2 or Y=2}.}
#' 
#' The \strong{Minimum Cell Number (MCN)} is defined as
#' \deqn{\textrm{MCN}(D,\epsilon) = \min_{XY<B(n)} \{\log(XY): M(D)_{X,Y} \geq (1-\epsilon)MIC(D)\}.}{MCN(D,\epsilon) = min_{XY<B(n)} {log(XY): M(D)_{X,Y} >= (1-\epsilon)MIC(D)}.}
#' More details are provided in the supplementary material (SOM) of the original paper.
#' 
#' The MINE statistics can be computed for two numeric vectors \code{x} and \code{y}. 
#' Otherwise a matrix (or data frame) can be provided and two options are available 
#' according to the value of \code{master}. If \code{master} is a column identifier, 
#' then the MINE statistics are computed for the \emph{master} variable versus the 
#' other matrix columns. If \code{master} is a set of column identifiers, then all 
#' mutual MINE statistics are computed among the column subset.
#' \code{master}, \code{alpha}, and \code{C} refers respectively to the \emph{style},
#' \emph{exp}, and \emph{c} parameters of the original \emph{java} code. 
#' In the original article, the authors state that the default value \eqn{\alpha=0.6} 
#' (which is the exponent of the search-grid size \eqn{B(n)=n^{\alpha}}) has been 
#' empirically chosen. It is worthwhile noting that \code{alpha} and \code{C} are 
#' defined to obtain an heuristic approximation in a reasonable amount of time. In case
#' of small sample size (\emph{n}) it is preferable to increase \code{alpha} to 1 to 
#' obtain a solution closer to the theoretical one.
#' @return The Maximal Information-Based Nonparametric Exploration (MINE) statistics 
#' provide quantitative evaluations of different aspects of the relationship 
#' between two variables. 
#' In particular \code{mine} returns a list of 5 statistics: 
#'   \item{MIC}{ \strong{Maximal Information Coefficient.} \cr 
#'     It is related to the relationship strenght and it can be interpreted as a 
#'     correlation measure. It is symmetric and it ranges in [0,1], where it 
#'     tends to 0 for statistically independent data and it approaches 1 in 
#'     probability for noiseless functional relationships (more details can 
#'                                                         ben found in the original paper). }
#' \item{MAS}{ \strong{Maximum Asymmetry Score.} \cr 
#'   It captures the deviation from monotonicity. Note that 
#'   \eqn{\textrm{MAS} < \textrm{MIC}}{MAS < MIC}. \cr
#'   \emph{Note:} it can be useful for detecting periodic relationships 
#'   (unknown frequencies). }
#' \item{MEV}{ \strong{Maximum Edge Value.} \cr 
#'   It measures the closeness to being a function. Note that 
#'   \eqn{\textrm{MEV} \leq \textrm{MIC}}{MEV <= MIC}. }
#' \item{MCN}{ \strong{Minimum Cell Number.} \cr 
#'   It is a complexity measure. }
#' \item{MIC-R2}{It is the difference between the MIC value and the Pearson
#'   correlation coefficient. } \cr
#' 
#' When computing \code{mine} between two numeric vectors \code{x} and \code{y}, 
#' the output is a list of 5 numeric values. When \code{master} is provided, 
#' \code{mine} returns a list of 5 matrices having \code{ncol} equal to 
#' \emph{m}. In particular, if \code{master} is a single value, 
#' then \code{mine} returns a list of 5 matrices having 1 column, 
#' whose rows correspond to the MINE measures between the \emph{master} 
#' column versus all. Instead if \code{master} is a vector of \emph{m} indices,
#' then \code{mine} output is a list of 5 \emph{m-by-m} matrices, whose element 
#' \emph{i,j} corresponds to the MINE statistics computed between the \emph{i} 
#' and \emph{j} columns of \code{x}.
#' @references 
#' D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P. Turnbaugh, 
#' E. Lander, M. Mitzenmacher, P. Sabeti. (2011)
#' \emph{Detecting novel associations in large datasets}. 
#' Science 334, 6062\cr
#' \url{http://www.exploredata.net}\cr
#' (SOM: Supplementary Online Material at
#'   \url{http://www.sciencemag.org/content/suppl/2011/12/14/334.6062.1518.DC1})
#' 
#' D. Albanese, M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman,
#' C. Furlanello. 
#' \emph{minerva and minepy: a C engine for the MINE suite and its R, Python and MATLAB wrappers}. 
#' Bioinformatics (2013) 29(3): 407-408, \doi{doi:10.1093/bioinformatics/bts707}.\cr
#' 
#' \emph{minepy. Maximal Information-based Nonparametric Exploration in C and Python.}\cr 
#' \url{http://minepy.sourceforge.net}
#' 
#' @author 
#' Michele Filosi and Roberto Visintainer 
#' @examples 
#' A <- matrix(runif(50),nrow=5)
#' mine(x=A, master=1)
#' mine(x=A, master=c(1,3,5,7,8:10))
#' 
#' x <- runif(10); y <- 3*x+2; plot(x,y,type="l")
#' mine(x,y)
#' # MIC = 1 
#' # MAS = 0
#' # MEV = 1
#' # MCN = 2
#' # MIC-R2 = 0
#' 
#' set.seed(100); x <- runif(10); y <- 3*x+2+rnorm(10,mean=2,sd=5); plot(x,y)
#' mine(x,y)
#' # rounded values of MINE statistics
#' # MIC = 0.61
#' # MAS = 0
#' # MEV = 0.61
#' # MCN = 2
#' # MIC-R2 = 0.13
#' 
#' t <-seq(-2*pi,2*pi,0.2); y1 <- sin(2*t); plot(t,y1,type="l")
#' mine(t,y1)
#' # rounded values of MINE statistics
#' # MIC = 0.66 
#' # MAS = 0.37
#' # MEV = 0.66
#' # MCN = 3.58
#' # MIC-R2 = 0.62
#' 
#' y2 <- sin(4*t); plot(t,y2,type="l")
#' mine(t,y2)
#' # rounded values of MINE statistics
#' # MIC = 0.32 
#' # MAS = 0.18
#' # MEV = 0.32
#' # MCN = 3.58
#' # MIC-R2 = 0.31
#' 
#' # Note that for small n it is better to increase alpha
#' mine(t,y1,alpha=1)
#' # rounded values of MINE statistics
#' # MIC = 1 
#' # MAS = 0.59
#' # MEV = 1
#' # MCN = 5.67
#' # MIC-R2 = 0.96
#' 
#' mine(t,y2,alpha=1)
#' # rounded values of MINE statistics
#' # MIC = 1 
#' # MAS = 0.59
#' # MEV = 1
#' # MCN = 5
#' # MIC-R2 = 0.99
#' 
#' # Some examples from SOM
#' x <- runif(n=1000, min=0, max=1)
#' 
#' # Linear relationship
#' y1 <- x; plot(x,y1,type="l"); mine(x,y1)
#' # MIC = 1 
#' # MAS = 0
#' # MEV = 1
#' # MCN = 4
#' # MIC-R2 = 0
#' 
#' # Parabolic relationship
#' y2 <- 4*(x-0.5)^2; plot(sort(x),y2[order(x)],type="l"); mine(x,y2)
#' # rounded values of MINE statistics
#' # MIC = 1 
#' # MAS = 0.68
#' # MEV = 1
#' # MCN = 5.5
#' # MIC-R2 = 1
#' 
#' # Sinusoidal relationship (varying frequency)
#' y3 <- sin(6*pi*x*(1+x)); plot(sort(x),y3[order(x)],type="l"); mine(x,y3)
#' # rounded values of MINE statistics
#' # MIC = 1 
#' # MAS = 0.85
#' # MEV = 1
#' # MCN = 4.6
#' # MIC-R2 = 0.96
#' 
#' # Circle relationship
#' t <- seq(from=0,to=2*pi,length.out=1000)
#' x4 <- cos(t); y4 <- sin(t); plot(x4, y4, type="l",asp=1)
#' mine(x4,y4)
#' # rounded values of MINE statistics
#' # MIC = 0.68 
#' # MAS = 0.01
#' # MEV = 0.32
#' # MCN = 5.98
#' # MIC-R2 = 0.68
#' 
#' data(Spellman)
#' res <- mine(Spellman,master=1,n.cores=1)
#' 
#' \dontrun{## example of multicore computation
#' res <- mine(Spellman,master=1,n.cores=parallel::detectCores()-1)}
#' @export
mine <- function(x, y=NULL, master=NULL, alpha=0.6, C=15, n.cores=1, var.thr=1e-5, eps=NULL, est="mic_approx",
                 na.rm=FALSE, use="all.obs", ...){

  ## Control on input arguments
  checked <- check.inputs(x,y,alpha,C,n.cores,var.thr,eps, est, na.rm, use)
  x <- checked[[1]]
  y <- checked[[2]]
  alpha <- checked[[3]]
  C <- checked[[4]]
  n.cores <- checked[[5]]
  eps <- checked[[6]]
  est <- checked[[7]]
  var.idx <- checked[[8]]
  na.rm <- checked[[9]]
  use <- checked[[10]]
  
  ## only one matrix given
  if (is.null(y)){
    s <- dim(x)[1]
    f <- dim(x)[2]

    ## Check for na and overwrite the input
    ## No need to store the input
    if (na.rm | use %in% c(2L, 3L)){
      x <- na.omit(x)
    }
    
    ## If a master variable is passed check the type
    if (!is.null(master)){
      if (!is.numeric(master))
        stop("'master' must be numeric")
      if (max(master)>f)
        stop("Subscript out of bound!\nMaximum value allowed: ",f)
      if (min(master)<1)
        stop("Subscript out of bound!\nMinimum value allowed: ",1)
    }
    if (is.null(master)){
      if (n.cores>1){
        ## Launch parallel
        res <- .allvsallparall(x,alpha,C,n.cores,eps,est)
        ## return(.allvsallparall(x,alpha,C,n.cores,eps))
      } else{
        res <- .allvsall(x,alpha,C,eps,est)
        ## return(.allvsall(x,alpha,C,eps))
      }
    } else {
      if (length(master)==1){
        if (n.cores>1){
          res <- .onevsallparall(x,master,alpha,C,n.cores,eps,est)
          ## return(.onevsallparall(x,master,alpha,C,n.cores,eps))
        }
        else{
          res <- .onevsall(x,master,alpha,C,exclude=FALSE,eps=eps, est=est)
          ## return(.onevsall(x,master,alpha,C,exclude=FALSE,eps))
        }
      }
      if (length(master)>1){
        newdata <- x[,master]
        if (n.cores>1){
          ## Launch parallel
          res <- .allvsallparall(newdata,alpha,C,n.cores,eps,est)
          ## return(.allvsallparall(newdata,alpha,C,n.cores,eps))
        } else {
          res <- .allvsall(newdata,alpha,C,eps,est)
          ## return(.allvsall(newdata,alpha,C,eps))
        }
      }
    }
  } else { ## y is given

    ## pairwise.complete.obs // complete.obs
    if (na.rm | use%in% c(2L,3L)){
      xidx <- (apply(!is.na(x), 1, sum))
      yidx <- !is.na(y)
      idx <- xidx & yidx
      x <- as.matrix(x[idx,])
      y <- as.matrix(y[idx,])
    }
    
    ## two variables given
    if (ncol(x) == 1 & ncol(y) == 1){
      res <- .Call("mineRonevar",as.double(x),as.double(y),alpha=alpha,C=C,eps=eps, est=est)
      names(res) <- c("MIC","MAS","MEV","MCN","MIC-R2", "GMIC", "TIC")
      res <- as.list(res)
    } else {
      newdata <- cbind(x,y)
      colnames(newdata)[ncol(newdata)] <- "Y"
      res <- .onevsall(newdata,ncol(newdata),alpha,C,exclude=TRUE,eps=eps,est=est)
      res <- lapply(res, function(x){dimnames(x) <- list(NULL, colnames(x)); return(x)})
    }
  }
  
  ## Set NA variables with nearly 0 variance
  ## TODO: set NA with nearly 0 variance when y is passed as arguments
  if (!is.null(var.idx[["x"]]))
    res <- lapply(res,
                  function(x,var.idx){
                    if(is.null(dim(x))){
                      x <- 0.0
                    } else {
                      if (length(master)==1) {
                        if (any(var.idx==master)){
                          x[,1] <- 0.0
                        } else {
                          x[var.idx,] <- 0.0
                        }
                      } else {
                        x[var.idx,] <- 0.0
                      }
                    }
                    return(x)},
                  var.idx=var.idx[["x"]])

  ## Return results
  return(res)
}

##--------------------------------------------------
## Functions for internal use!
##--------------------------------------------------

## Checking input integrity
## x should be a matrix or a vector
##   if x is a vector y should be given
## y should be a one dimensional vector
check.inputs <- function(x,y,alpha,C,n.cores,var.thr,eps,est,na.rm,use) {

  ## MINE parameters check!
  if (alpha<=0.0 || alpha>1.0 || !is.numeric(alpha))
    stop("'alpha' must be in (0.0, 1.0]",call.=FALSE)
  if(C<=0.0 || !is.numeric(C))
    stop("'C' must be > 0.0",call.=FALSE)

  ## use argument
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything"))
  if (is.na(na.method))
    stop("invalid 'use' argument")
  
  ## Data check!
  if (is.data.frame(x))
    x <- as.matrix(x)

  ## Check for NA/missing values
  if (any(is.na(x)) & !(na.rm | na.method %in% c(2L,3L))){
    stop("Missing values present in input variable 'x'. Consider using use = 'pairwise.complete.obs'.", call.=FALSE)
  }
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'", call.=FALSE)
  if (!is.numeric(x))
    stop("'x' must be numeric", call.=FALSE)
  stopifnot(is.atomic(x))
  x <- as.matrix(x)

  ## Check variance
  var.idx <- list()
  myvar <- apply(x,2,var,na.rm=TRUE)
  if (sum(myvar<var.thr)!=0)
    var.idx[["x"]] <- which(myvar<var.thr)
  
  ## y not empty
  if (!is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric", call.=FALSE)
    stopifnot(is.atomic(y))
    y <- as.matrix(y)

    ## Check for NA/missing values
    if (any(is.na(y)) & !(na.rm | na.method %in% c(2L,3L))){
      stop("Missing values present in input variable 'y'. Consider using use = 'pairwise.complete.obs'.", call.=FALSE)
    }
    
    ## Check variance on argument y
    if (var(y, na.rm=TRUE)<var.thr)
      var.idx[["y"]] <- c(var.idx,1)

    ## Control dimensions of y
    if (dim(y)[2] != 1)
      stop("'y' must be a 1-dimensional object", call.=FALSE)

    ## Consistency check between x and y
    if (nrow(y) != nrow(x))
      stop ("'x' and 'y': incompatible dimensions", call.=FALSE)
  }
  
  ## Check for parallel computing
  if (!is.numeric(n.cores))
    stop("'n.cores' must be numeric")

  if (!exists("detectCores")){
    n.cores <- 1
    warning("It seems you have an old version of R (<2.14).\nMulticore computing using 'parallel' is not available.\n'n.cores' has been set to 1",
            immediate.=TRUE)
  }
  
  if (n.cores > 1){
    if (detectCores()==1){
      warning("Only 1 core available on this machine:'n.cores' will be ignored.\nPlease consider using a cluster.",
              immediate.=TRUE)
      n.cores <- 1
    }
    if (n.cores>detectCores())
      stop("Query for more cores then the available ones.\n Max number of cores for this machine is: ",detectCores())
  }
  if (n.cores <= 0){
    stop("You are trying to compute mic using ",n.cores," cores.. are you sure?", call.=FALSE)
  }
  
  ## Check epsilon parameter for MINE
  if(!is.null(eps)){
    if( eps<0.0 || eps>1 || !is.numeric(eps) )
      stop("'eps' must be > 0.0 and < 1.0",call.=FALSE)
  }  

  ## Check est parameter for MIC
  ## mic_approx -> 0
  ## mic_e -> 1
    
  EST <- 0
  if (est == "mic_e"){
      EST <- 1
  }
       
  ## Send a warning if variable with 0 variance have been found
  if (!is.null(var.idx[["x"]]) || !is.null(var.idx[["y"]])){
    warning("Found variables with nearly 0 variance")
  }

  ## Check if na.rm is a boolean or a integer
  if (length(na.rm)>1 || (!is.logical(na.rm) && !is.integer(na.rm)))
    na.rm <- FALSE
  
  return(list(x=x, y=y, alpha=alpha, C=C, n.cores=n.cores, eps=eps, est=EST, var.idx=var.idx, na.rm=na.rm, use=na.method))
}


## Calling all features vs all features using C implementation
## For the source see src/mine_interface.c
.allvsall <- function(x, alpha, C,eps, est){
  tmp <- .Call("mineRall",x,nrow(x),ncol(x),alpha,C, eps, est)
  tmp <- lapply(tmp,function(y,n){colnames(y) <- rownames(y) <- n
                                  return(y)}, n <- colnames(x))
  return(tmp)
}

## Calling feature x[,idx] vs all other features
## Using C implementation feature vs feature
## For the source see src/mine_interface.c
.onevsall <- function(x,idx,alpha,C,eps,exclude,est,diagonal=FALSE){
  if (exclude)
    f <- dim(x)[2]-1
  else
    f <- dim(x)[2]
  if (diagonal)
    start <- idx
  else
    start <- 1
  
  Mat1 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[idx]))
  Mat2 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[idx]))
  Mat3 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[idx]))
  Mat4 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[idx]))
  Mat5 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[idx]))
  Mat6 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[idx]))
  Mat7 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[idx]))
  
  for (i in start:f){
    res <- .Call("mineRonevar",as.double(x[,idx]),as.double(x[,i]),
                 alpha=alpha,C=C,eps=eps, est=est)
    names(res) <- c("MIC","MAS","MEV","MCN","MIC-R2", "GMIC", "TIC")
    Mat1[i,1] <- res["MIC"]
    Mat2[i,1] <- res["MAS"]
    Mat3[i,1] <- res["MEV"]
    Mat4[i,1] <- res["MCN"]
    Mat5[i,1] <- res["MIC-R2"]
    Mat6[i,1] <- res["GMIC"]
    Mat7[i,1] <- res["TIC"]
  }
  return(list(MIC=Mat1,MAS=Mat2,MEV=Mat3,MCN=Mat4,MICR2=Mat5,GMIC=Mat6, TIC=Mat7))
}

## Parallel implementation of one vs all function
## NB using 'parallel' package from CRAN for R >= 2.14
## If older version of R install multicore package
.onevsallparall <- function(x,master,alpha,C,n.cores,eps, est){
  f <- dim(x)[2]
  cl <- makeCluster(n.cores)
  res <- parLapply(cl,1:f,function(i,master,alpha,C,data,eps,est){
    return(.Call("mineRonevar",as.double(data[,master]),
                 as.double(data[,i]),alpha=alpha,C=C,eps=eps,est=est))},
                   master=master,alpha=alpha,C=C,eps=eps,est=est,data=x)
  stopCluster(cl)
  
  Mat1 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat2 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat3 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat4 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat5 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat6 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat7 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  
  for (i in 1:f){
    Mat1[i,1] <- res[[i]][1]
    Mat2[i,1] <- res[[i]][2]
    Mat3[i,1] <- res[[i]][3]
    Mat4[i,1] <- res[[i]][4]
    Mat5[i,1] <- res[[i]][5]
    Mat6[i,1] <- res[[i]][6]
    Mat7[i,1] <- res[[i]][7]
  }
  return(list(MIC=Mat1,MAS=Mat2,MEV=Mat3,MCN=Mat4,MICR2=Mat5,GMIC=Mat6,TIC=Mat7))
}

## Parallel implementation of all vs all function
## NB using 'parallel' package from CRAN for R >= 2.14
## If older version of R install multicore package
.allvsallparall <- function(x, alpha, C, n.cores,eps,est){
  f <- dim(x)[2]
  cl <- makeCluster(n.cores)
  res <- parLapply(cl,1:f,function(y,data,alpha,C,eps,est){
    return(.onevsall(x=data,idx=y,alpha=alpha,C=C,eps=eps,est=est,exclude=FALSE,diagonal=TRUE))},
                   data=x,alpha=alpha,C=C,eps=eps,est=est)
  
  stopCluster(cl)
  Mat1 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat2 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat3 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat4 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat5 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat6 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat7 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  
  for (i in seq(length(res))){

    Mat1[i,i:f] <- res[[i]][[1]][i:f,]
    Mat1[i:f,i] <- res[[i]][[1]][i:f,]
    
    Mat2[i,i:f] <- res[[i]][[2]][i:f,]
    Mat2[i:f,i] <- res[[i]][[2]][i:f,]

    Mat3[i,i:f] <- res[[i]][[3]][i:f,]
    Mat3[i:f,i] <- res[[i]][[3]][i:f,]
    
    Mat4[i,i:f] <- res[[i]][[4]][i:f,]
    Mat4[i:f,i] <- res[[i]][[4]][i:f,]

    Mat5[i,i:f] <- res[[i]][[5]][i:f,]
    Mat5[i:f,i] <- res[[i]][[5]][i:f,]

    Mat6[i,i:f] <- res[[i]][[6]][i:f,]
    Mat6[i:f,i] <- res[[i]][[6]][i:f,]

    Mat7[i,i:f] <- res[[i]][[7]][i:f,]
    Mat7[i:f,i] <- res[[i]][[7]][i:f,]

  }
  return(list(MIC=Mat1,MAS=Mat2,MEV=Mat3,MCN=Mat4,MICR2=Mat5,GMIC=Mat6,TIC=Mat7))
}


# .onUnload <- function (libpath) {
#   library.dynam.unload("minerva", libpath)
# }