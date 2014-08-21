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

mine <- function(x, y=NULL, master=NULL, alpha=0.6, C=15, n.cores=1, var.thr=1e-5, eps=NULL, ...){
  ## Input controls
  checked <- check.inputs(x,y,alpha,C,n.cores,var.thr,eps)
  x <- checked[[1]]
  y <- checked[[2]]
  alpha <- checked[[3]]
  C <- checked[[4]]
  n.cores <- checked[[5]]
  eps <- checked[[6]]
  var.idx <- checked[[7]]
  
  ## only one matrix given
  if (is.null(y)){
    s <- dim(x)[1]
    f <- dim(x)[2]
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
        res <- .allvsallparall(x,alpha,C,n.cores,eps)
        ## return(.allvsallparall(x,alpha,C,n.cores,eps))
      } else{
        res <- .allvsall(x,alpha,C,eps)
        ## return(.allvsall(x,alpha,C,eps))
      }
    } else {
      if (length(master)==1){
        if (n.cores>1){
          res <- .onevsallparall(x,master,alpha,C,n.cores,eps)
          ## return(.onevsallparall(x,master,alpha,C,n.cores,eps))
        }
        else{
          res <- .onevsall(x,master,alpha,C,exclude=FALSE,eps)
          ## return(.onevsall(x,master,alpha,C,exclude=FALSE,eps))
        }
      }
      if (length(master)>1){
        newdata <- x[,master]
        if (n.cores>1){
          ## Launch parallel
          res <- .allvsallparall(newdata,alpha,C,n.cores,eps)
          ## return(.allvsallparall(newdata,alpha,C,n.cores,eps))
        } else {
          res <- .allvsall(newdata,alpha,C,eps)
          ## return(.allvsall(newdata,alpha,C,eps))
        }
      }
    }
  } else {
    ## two variables given
    if (ncol(x) == 1 && ncol(y) == 1){
      res <- .Call("mineRonevar",as.double(x),as.double(y),alpha=alpha,C=C,eps=eps)
      names(res) <- c("MIC","MAS","MEV","MCN","MIC-R2")
      res <- as.list(res)
      ## return(as.list(res))
    } else {
      newdata <- cbind(x,y)
      colnames(newdata)[ncol(newdata)] <- "Y"
      res <- .onevsall(newdata,ncol(newdata),alpha,C,exclude=TRUE,eps)
      ## return(.onevsall(newdata,ncol(newdata),alpha,C,exclude=TRUE,eps))
    }
  }
  ## Set NA variables with nearly 0 variance
  if (!is.null(var.idx))
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
                        x[var.idx,] <- x[,var.idx] <- 0.0
                      }
                    }
                    return(x)},
                  var.idx=var.idx)
  
  return(res)
}

##--------------------------------------------------
## Function for internal use!
##--------------------------------------------------

## Checking input integrity
## x should be a matrix or a vector
##   if x is a vector y should be given
## y should be a one dimensional vector
check.inputs <- function(x,y,alpha,C,n.cores,var.thr,eps) {

  ## MINE parameters check!
  if (alpha<=0.0 || alpha>1.0 || !is.numeric(alpha))
    stop("'alpha' must be in (0.0, 1.0]",call.=FALSE)
  if(C<=0.0 || !is.numeric(C))
    stop("'C' must be > 0.0",call.=FALSE)
  
  ## Data check!
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (any(is.na(x))){
    nas <- sum(is.na(x))
    stop(nas," NAs found in 'x', please, consider imputing or remove them.", call.=FALSE)
  }
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'", call.=FALSE)
  if (!is.numeric(x))
    stop("'x' must be numeric", call.=FALSE)
  stopifnot(is.atomic(x))
  x <- as.matrix(x)

  ## Check variance
  ## NB modified to return NA -> check on the mine function
  var.idx <- NULL
  if (sum(apply(x,2,var)<var.thr)!=0){
    var.idx <- which(apply(x,2,var)<var.thr)
    ## write.table(var.idx,"var_thr_x.log",col.names=FALSE,row.names=FALSE,quote=FALSE,sep=",")
    ## stop("Found columns having variance < ", var.thr,"; in 'x';\nread file var_thr_x.log for more information.\n")
  }
  if (!is.null(y)) {
    if (!is.numeric(y))
      stop("'y' must be numeric", call.=FALSE)
    stopifnot(is.atomic(y))
    y <- as.matrix(y)
    if (any(is.na(y))){
      nas <- sum(is.na(y))
      stop(nas," NAs found in 'y', please, consider imputing or remove them.", call.=FALSE)
    }
    ## NB remove the following comments to get back to the previous version
    if (var(y)<var.thr)
      var.idx <- 1
      ## stop("'y' has variance < ", var.thr,"\n")
    
    if (dim(y)[2] != 1)
      stop("'y' must be a 1-dimensional object", call.=FALSE)
    if (nrow(y) != nrow(x))
      stop ("'x' and 'y': incompatible dimensions", call.=FALSE)
  }

  ## Parallel computation check!
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
  
  if(!is.null(eps)){
    if( eps<0.0 || eps>1 || !is.numeric(eps) )
      stop("'eps' must be > 0.0 and < 1.0",call.=FALSE)
  }  

  if (!is.null(var.idx)){
    warning("Found variables with nearly 0 variance")
  }
  
  return(list(x, y, alpha, C, n.cores, eps, var.idx))
}


## Calling all features vs all features using C implementation
## For the source see src/mine_interface.c
.allvsall <- function(x, alpha, C,eps){
  tmp <- .Call("mineRall",x,nrow(x),ncol(x),alpha,C, eps)
  tmp <- lapply(tmp,function(y,n){colnames(y) <- rownames(y) <- n
                                  return(y)}, n <- colnames(x))
  return(tmp)
}

## Calling feature x[,idx] vs all other features
## Using C implementation feature vs feature
## For the source see src/mine_interface.c
.onevsall <- function(x,idx,alpha,C,eps,exclude,diagonal=FALSE){
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
  
  for (i in start:f){
    res <- .Call("mineRonevar",as.double(x[,idx]),as.double(x[,i]),
                 alpha=alpha,C=C,eps=eps,package="minerva")
    names(res) <- c("MIC","MAS","MEV","MCN","MIC-R2")
    Mat1[i,1] <- res["MIC"]
    Mat2[i,1] <- res["MAS"]
    Mat3[i,1] <- res["MEV"]
    Mat4[i,1] <- res["MCN"]
    Mat5[i,1] <- res["MIC-R2"]
  }
  return(list(MIC=Mat1,MAS=Mat2,MEV=Mat3,MCN=Mat4,MICR2=Mat5))
}

## Parallel implementation of one vs all function
## NB using 'parallel' package from CRAN for R >= 2.14
## If older version of R install multicore package
.onevsallparall <- function(x,master,alpha,C,n.cores,eps){
  f <- dim(x)[2]
  cl <- makeCluster(n.cores)
  res <- parLapply(cl,1:f,function(i,master,alpha,C,data,eps){
    return(.Call("mineRonevar",as.double(data[,master]),
                 as.double(data[,i]),alpha=alpha,C=C,eps=eps,package="minerva"))},
                   master=master,alpha=alpha,C=C,eps=eps,data=x)
  stopCluster(cl)
  
  Mat1 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat2 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat3 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat4 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))
  Mat5 <- matrix(0,nrow=f,ncol=1,dimnames=list(colnames(x)[1:f],colnames(x)[master]))

  for (i in 1:f){
    Mat1[i,1] <- res[[i]][1]
    Mat2[i,1] <- res[[i]][2]
    Mat3[i,1] <- res[[i]][3]
    Mat4[i,1] <- res[[i]][4]
    Mat5[i,1] <- res[[i]][5]
  }
  return(list(MIC=Mat1,MAS=Mat2,MEV=Mat3,MCN=Mat4,MICR2=Mat5))
}

## Parallel implementation of all vs all function
## NB using 'parallel' package from CRAN for R >= 2.14
## If older version of R install multicore package
.allvsallparall <- function(x, alpha, C, n.cores,eps){
  f <- dim(x)[2]
  cl <- makeCluster(n.cores)
  res <- parLapply(cl,1:f,function(y,data,alpha,C,eps){
    return(.onevsall(x=data,idx=y,alpha=alpha,C=C,eps=eps,exclude=FALSE,diagonal=TRUE))},
                   data=x,alpha=alpha,C=C,eps=eps)
  
  stopCluster(cl)
  Mat1 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat2 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat3 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat4 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  Mat5 <- matrix(0,ncol=f,nrow=f,dimnames=list(colnames(x),colnames(x)))
  
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
    
  }
  return(list(MIC=Mat1,MAS=Mat2,MEV=Mat3,MCN=Mat4,MICR2=Mat5))
}
