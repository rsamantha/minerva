context("Check MINE scores consistency")

## Create constant vector
const.create <- function(n){
    x <- seq(0, 1, length.out=n)
    y <- rep(0, n)
    retval <- data.frame(X=x, Y=y)
    return(retval)
}

## Create two linear vector
lin.create <- function(n){
    x <- seq(0, 1, length.out=n)
    retval <- data.frame(X=x, Y=x)
    return(retval)
}


## Create  SIN
sin.create <- function(n){
    x <- seq(0, 1, length.out=n)
    y <- sin(8 * pi * x)
    retval <- data.frame(X=x, Y=y)
    return(retval)
}

## Create exp
exp.create <- function(n){
    x <- seq(0, 1, length.out=n)
    y <- 2**x
    retval <- data.frame(X=x, Y=y)
    return(retval)
}

test_that("Test constant:", {
    mydata <- const.create(1000)
    expect_warning(mm <- mine(mydata$X, mydata$Y))
    expect_equal(mm$MIC, 0, tolerance=1e-4)
    expect_equal(mm$MAS, 0, tolerance=1e-4)
    expect_equal(mm$MEV, 0, tolerance=1e-4)
    expect_equal(mm$MCN, 2, tolerance=1e-4)
}
)

test_that("Test linear:", {
    mydata <- lin.create(1000)
    mm <- mine(mydata$X, mydata$Y)
    expect_equal(mm$MIC, 1, tolerance=1e-4)
    expect_equal(mm$MAS, 0, tolerance=1e-4)
    expect_equal(mm$MEV, 1, tolerance=1e-4)
    expect_equal(mm$MCN, 2, tolerance=1e-4)
}
)

test_that("Test sin:", {
    mydata <- sin.create(1000)
    mm <- mine(mydata$X, mydata$Y)
    expect_equal(mm$MIC, 1, tolerance=1e-4)
    expect_equal(mm$MAS, 0.875, tolerance=1e-4)
    expect_equal(mm$MEV, 1, tolerance=1e-4)
    expect_equal(mm$MCN, 4, tolerance=1e-4)
}
)

test_that("Test exp:", {
    mydata <- exp.create(1000)
    mm <- mine(mydata$X, mydata$Y)
    expect_equal(mm$MIC, 1, tolerance=1e-4)
    expect_equal(mm$MAS, 0, tolerance=1e-4)
    expect_equal(mm$MEV, 1, tolerance=1e-4)
    expect_equal(mm$MCN, 2, tolerance=1e-4)
}
)

## test_that("Test bug 1:", {
##     a <- mine(x = cbind(matrix(c(rnorm(6), 0, 0, 0), nrow = 3), rnorm(3)), master = 4)
##     b <- mine(x = matrix(c(rnorm(6), 0, 0, 0), nrow = 3), y = rnorm(3))
##     expect_equal(a$MIC, b$MIC, tolerance=1e-4)
## }
## )


context("Check Mictools pipeline")

test_that("Mictools pipeline:", {
  mydata <- matrix(rep(1:10, 4), ncol=4, nrow=10) + matrix(rep((0:3), each=10), ncol=4, nrow=10)
  nperm <- 100
  rr <- mictools(mydata,nperm=nperm, seed=0)
  lr <- length(rr)
  expect_equal(lr, 5, tolerance=1e-4)
  expect_equal(length(rr$tic), nperm)
})

test_that("mic_strength pval.col parameters",{
  mydata <- matrix(rep(1:10, 4), ncol=4, nrow=10) + matrix(rep((0:3), each=10), ncol=4, nrow=10)
  rr <- mictools(mydata, nperm=100, seed=10)
  expect_error(mic_strength(mydata, pval=rr$pval, pval.col=6))
})

test_that("mic_strength high correlation == 1",{
  mydata <- matrix(rep(1:10, 4), ncol=4, nrow=10) + matrix(rep((0:3), each=10), ncol=4, nrow=10)
  rr <- mictools(mydata, nperm=100, seed=10)
  ms <- mic_strength(mydata, pval=rr$pval, pval.col=c(6, 2, 3))
  expect_equal(ms$MIC-1, rep(0, nrow(ms)), tolerance=1e-10)
})


context("Rcpp Interface")

test_that("Test rcpp interface:", {
  mydata <- lin.create(1000)
  mm <- mine(mydata$X, mydata$Y, alpha=0.6, C=15)
  mm2 <- mine_stat(mydata$X, mydata$Y, alpha=0.6, C=15)
  expect_equal(mm$MIC, mm2, tolerance=1e-6)
})

test_that("Parameters error through Rcpp interface",{
  mydata <- lin.create(1000)
  expect_error(mine_stat(mydata[, 1], mydata[, 2], alpha=-5), NULL)
  expect_error(mine_stat(mydata[, 1], mydata[, 2], C=-3), NULL)
  expect_error(mine_stat(mydata[, 1], mydata[, 2], est='ciccio'), NULL)
  expect_error(mine_stat(mydata[, 1], mydata[, 2], eps=10, measure="mcn"))
})

test_that("Array dimension in cstats", {
  x <- matrix(rnorm(320*8), nrow=320, ncol=8)
  y <- matrix(rnorm(320*4), nrow=320, ncol=4)
  y2 <- matrix(rnorm(240*8), nrow=240, ncol=8)
  mictic <- cstats(x, y, alpha=9, C=5, est="mic_e")
  expect_error(cstats(x, y2, alpha=9, C=5, est="mic_e"))
  expect_equal(ncol(mictic), 4)
  expect_equal(nrow(mictic), 8*4)
})

