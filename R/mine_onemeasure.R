

# @param x Numeric Vector of size `n`
# @param y Numeric Vector of size `n`
# @param alpha an optional number of cells allowed in the X-by-Y search-grid. Default value is 0.6 (see Details of [mine()] function)
# @param C an optional number determining the starting point of the _X-by-Y_ search-grid. When trying to partition the _x_-axis into
# _X_ columns, the algorithm will start with at most `C` _X_
# _clumps_. Default value is 15 (see Details of [mine()] function).
# @param est Default value is "mic_approx". With `est="mic_approx"` the original MINE statistics will
# be computed, with `est="mic_e"` the equicharacteristic matrix is
# is evaluated and the `mic()` and `tic()` methods will return `MIC_e` and
# `TIC_e` values respectively.
# @param measure character indicating the measure to extract. Availiable values are: "mic", "mas", "mev", "mcn", "tic", "gmic".
# @param eps nteger in (0,1).  If `NULL` (default) it is set to 1-MIC.
# It can be set to zero for noiseless functions, but the default choice is the most appropriate parametrization
# for general cases (as stated in Reshef et al. SOM).
# It provides robustness.


#' Function to compute one measure at time
#' @aliases MIC mic MAS mas MEV mev MCN mcn GMIC gmic TIC tic
#' @param measure character indicating the measure to extract. Default value "mic". 
#' Availiable values are: "mic", "mas", "mev", "mcn", "tic", "gmic".
#' @param p probability for the generalized `mic`
#' @param norm boolean if require normalization between 0 and 1 for the `tic` statistic
#' @inheritParams mine
#' @examples 
#' x <- runif(10); y <- 3*x+2;
#' mine_stat(x,y, measure="mic")
#' @export
mine_stat <- function(x, y, alpha = 0.6, C = 15, est="mic_approx", measure="mic", eps=NULL, p=-1, norm=FALSE){
  mymeasures <- c(mic=1L, mas=2L, mev=3L, mcn=4L, tic=5L, gmic=6L)
  
  mm <- mymeasures[pmatch(measure, names(mymeasures))]
  if (is.na(mm))
    mm <- 1
  
  if (is.null(eps))
    eps <- NA
  
  mine_compute(x, y, alpha, C, est, mm, eps, p, norm)
}
