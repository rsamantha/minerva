\encoding{UTF-8}
\name{mine}
\alias{mine}
\alias{MINE}
\alias{MIC}
\alias{mic}
\alias{MAS}
\alias{mas}
\alias{MEV}
\alias{mev}
\alias{MCN}
\alias{mcn}
\alias{MIC-R2}
\alias{mic-r2}
\title{MINE family statistics}
\description{Maximal Information-Based Nonparametric Exploration (MINE) 
statistics. \code{mine} computes the MINE family measures between two variables. 
}
\usage{mine(x, y=NULL, master=NULL, alpha=0.6, C=15, n.cores=1,
var.thr=1e-5, eps=NULL, est="mic_approx", na.rm=FALSE, use="all.obs",  ...)}
\arguments{
  \item{x}{ a numeric vector (of size \emph{n}), matrix or data frame (which is coerced 
    to matrix). }
  \item{y}{  NULL (default) or a numeric vector of size \emph{n} (\emph{i.e.}, 
      with compatible dimensions to x). }
  \item{master}{ an optional vector of indices (numeric or character) to
    be given when \code{y} is not set, otherwise master is ignored. It can be 
    either one column index to be used as reference for the comparison 
    (versus all other columns) or a vector of column indices to be used for 
    computing all mutual statistics.
    If not specified it is set to \code{1:ncol(x)}. }
  \item{alpha}{ an optional number of cells allowed in the \emph{X-by-Y} 
    search-grid. Default value is 0.6 (see Details). }
  \item{C}{ an optional number determining the starting point of the 
    \emph{X-by-Y} search-grid. When trying to partition the \emph{x}-axis into 
    \emph{X} columns, the algorithm will start with at most \code{C}\emph{X} 
    \emph{clumps}. Default value is 15 (see Details). }
  \item{n.cores}{ optional number of cores to be used in the
    computations, when master is specified. It requires the
    \pkg{parallel} package, which provides support for parallel 
    computing, released with \R >= 2.14.0. Defaults is 1 (\emph{i.e.}, 
    not performing parallel computing). }
  \item{var.thr}{ minimum value allowed for the variance of the input
    variables, since \code{mine} can not be computed in case of variance
    close to 0. Default value is 1e-5. Information about failed check
    are reported in \emph{var_thr.log} file. }
  \item{eps}{ integer in [0,1].  If 'NULL' (default) it is set to
    1-MIC. It can be set to zero for noiseless functions,
    but the default choice is the most appropriate parametrization
    for general cases (as stated in Reshef et al. SOM).
    It provides robustness. }
  \item{est}{Default value is "mic_approx". With est="mic_approx" the original MINE statistics will
        be computed, with est="mic_e" the equicharacteristic matrix is
        is evaluated and the mic() and tic() methods will return MIC_e and
        TIC_e values respectively.}
  \item{na.rm}{boolean. This variable is passed directly to the
	\code{cor}-based functions. See \code{cor} for further details.}
  \item{use}{Default value is "all.obs". This variable is passed directly to the
    \code{cor}-based functions. See \code{cor} for further details.}
  
\item{\dots}{ currently not used. }
} 
\details{\code{mine} is an R wrapper for the C engine \emph{cmine} 
(\url{http://minepy.readthedocs.io/en/latest/}), 
an implementation of Maximal Information-Based Nonparametric Exploration (MINE) 
statistics. The MINE statistics were firstly detailed in 
D. Reshef et al. (2011) \emph{Detecting novel associations in large datasets}. 
Science 334, 6062 (\url{http://www.exploredata.net}).

Here we recall the main concepts of the MINE family statistics.
Let \eqn{D={(x,y)}} be the set of \emph{n} ordered pairs of elements of \code{x}
and \code{y}. The data space is partitioned in 
an \emph{X-by-Y} grid, grouping the \emph{x} and \emph{y} values 
in \emph{X} and \emph{Y} bins respectively.\cr

The \strong{Maximal Information Coefficient (MIC)} is defined as 
\deqn{\textrm{MIC}(D)=\max_{XY<B(n)} M(D)_{X,Y} = \max_{XY<B(n)} \frac{I^*(D,X,Y)}{log(\min{X,Y})},}{MIC(D)=max_{XY<B(n)} M(D)_{X,Y}=max_{XY<B(n)} I*(D,X,Y)/log(min(X,Y)),} where
\eqn{B(n)=n^{\alpha}} is the search-grid size,
\eqn{I^*(D,X,Y)}{I*(D,X,Y)}
is the maximum mutual information over all grids \emph{X-by-Y}, of the distribution induced by D on 
a grid having \emph{X} and \emph{Y} bins (where the probability mass on a cell 
of the grid is the fraction of points of D falling in that cell).
The other statistics of the MINE family are derived from the mutual information 
matrix achieved by an \emph{X-by-Y} grid on D.\cr

The \strong{Maximum Asymmetry Score (MAS)} is defined as
\deqn{\textrm{MAS}(D) = \max_{XY<B(n)} |M(D)_{X,Y} - M(D)_{Y,X}|.}{MAS(D) = max_{XY<B(n)} |M(D)_{X,Y} - M(D)_{Y,X}|.}

The \strong{Maximum Edge Value (MEV)} is defined as
\deqn{\textrm{MEV}(D) = \max_{XY<B(n)} \{M(D)_{X,Y}: X=2~or~Y=2\}.}{MEV(D) = max_{XY<B(n)} {M(D)_{X,Y}: X=2 or Y=2}.}

The \strong{Minimum Cell Number (MCN)} is defined as
\deqn{\textrm{MCN}(D,\epsilon) = \min_{XY<B(n)} \{\log(XY): M(D)_{X,Y} \geq (1-\epsilon)MIC(D)\}.}{MCN(D,\epsilon) = min_{XY<B(n)} {log(XY): M(D)_{X,Y} >= (1-\epsilon)MIC(D)}.}
More details are provided in the supplementary material (SOM) of the original paper.

The MINE statistics can be computed for two numeric vectors \code{x} and \code{y}. 
Otherwise a matrix (or data frame) can be provided and two options are available 
according to the value of \code{master}. If \code{master} is a column identifier, 
then the MINE statistics are computed for the \emph{master} variable versus the 
other matrix columns. If \code{master} is a set of column identifiers, then all 
mutual MINE statistics are computed among the column subset.
\code{master}, \code{alpha}, and \code{C} refers respectively to the \emph{style},
\emph{exp}, and \emph{c} parameters of the original \emph{java} code. 
In the original article, the authors state that the default value \eqn{\alpha=0.6} 
(which is the exponent of the search-grid size \eqn{B(n)=n^{\alpha}}) has been 
empirically chosen. It is worthwhile noting that \code{alpha} and \code{C} are 
defined to obtain an heuristic approximation in a reasonable amount of time. In case
of small sample size (\emph{n}) it is preferable to increase \code{alpha} to 1 to 
obtain a solution closer to the theoretical one.
}
\value{The Maximal Information-Based Nonparametric Exploration (MINE) statistics 
  provide quantitative evaluations of different aspects of the relationship 
  between two variables. 
  In particular \code{mine} returns a list of 5 statistics: 
  \item{MIC}{ \strong{Maximal Information Coefficient.} \cr 
    It is related to the relationship strenght and it can be interpreted as a 
    correlation measure. It is symmetric and it ranges in [0,1], where it 
    tends to 0 for statistically independent data and it approaches 1 in 
    probability for noiseless functional relationships (more details can 
    ben found in the original paper). }
  \item{MAS}{ \strong{Maximum Asymmetry Score.} \cr 
    It captures the deviation from monotonicity. Note that 
    \eqn{\textrm{MAS} < \textrm{MIC}}{MAS < MIC}. \cr
    \emph{Note:} it can be useful for detecting periodic relationships 
    (unknown frequencies). }
  \item{MEV}{ \strong{Maximum Edge Value.} \cr 
    It measures the closeness to being a function. Note that 
    \eqn{\textrm{MEV} \leq \textrm{MIC}}{MEV <= MIC}. }
  \item{MCN}{ \strong{Minimum Cell Number.} \cr 
    It is a complexity measure. }
  \item{MIC-R2}{It is the difference between the MIC value and the Pearson
    correlation coefficient. } \cr

  When computing \code{mine} between two numeric vectors \code{x} and \code{y}, 
  the output is a list of 5 numeric values. When \code{master} is provided, 
  \code{mine} returns a list of 5 matrices having \code{ncol} equal to 
  \emph{m}. In particular, if \code{master} is a single value, 
  then \code{mine} returns a list of 5 matrices having 1 column, 
  whose rows correspond to the MINE measures between the \emph{master} 
  column versus all. Instead if \code{master} is a vector of \emph{m} indices,
  then \code{mine} output is a list of 5 \emph{m-by-m} matrices, whose element 
  \emph{i,j} corresponds to the MINE statistics computed between the \emph{i} 
  and \emph{j} columns of \code{x}.
}
\references{
  D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P. Turnbaugh, 
  E. Lander, M. Mitzenmacher, P. Sabeti. (2011)
  \emph{Detecting novel associations in large datasets}. 
  Science 334, 6062\cr
  \url{http://www.exploredata.net}\cr
  (SOM: Supplementary Online Material at
  \url{http://www.sciencemag.org/content/suppl/2011/12/14/334.6062.1518.DC1})
  
  \emph{minepy. Maximal Information-based Nonparametric Exploration in C and Python.}\cr 
  \url{http://minepy.sourceforge.net}
  }
\author{Michele Filosi and Roberto Visintainer\cr
\cr
Special thanks to:\cr
Davide Albanese, Giuseppe Jurman, Samantha Riccadonna
}
\examples{

A <- matrix(runif(50),nrow=5)
mine(x=A, master=1)
mine(x=A, master=c(1,3,5,7,8:10))

x <- runif(10); y <- 3*x+2; plot(x,y,type="l")
mine(x,y)
# MIC = 1 
# MAS = 0
# MEV = 1
# MCN = 2
# MIC-R2 = 0

set.seed(100); x <- runif(10); y <- 3*x+2+rnorm(10,mean=2,sd=5); plot(x,y)
mine(x,y)
# rounded values of MINE statistics
# MIC = 0.61
# MAS = 0
# MEV = 0.61
# MCN = 2
# MIC-R2 = 0.13

t <-seq(-2*pi,2*pi,0.2); y1 <- sin(2*t); plot(t,y1,type="l")
mine(t,y1)
# rounded values of MINE statistics
# MIC = 0.66 
# MAS = 0.37
# MEV = 0.66
# MCN = 3.58
# MIC-R2 = 0.62

y2 <- sin(4*t); plot(t,y2,type="l")
mine(t,y2)
# rounded values of MINE statistics
# MIC = 0.32 
# MAS = 0.18
# MEV = 0.32
# MCN = 3.58
# MIC-R2 = 0.31

# Note that for small n it is better to increase alpha
mine(t,y1,alpha=1)
# rounded values of MINE statistics
# MIC = 1 
# MAS = 0.59
# MEV = 1
# MCN = 5.67
# MIC-R2 = 0.96

mine(t,y2,alpha=1)
# rounded values of MINE statistics
# MIC = 1 
# MAS = 0.59
# MEV = 1
# MCN = 5
# MIC-R2 = 0.99

# Some examples from SOM
x <- runif(n=1000, min=0, max=1)

# Linear relationship
y1 <- x; plot(x,y1,type="l"); mine(x,y1)
# MIC = 1 
# MAS = 0
# MEV = 1
# MCN = 4
# MIC-R2 = 0

# Parabolic relationship
y2 <- 4*(x-0.5)^2; plot(sort(x),y2[order(x)],type="l"); mine(x,y2)
# rounded values of MINE statistics
# MIC = 1 
# MAS = 0.68
# MEV = 1
# MCN = 5.5
# MIC-R2 = 1

# Sinusoidal relationship (varying frequency)
y3 <- sin(6*pi*x*(1+x)); plot(sort(x),y3[order(x)],type="l"); mine(x,y3)
# rounded values of MINE statistics
# MIC = 1 
# MAS = 0.85
# MEV = 1
# MCN = 4.6
# MIC-R2 = 0.96

# Circle relationship
t <- seq(from=0,to=2*pi,length.out=1000)
x4 <- cos(t); y4 <- sin(t); plot(x4, y4, type="l",asp=1)
mine(x4,y4)
# rounded values of MINE statistics
# MIC = 0.68 
# MAS = 0.01
# MEV = 0.32
# MCN = 5.98
# MIC-R2 = 0.68

data(Spellman)
res <- mine(Spellman,master=1,n.cores=1)

\dontrun{## example of multicore computation
res <- mine(Spellman,master=1,n.cores=parallel::detectCores()-1)}
}
