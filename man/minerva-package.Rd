\encoding{UTF-8}
\name{minerva-package}
\alias{minerva-package}
\alias{minerva}
\docType{package}
\title{The minerva package}
\description{
  Maximal Information-Based Nonparametric Exploration R Package for 
  Variable Analysis. The package provides the \code{\link{mine}} function
  allowing the computation of Maximal Information-based Nonparametric
  Exploration statistics, firstly introduced in D. Reshef et al. (2011)
  \emph{Detecting novel associations in large datasets}. Science 334,
  6062 (\url{http://www.exploredata.net}). In particular, the package
  is an R wrapper for the C engine \emph{cmine}
  (\url{http://minepy.readthedocs.io/en/latest/}).  }
\details{
Summary:
  \tabular{ll}{
  Package: \tab minerva\cr
  Version: \tab 1.4.3\cr
  Date: \tab 2014-10-08\cr
  Depends: \tab R >= (2.14.0)\cr
  Enshances: \tab parallel\cr
  URL: \tab http://www.r-project.org, \cr
  \tab http://minepy.readthedocs.io/en/latest/, \cr
  \tab http://www.exploredata.net\cr
  License: \tab GPL-3\cr
  }
  Index:
  \tabular{ll}{
  \code{\link{Spellman}} \tab Yeast Gene Expression Dataset\cr
  \code{\link{mine}} \tab MINE-family statistics\cr
  \code{\link{minerva-package}} \tab The minerva package\cr
  }
}
\author{
  Michele Filosi [aut, cre],
  Roberto Visintainer [aut],
  Davide Albanese [aut],
  Samantha Riccadonna [ctb],
  Giuseppe Jurman [ctb],
  Cesare Furlanello [ctb]
  
  Maintainer: Michele Filosi <filosi@fbk.eu>

}
\references{
  D. Reshef, Y. Reshef, H. Finucane, S. Grossman, G. McVean, P. Turnbaugh, 
  E. Lander, M. Mitzenmacher, P. Sabeti. (2011)
  \emph{Detecting novel associations in large datasets}. 
  Science 334, 6062  
  (\href{http://www.exploredata.net}{http://www.exploredata.net}).

  D. Albanese, M. Filosi, R. Visintainer, S. Riccadonna, G. Jurman,
  C. Furlanello. \emph{cmine, minerva & minepy: a C engine for the
    MINE suite an its R and Python wrappers.}
  \url{http://minepy.readthedocs.io/en/latest/}  

  \emph{minepy. Maximal Information-based Nonparametric Exploration
    in C and Python.}\cr 
  (\url{http://minepy.sourceforge.net})
}
\keyword{ package }
