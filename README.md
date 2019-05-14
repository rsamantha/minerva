# minerva


R package for Maximal Information-Based Nonparametric Exploration computation

* Minepy [Homepage](http://minepy.readthedocs.io/en/latest/)
* Minepy [Github](https://github.com/minepy/minepy)


## Install
* Latest cran release
```
install.packages("minerva")
```
* Development version
```
devtools::install_github('rsamantha/minerva')
```

## Usage

-  Basic usage with helper function `mine`.

```
library(minerva)

x <- 0:200 / 200
y <- sin(10 * pi * x) + x
mine(x,y, n.cores=1)
```
-  Compute a single measure from the MINE suite using `mine_stat`.
   +  Available mesures are: *mic*, *mas*, *mev*, *mcn*, *tic*, *gmic*

```
x <- 0:200 / 200
y <- sin(10 * pi * x) + x
mine_stat(x, y, measure="mic")
```
### Compute statistic on matrices

-  All features in a single matrix (`mine_compute_pstat`).
-  All possible combination of features between two matrices (`mine_compute_cstat`). 
   + When comparing two matrices the function check for euquality of number of rows of the two matrices. If the matrices have different number of rows then an error is thrown.

```
x <- as.matrix(rnorm(1000), ncol=10)
y <- as.matrix(rnorm(1000), ncol=10)

## Compare feature of the same matrix
mine_compute_pstat(x)

## Compare features of matrix x with feature in matrix y
mine_compute_cstat(x, y)
```

## Tests

```
devtools::test()
```

## Citing minepy/minerva 
Davide Albanese, Michele Filosi, Roberto Visintainer, Samantha
Riccadonna, Giuseppe Jurman and Cesare Furlanello. minerva and minepy:
a C engine for the MINE suite and its R, Python and MATLAB
wrappers. Bioinformatics (2013) 29(3): 407-408 first published online
December 14, 2012 [doi:10.1093/bioinformatics/bts707](doi:10.1093/bioinformatics/bts707).

