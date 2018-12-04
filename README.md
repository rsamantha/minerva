# minerva


R package for Maximal Information-Based Nonparametric Exploration computation

* [Homepage](http://minepy.readthedocs.io/en/latest/)
* [Github](https://github.com/minepy/minepy)


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
```
library(minerva)

x <- 0:200 / 200
y <- sin(10 * pi * x) + x
mine(x,y, n.cores=1)
```
A new function has been added to extract a single statistic from the MINE suite.

```
x <- 0:200 / 200
y <- sin(10 * pi * x) + x
mine_stat(x, y)
```


## Citing minepy/minerva 
Davide Albanese, Michele Filosi, Roberto Visintainer, Samantha
Riccadonna, Giuseppe Jurman and Cesare Furlanello. minerva and minepy:
a C engine for the MINE suite and its R, Python and MATLAB
wrappers. Bioinformatics (2013) 29(3): 407-408 first published online
December 14, 2012 doi:10.1093/bioinformatics/bts707.

