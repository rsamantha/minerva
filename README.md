# minerva


R package for Maximal Information-Based Nonparametric Exploration computation

* Minepy [Homepage](http://minepy.readthedocs.io/en/latest/)
* Minepy [Github](https://github.com/minepy/minepy)


## Install
* Latest cran release
```r
install.packages("minerva")
```
* Development version
```r
devtools::install_github('rsamantha/minerva')
```

## Usage

-  Basic usage with helper function `mine`.

```r
library(minerva)

x <- 0:200 / 200
y <- sin(10 * pi * x) + x
mine(x,y, n.cores=1)
```
-  Compute a single measure from the MINE suite using `mine_stat`.
   +  Available mesures are: *mic*, *mas*, *mev*, *mcn*, *tic*, *gmic*

```r
x <- 0:200 / 200
y <- sin(10 * pi * x) + x
mine_stat(x, y, measure="mic")
```
### Compute statistic on matrices

-  All features in a single matrix (`mine_compute_pstat`).
-  All possible combination of features between two matrices (`mine_compute_cstat`). 
   + When comparing two matrices the function check for euquality of number of rows of the two matrices. If the matrices have different number of rows then an error is thrown.

```r
x <- matrix(rnorm(1000), ncol=10, nrow=10)
y <- as.matrix(rnorm(1000), ncol=10, nrow=20)

## Compare feature of the same matrix
pstats(x)

## Compare features of matrix x with feature in matrix y
cstats(x, y)
```

----

## Mictools pipeline
This is inspired to the original implementation by Albanese et al. available in python here:
[https://github.com/minepy/mictools](https://github.com/minepy/mictools).


#### Reading the data from mictool repository
```r
datasaurus <- read.table("https://raw.githubusercontent.com/minepy/mictools/master/examples/datasaurus.txt", 
header=TRUE, row.names=1, as.is=TRUE, stringsAsFactors=FALSE)
datasaurus.m <- t(datasaurus)
```

### Compute null distribution for tic_e

Automatically compute:

-  `tic_e` null distribution based on permutations.
-  histogram of the distribution with cumulative distribution.
-  Observed values of `tic_e` for each pair of variable in `datasaurus`.
-  Observed distribution of `tic_e`.
-  P-value for each variable pair association.

```r
ticnull <- mictools(datasaurus.m, nperm=100000, seed=1234)

## Get the names of the named list
names(ticnull)
##[1]  "tic"      "nulldist" "obstic"   "obsdist"  "pval"
str(ticnull)

ticnull$nulldist
ticnull$obsdist
```

##### Null Distribution

| BinStart| BinEnd| NullCount| NullCumSum|
|--------:|------:|---------:|----------:|
|    0e+00|      0|         0|      1e+05|
|    1e-04|      0|         0|      1e+05|
|    2e-04|      0|         0|      1e+05|
|    3e-04|      0|         0|      1e+05|
|    4e-04|      0|         0|      1e+05|
|    5e-04|      0|         0|      1e+05|
| ...     | ...   |    ....  |   ....    |



##### Observed distribution

| BinStart| BinEnd| Count| CountCum|
|--------:|------:|-----:|--------:|
|    0e+00|      0|     0|      325|
|    1e-04|      0|     0|      325|
|    2e-04|      0|     0|      325|
|    3e-04|      0|     0|      325|
|    4e-04|      0|     0|      325|
|    5e-04|      0|     0|      325|
| ...     | ...   | .... | ....    |


Plot `tic_e` and pvalue distribution.

```r
hist(ticnull$tic)

hist(ticenull$pval, breaks=50, freq=FALSE)
```

Use `p.adjust.method` to use a different pvalue correction method, or use the `qvalue` package to use Storey's qvalue.

```r
## Correct pvalues using qvalue
qobj <- qvalue(ticnull$pval$pval)

## Add column in the pval data.frame
ticnull$pval$qvalue <- qobj$qvalue
ticnull$pval
```

Same table as above with the qvalue column added at the end.

|   pval| I1| I2|Var1   |Var2         | adj.P.Val| qvalue|
|------:|--:|--:|:------|:------------|---------:|------:|
| 0.5202|  1|  2|away_x |bullseye_x   |      0.95|      1|
| 0.9533|  1|  3|away_x |circle_x     |      0.99|      1|
| 0.0442|  1|  4|away_x |dino_x       |      0.52|      0|
| 0.6219|  1|  5|away_x |dots_x       |      0.95|      1|
| 0.8922|  1|  6|away_x |h_lines_x    |      0.98|      1|
| 0.3972|  1|  7|away_x |high_lines_x |      0.91|      1|
| ...   |...|...| ...   | ...         | ...      | ....  | 

### Strenght of the association (MIC)

```r
## Use columns of indexes and FDR adjusted pvalue 
micres <- mic_strength(datasaurus.m, ticnull$pval, pval.col=c(6, 2, 3))
```

| TicePval|  MIC| I1| I2|
|:--------|:----|:--|:--|
|   0.0457| 0.42|  2| 15|
|   0.0000| 0.63|  3| 16|
|   0.0196| 0.50|  5| 18|
|   0.0162| 0.36|  9| 22|
|   0.0000| 0.63| 10| 23|
|   0.0000| 0.57| 13| 26|
| ...     | ... | ...|...|


Association strength computed based on the `qvalue` adjusted pvalue

```r
## Use qvalue adjusted pvalue 
micresq <- mic_strength(datasaurus.m, ticnull$pval, pval.col=c("qvalue", "Var1", "Var2"))
```

| TicePval|  MIC|I1         |I2         |
|:--------|:----|:----------|:----------|
|   0.0401| 0.42|bullseye_x |bullseye_y |
|   0.0000| 0.63|circle_x   |circle_y   |
|   0.0172| 0.50|dots_x     |dots_y     |
|   0.0143| 0.36|slant_up_x |slant_up_y |
|   0.0000| 0.63|star_x     |star_y     |
|   0.0000| 0.57|x_shape_x  |x_shape_y  |
| ...     | ... | ...       |...        | 


----


## Citing minepy/minerva and mictools


 |||
 |:----------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
 | [minepy2013](https://doi.org/10.1093/bioinformatics/bts707)                     | Davide Albanese, Michele Filosi, Roberto Visintainer, Samantha Riccadonna, Giuseppe Jurman and Cesare Furlanello. *minerva and minepy:a C engine for the MINE suite and its R, Python and MATLAB wrappers*. Bioinformatics (2013) 29(3): 407-408 first published online December 14, 2012 |
| [mictools2018](https://academic.oup.com/gigascience/article/7/4/giy032/4958979) | Davide Albanese, Samantha Riccadonna, Claudio Donati, Pietro Franceschi. *A practical tool for maximal information coefficient analysis*. GigaScience (2018)                                                                                                                              |

