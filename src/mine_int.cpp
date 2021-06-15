#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include <limits>
#include <map>
#include "mine.h"
#include "mine_int.h"

/* DEFINE CONSTANT MAPS FOR MEASURE AND EST */
//const std::map<std::string, int> MEASURE=create_measure_map();
//({{"mic", 1}, {"mas", 2}, {"mev", 3}, {"mcn", 4}, {"tic", 5}, {"gmic", 6}});
//const std::map<std::string, int> EST=create_est_map();
//({{"mic_approx", 0}, {"mic_e", 1}});


const std::map<std::string, int> create_measure_map()
{
  std::map<std::string, int> MEASURE;
  // {{"mic", 1}, {"mas", 2}, {"mev", 3}, {"mcn", 4}, {"tic", 5}, {"gmic", 6}}
  MEASURE["mic"]=1;
  MEASURE["mas"]=2;
  MEASURE["mev"]=3;
  MEASURE["mcn"]=4;
  MEASURE["tic"]=5;
  MEASURE["gmic"]=6;
  MEASURE["MIC"]=1;
  MEASURE["MAS"]=2;
  MEASURE["MEV"]=3;
  MEASURE["MCN"]=4;
  MEASURE["TIC"]=5;
  MEASURE["GMIC"]=6;
  
  return(MEASURE);
}

const std::map<std::string, int> create_est_map()
{
  std::map<std::string, int> EST;
  EST["mic_approx"]=0;
  EST["mic_e"]=1;
  
  return(EST);
}

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/* Switch between available measure :
 * mic; mas; mev; mcn; tic; gmic
 */
int switch_measure(String measure)
{
  int res=0;
  
  if (MEASURE.find(measure)!=MEASURE.end())
    res = MEASURE.find(measure)->second;
  
  return(res);
}

/* Switch between possible estimators 
 * est = 0 -> est = 'mic_approx'
 * est = 1 -> est = 'mic_e'
 */
int switch_est(String est)
{
  int res = -1;
  
  if (EST.find(est)!=EST.end())
    res = EST.find(est)->second;
  
  return(res);
}

/* Handle eps parameter 
 * if eps not in (0, 1) throws an error...
 */
char *check_eps(double eps)
{
  if ((eps < 0.0) | (eps > 1.0))
    return (char *)"'eps' must be > 0.0 and < 1.0";
  
  return NULL;
}

/* Helper funtion to convert from R Matrix to mine_matrix structure
 * Not needed anymore. All the computation handled by pointers 
 * between R memory and C memory
 */ 
mine_matrix *rMattomine(NumericMatrix x)
{
  mine_matrix *X;
  
  X = (mine_matrix *) R_Calloc(1, mine_matrix);
  X->data=REAL(x);
  X->n=x.ncol();
  X->m=x.nrow();
  
  return(X);
}

//' This is an helper function to compute one \code{mine} statistic.
//' It take two vectors of the same dimension as an input.
//' 
//' @param x Numeric Vector of size \code{n}
//' @param y Numeric Vector of size \code{n}
//' @param alpha numeric value representing parameter for the mine statistic see \code{\link[minerva]{mine}}
//' @param C c parameter for the mine statistic see \code{\link[minerva]{mine}}
//' @param est character estimation parameter for the mine statistic.
//' Possible values are \code{"mic_approx"} or \code{"mic_e"}
//' @param measure integer indicating which measure to return
//' available measures are: \code{mic, mas, mev, mcn, tic, gmic}. The string could be also uppercase.
//' For measure \code{mic-r2} see details.
//' @param eps eps value for MCN statistic should be in (0,1). If NA (default) is passed then the normal MCN statistic is returned.
//' @param p probability for the generalized mic
//' @param norm boolean if require normalization between 0 and 1 for the \code{tic} statistic
//' @details This is a wrapper function to compute the mine statistic between two variables.
//' for more details on the available measure and the meaning of the other parameters see also the 
//' documentation for the \code{\link[minerva]{mine}} function.
//'
//' For measure \code{mic-r2} use the Pearson R coefficient score \code{\link[stats]{cor}} and the measure \code{mic}. 
//' See the example below.
//' @seealso \code{\link[minerva]{mine}}
//' @examples 
//' x <- runif(10); y <- 3*x+2;
//' mine_stat(x,y, measure="mic")
//' 
//' ## Measure mic-r2
//' x <- matrix(rnorm(20), ncol=2, nrow=10)
//' mmic <- mine_stat(x[,1], x[,2], measure="mic")
//' r2 <- cor(x[,1], x[,2])
//' 
//' mmic - r2**2
//' 
//' @export
// [[Rcpp::export]]
double mine_stat(NumericVector x, NumericVector y, double alpha =0.6, double C=15, String est="mic_approx", String measure="mic", double eps=NA_REAL, double p=-1, bool norm=false)
{
  
  mine_problem prob;
  mine_score *minescore;
  mine_parameter param;
  char  *err;
  double res;
  int m, lest;
  // int nnorm = norm==true;
  
  /* Convert measure and est parameter */
  lest = switch_est(est);
  m = switch_measure(measure);
  
  /* Allocate parameter class */
  // param = (mine_parameter *) Calloc(1,mine_parameter);
  param.alpha=alpha;
  param.c=C;
  param.est=lest;
  
  /* Check parameteres for MINE statistic */
  err = mine_check_parameter(&param);
  if (err)
    stop(err);
  
  /* Check vector dimension compatibility */
  if (x.length() != y.length())
    stop("Not conformable arrays!");
  
  /* Allocate problem structure for MINE */
  // prob = (mine_problem *) Calloc(1,mine_problem);
  prob.n=x.size();
  prob.x=x.begin();
  prob.y=y.begin();
  
  /* Compute MIC score */
  minescore = mine_compute_score(&prob, &param);
  
  /* Check eps for MCN */
  err = check_eps(eps);
  
  /* Switch between measure to return */
  switch(m){
  case 1: res = mine_mic(minescore);
    break;
  case 2: res = mine_mas(minescore);
    break;
  case 3: res = mine_mev(minescore);
    break;
  case 4:
    if (err){
      stop(err);
    } else {
      if (!std::isnan(eps))
        res = mine_mcn(minescore, eps);
      else
        res = mine_mcn_general(minescore);
    }
    break;
  case 5: res = mine_tic(minescore, norm);
    break;
  case 6: res = mine_gmic(minescore, p);
    break;
  case 0: res = NA_REAL;
    break;
  default: res = NA_REAL;
  };
  
  /* Free memory */
  mine_free_score(&minescore);

  return(res);
}


//' Compute pairwise statistics (MIC and normalized TIC) between variables
//' (convenience function).
//'
//' For each statistic, the upper triangle of the matrix is stored by row
//' (condensed matrix). If m is the number of variables, then for i < j < m, the
//' statistic between (col) i and j is stored in k = m*i - i*(i+1)/2 - i - 1 + j.
//' The length of the vectors is n = m*(m-1)/2.
//' @inheritParams cstats
//' @return
//' A matrix of (n x (n-1)/2) rows and 4 columns. The first and second column are
//' the indexes relative to the columns in the input matrix \code{x} for which the statistic is computed for.
//' Column 3 contains the MIC statistic, while column 4 contains the normalized TIC statistic. 
//' @examples
//' ## Create a matrix of random numbers
//' ## 10 variables x 100 samples
//' x <- matrix(rnorm(1000), ncol=10)
//' res <- pstats(x)
//' 
//' head(res)
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix pstats(NumericMatrix x, double alpha=0.6, double C=15, String est="mic_approx")
{
  mine_parameter param;
  mine_matrix mmat;
  mine_pstats *pstatscore;
  int i, j;
  int k=0;
  int nr=x.nrow();
  int nc=x.ncol();
  int lest;
  char *err;

  /* Set EST parameter */
  lest = switch_est(est);
  
  /* Allocate parameters for computation */
  // param = (mine_parameter *) Calloc(1,mine_parameter);
  param.alpha=alpha;
  param.c=C; 
  param.est=lest;
  
  /* Check parameteres for MINE statistic */
  err = mine_check_parameter(&param);
  if (err)
    stop(err);
  
  /*Prepare data structure */
  // mmat = (mine_matrix *) Calloc (1, mine_matrix);
  mmat.data = x.begin();
  mmat.m = nr;
  mmat.n = nc;
  
  pstatscore = mine_compute_pstats(&mmat, &param);
  
  NumericMatrix mres(pstatscore->n, 4);
  for (i=0; i<pstatscore->n; i++)
  {
    // Rcout << pstats->mic[i] << std::endl; 
    mres(i,2) = pstatscore->mic[i];
    mres(i,3) = pstatscore->tic[i];
  }
  
  /* Fill index of variable computed 1-based*/
  for (i=0; i<mmat.n-1; i++)
  {
    for (j=i+1; j<mmat.n; j++)
    {
      mres(k, 0) = i+1;
      mres(k, 1) = j+1;
      k++;
    }
  }
  
  /* Set colnames of the matrix */
  StringVector resnames = StringVector::create("Var1", "Var2", "MIC", "TIC");
  colnames(mres) = resnames;
  
  return(mres);
}

//' Compute statistics (MIC and normalized TIC) between each pair of the two
//' collections of variables (convenience function).
//' If n and m are the number of variables in X and Y respectively, then the
//'   statistic between the (row) i (for X) and j (for Y) is stored in \code{mic[i, j]}
//' and \code{tic[i, j]}.
//' @param x Numeric Matrix of m-by-n with n variables and m samples.
//' @param y Numeric Matrix of m-by-p with p variables and m samples.
//' @param alpha number (0, 1.0] or >=4 if alpha is in (0,1] then B will be max(n^alpha, 4) where n is the
//' number of samples. If alpha is >=4 then alpha defines directly the B
//' parameter. If alpha is higher than the number of samples (n) it will be
//' limited to be n, so B = min(alpha, n).
//' @param C number (> 0) determines how many more clumps there will be than columns in
//' every partition. Default value is 15, meaning that when trying to
//' draw x grid lines on the x-axis, the algorithm will start with at
//' most 15*x clumps.
//' @param est string ("mic_approx", "mic_e") estimator. 
//' With est="mic_approx" the original MINE statistics will
//' be computed, with est="mic_e" the equicharacteristic matrix is
//' is evaluated and MIC_e and TIC_e are returned.
//' 
//' @return list of two elements:
//' MIC: the MIC statistic matrix (n x p).
//' TIC: the normalized TIC statistic matrix (n x p).
//' @examples
//' x <- matrix(rnorm(2560), ncol=8, nrow=320)
//' y <- matrix(rnorm(1280), ncol=4, nrow=320)
//' 
//' mictic <- cstats(x, y, alpha=9, C=5, est="mic_e")
//' head(mictic)
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix cstats(NumericMatrix x, NumericMatrix y, double alpha=0.6, double C=15, String est="mic_approx")
{
  mine_parameter param;
  mine_matrix x_, y_;
  mine_cstats *cres;
  char *err;
  int i, j, l;
  int lest;
    
  /* Set EST parameter */
  lest = switch_est(est);

  /* Allocate parameters for computation */
  param.alpha=alpha;
  param.c=C;
  param.est=lest;
  
  /* Check parameteres for MINE statistic */
  err = mine_check_parameter(&param);
  if (err)
    stop(err);
 
  /*Prepare data structure */
  x_.data = x.begin();
  x_.m = x.nrow();
  x_.n = x.ncol();
  
  y_.data = y.begin();
  y_.m = y.nrow();
  y_.n = y.ncol();
  
  /* Compute pairwise statistic */
  cres = mine_compute_cstats(&x_, &y_, &param);
  
  if (!cres)
    stop("Not conformable arrays");
  
  /* Matrix for results */
  NumericMatrix rres(x_.n * y_.n, 4);
  
  /* Put result in the result matrix to pass to R */
  l = (cres->n * cres->m);
  for (i=0; i<l; i++)
    {
      rres(i, 2) = cres->mic[i];
      rres(i, 3) = cres->tic[i];
    }
   
  /* Fill Indices */
  l = 0;
  for (i=0; i<cres->n; i++)
  {
    for(j=0; j<cres->m; j++)
      {
        rres(l, 0) = i+1; 
        rres(l, 1) = j+1;
        l++;
      }
  }
  /* Set colnames of the matrix */
  StringVector resnames = StringVector::create("VarX", "VarY", "MIC", "TIC");
  colnames(rres) = resnames;
  
  return(rres);
}


/* 
 Do we need this function? 
 guess not...you can use mine
*/
NumericMatrix mine_allvar_onemeasure(NumericMatrix x, double alpha=0.6, double C=15, String est="mic_approx", String measure="mic", double eps=0.0, double p=-1, bool norm=false)
{

  int v;
  int i, j, it;
  double res;

  /* Set matrix dimension */
  v = x.ncol();
  
  NumericMatrix ret((v * (v-1))/2, 3);
  
  it=0;
  for (i=0; i<(v-1); i++)
  {
    for (j=i+1; j<v; j++)
    {
      /* Compute score */
      res = mine_stat(x(_,i), x(_, j), alpha, C, est, measure, eps, p, norm);
      
      /* Store results */
      ret(it, 0) = res;
      ret(it, 1) = i + 1;
      ret(it, 2) = j + 1;
      
      /* Increment iterator */
      it++;
    }
  }
  
  return(ret);
}


/* This function is needed to set the random seed */
void set_seed(unsigned int seed)
{
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

//' This set of functions are helper function to compute null distribution of the \code{tic_e} and 
//' \code{tic_e} observed distribution from a matrix
//'
//' @inheritParams mine_stat
//' @param x matrix N x M with M variables and N samples
//' @param nperm numper of permutation
//' @param seed integer to set the starting seed for random number generation (for reproducibility).
//' @describeIn mictools_null compute the \code{tic_e} null distribution
//' @return It returns a vector of \code{nperm} \code{tic_e} values.
//' @seealso \code{\link[minerva]{mictools}}
//' @export
// [[Rcpp::export]]
NumericVector mictools_null(NumericMatrix x, double alpha=9, double C=5, int nperm=200000, int seed=0)
{
  
  int i;
  int n = x.nrow();
  int v = x.ncol();
  char *err;

  IntegerVector idxs = seq(0, v-1);
  NumericVector restic(nperm);
  
  mine_problem *prob;
  mine_parameter *param;
  mine_score *minescore;
  
  /* Set parameters for mine computation */
  param = (mine_parameter *) R_Calloc(1,mine_parameter);
  param->alpha=alpha;
  param->c=C;
  param->est=1; // est == "mic_e"
  
  /* Check Parameters */
  err = mine_check_parameter(param);
  if (err){
    stop(err);
  }
  
  /* Set up mine problem */
  prob = (mine_problem *) R_Calloc(1,mine_problem);
  prob->n=n;
  
  /* set seed */
  set_seed(seed);
  
  /* Run permutation */
  for (i=0; i<nperm; i++)
  {
    /* Sample variables */
    IntegerVector ii = RcppArmadillo::sample(idxs, 2, false);
    
    /* Extract variables and permute one */
    NumericVector x1 = x(_,ii[0]);
    NumericVector x2 = x(_,ii[1]);
    x2 = RcppArmadillo::sample(x2, n, false);
    
    /* Set problem structure */
    prob->x=REAL(x1);
    prob->y=REAL(x2);
    
    /* Compute score */
    minescore = mine_compute_score(prob, param);
    restic[i] = mine_tic(minescore, true);
    
    /* Free score struct */
    mine_free_score(&minescore);
  }
  
  /* Free memory */
  R_Free(param);
  R_Free(prob);
  
  return(restic);
}



double corC(NumericVector x, NumericVector y)
{
  int i;
  int n = x.size();
  double xm=std::accumulate(x.begin(), x.end(), 0.0) / n;
  double ym=std::accumulate(y.begin(), y.end(), 0.0) / n;
  double num=0.0;
  
  double sx, sy, rr;
  sx=0.0;
  sy=0.0;
  for (i=0; i<n; i++)
  {

    num += (x[i] - xm) * (y[i] - ym);
    sx += pow((x[i] - xm), 2);
    sy += pow((y[i] - ym), 2);
  }
  
  rr =  num  / sqrt(sx * sy);
  return(rr);
}




