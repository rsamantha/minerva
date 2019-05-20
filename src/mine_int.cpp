#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include "mine.h"

/* DEFINE CONSTANT MAPS FOR MEASURE AND EST */
const std::map<std::string, int> MEASURE{{"mic", 1}, {"mas", 2}, {"mev", 3}, {"mcn", 4}, {"tic", 5}, {"gmic", 6}};
const std::map<std::string, int> EST{{"mic_approx", 0}, {"mic_e", 1}};

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

int switch_measure(String measure)
{
  int res=0;
  
  if (MEASURE.find(measure)!=MEASURE.end())
    res = MEASURE.find(measure)->second;
  
  return(res);
}

int switch_est(String est)
{
  int res = -1;
  
  if (EST.find(est)!=EST.end())
    res = EST.find(est)->second;
  
  return(res);
}

char *check_eps(double eps)
{
  if ((eps < 0.0) | (eps > 1.0))
    return "'eps' must be > 0.0 and < 1.0";
  
  return NULL;
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
//' @param eps eps value for MCN statistic should be in (0,1). If NA is passed then the normal MCN statistic is returned.
//' @param p probability for the generalized mic
//' @param norm boolean if require normalization between 0 and 1 for the tic statistic
//' @details This is a wrapper function to compute the mine statistic between two variables.
//' for more details on the available measure and the meaning of the other parameters see also the 
//' documentation for the \code{\link[minerva]{mine}} function.
//' @seealso \code{\link[minerva]{mine}}
//' @examples 
//' x <- runif(10); y <- 3*x+2;
//' mine_stat(x,y, measure="mic")
//' @export
// [[Rcpp::export]]
double mine_stat(NumericVector x, NumericVector y, double alpha=0.6, double C=15, String est="mic_approx", String measure="mic", double eps=0.0, double p=-1, bool norm=false)
{
  mine_problem *prob;
  mine_score *minescore;
  mine_parameter *param;
  char  *err;
  double res;
  int m, lest;

  /* Convert measure and est parameter */
  lest = switch_est(est);
  m = switch_measure(measure);
  
  /* Allocate parameter class */
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=alpha;
  param->c=C;
  param->est=lest;
  
  /* Check parameteres for MINE statistic */
  err = mine_check_parameter(param);
  if (err)
    stop(err);
  
  /* Check eps */
  err = check_eps(eps);
  if (err)
    stop(err);
  
  /* Check vector dimension compatibility */
  if (x.length() != y.length())
    stop("Not conformable arrays!");
  
  /* Allocate problem structure for MINE */
  prob = (mine_problem *) Calloc(1,mine_problem);
  prob->n=x.size();
  prob->x=REAL(x);
  prob->y=REAL(y);
  
  /* Compute MIC score */
  minescore = mine_compute_score(prob, param);
  
  /* Switch between measure to return */
  switch(m){
  case 1: res = mine_mic(minescore);
    break;
  case 2: res = mine_mas(minescore);
    break;
  case 3: res = mine_mev(minescore);
    break;
  case 4:
    if (!std::isnan(eps))
      res = mine_mcn(minescore, eps);
    else
      res = mine_mcn_general(minescore);
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
  Free(param);
  Free(prob);
  
  return(res);
}

//' @export
// [[Rcpp::export]]
NumericMatrix mine_allvar_onemeasure(NumericMatrix x, double alpha=0.6, double C=15, String est="mic_approx", String measure="mic", double eps=0.0, double p=-1, bool norm=false)
{

  int n, v;
  int i, j, it;
  double res;

  /* Set matrix dimension */
  n = x.nrow();
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
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=alpha;
  param->c=C;
  param->est=1; // est == "mic_e"
  
  /* Check Parameters */
  err = mine_check_parameter(param);
  if (err){
    stop(err);
  }
  
  /* Set up mine problem */
  prob = (mine_problem *) Calloc(1,mine_problem);
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
  Free(param);
  Free(prob);
  
  return(restic);
}


/* Helper funtion to convert from R Matrix to mine_matrix structure */
mine_matrix *rMattomine(NumericMatrix x)
{
  mine_matrix *X;
  
  X = (mine_matrix *) Calloc(1, mine_matrix);
  X->data=REAL(x);
  X->n=x.ncol();
  X->m=x.nrow();
  
  return(X);
}


//' @inheritParams mictools_null
//' @param est estimation parameter for the mine statistic
//' @describeIn mictools_null Computing the tic and mic statistics
//' @export
// [[Rcpp::export]]
NumericMatrix mictools_pstats(NumericMatrix x, double alpha=0.6, int C=15, String est="mic_approx")
{
  int i;
  mine_matrix *X;
  mine_parameter *param;
  mine_pstats *pstat;
  char *err;
  int EST=-1;
  
  if (est == "mic_approx")
    EST = 0;
  
  if (est == "mic_e")
    EST = 1;
  
  /* Set parameter */
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=alpha;
  param->c=C;
  param->est=EST;
  
  /* Check Error in parameters */
  err = mine_check_parameter(param);
  if (err){
    stop(err);
  }
  
  X = rMattomine(x);
  
  pstat = mine_compute_pstats(X, param);
  
  NumericMatrix stats(pstat->n, 2);
  for (i=0; i<pstat->n; i++)
  {
    stats(i,0) = pstat->mic[i];
    stats(i,1) = pstat->tic[i];
  }
  /* Free */
  Free(param);
  /* free(pstat); */
  
  return(stats);
}
