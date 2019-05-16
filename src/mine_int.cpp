#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>
#include "mine.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Function to compute one statistic at time
//' 
//' @param x Numeric Vector
//' @param y Numeric Vector
//' @param alpha alpha parameter for the mine statistic
//' @param C c parameter for the mine statistic
//' @param est estimation parameter for the mine statistic
//' @param measure character which measure to return
//' @param eps eps value for MCN statistic
//' @param p probability for the generalized mic
//' @param norm boolean if require normalization between 0 and 1 for the tic statistic
//' @export
// [[Rcpp::export]]
double mine_compute(NumericVector x, NumericVector y, double alpha=0.6, double C=15, String est="mic_approx", int measure=1, double eps=0.0, double p=-1, bool norm=false)
{
  mine_problem *prob;
  mine_score *minescore;
  mine_parameter *param;
  double res;
  
  int EST=0;
  if (est == "mic_e")
    {
      EST = 1;
    }
  
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=alpha;
  param->c=C;
  param->est=EST;
  
  prob = (mine_problem *) Calloc(1,mine_problem);
  prob->n=x.size();
  prob->x=REAL(x);
  prob->y=REAL(y);
  
  /* Compute MIC score */
  minescore = mine_compute_score(prob, param);
    
  switch(measure){
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
  default: res = NA_REAL;
  };
  
  /* Free score struct */
  mine_free_score(&minescore);
  Free(param);
  Free(prob);
  
  return(res);
}


void set_seed(unsigned int seed)
{
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

//' Function to compute mictools null distribution
//' 
//' @inheritParams mine_compute
//' @param nperm numper of permutation
//' @param seed character which measure to return
//' 
//' @export
// [[Rcpp::export]]
NumericVector mictools_null(NumericMatrix x, double alpha=9, double C=5, int nperm=250000, int seed=0)
{
  
  int i;
  int n = x.nrow();
  int v = x.ncol();

  IntegerVector idxs = seq(0, v-1);
  NumericVector restic(nperm);
  
  mine_problem *prob;
  mine_parameter *param;
  mine_score *minescore;
  
  // Set parameters for mine computation /
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=alpha;
  param->c=C;
  // est == "mic_e"
  param->est=1;
  
  // Set up mine problem //
  prob = (mine_problem *) Calloc(1,mine_problem);
  prob->n=n;
  
  // set seed //
  set_seed(seed);
  
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
  
  return restic;
}


/* Function needed to convert from R Matrix to mine_matrix structure */
mine_matrix *rMattomine(NumericMatrix x)
{
  mine_matrix *X;
  
  X = (mine_matrix *) Calloc(1, mine_matrix);
  X->data=REAL(x);
  X->n=x.ncol();
  X->m=x.nrow();
  
  return(X);
}


//' Function to compute the pvalue for the mictools pipeline
//' 
//' @inheritParams mictools_null
//' @param est estimation parameter for the mine statistic
//' @export
// [[Rcpp::export]]
NumericMatrix mictools_pstats(NumericMatrix x, double alpha=0.6, int C=15, String est="mic_approx")
{
  int i;
  mine_matrix *X;
  mine_parameter *param;
  mine_pstats *pstat;
  int EST=0;
  
  if (est == "mic_e")
  {
    EST = 1;
  }
  
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=alpha;
  param->c=C;
  param->est=EST;
  
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
  /* free(pstat->tic); */
  /* free(pstat->mic); */
  /* free(pstat); */
  
  return(stats);
}


/*
//' Function to compute mic strength based on pvalue
//' 
//' @inheritParams mictools_null
//' @export
// [[Rcpp::export]]
 */
/* NumericMatrix mictools_strength(NumericMatrix x, double alpha=0.6, int C=5, String est = "mic_approx")
{
  NumericVector binpoints ={1, 25, 50, 250, 500, 1000, 2500, 5000, 10000, 40000};
  NumericVector alphas = {0.85, 0.80, 0.75, 0.70, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4};
  
} */