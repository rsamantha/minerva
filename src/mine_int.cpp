#include <Rcpp.h>
#include "mine.h"

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
  
  /* Check size consistency of input vectors */
  if (x.size() != y.size())
    return(NA_REAL);
  
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

  return(res);
}


//' Compute pairwise statistics (MIC and normalized TIC) between variables
//' (convenience function).
//'
//' For each statistic, the upper triangle of the matrix is stored by row
//' (condensed matrix). If m is the number of variables, then for i < j < m, the
//' statistic between (col) i and j is stored in k = m*i - i*(i+1)/2 - i - 1 + j.
//' The length of the vectors is n = m*(m-1)/2.
//' @param x Numeric matrix of m-by-n of n variables and m samples
//' @param alpha alpha parameter for the mine statistic
//' @param C c parameter for the mine statistic
//' @param est estimation parameter for the mine statistic
//' 
//' @return
//' Matrix (n x (n-1)/2) by 4. The first and second column indicate the indexes relative of the columns 
//' in the inut matrix the statistic is computed for.
//' Column 3 contains the MIC statistic, while column 4 contains the normalized TIC statistic.
//' @export
// [[Rcpp::export]]
NumericMatrix mine_compute_pstats(NumericMatrix x, double alpha=0.6, double C=15, String est="mic_approx")
{
  mine_parameter param;
  mine_matrix mmat;
  mine_pstats *pstats;
  int i, j;
  int k=0;
  int nr=x.nrow();
  int nc=x.ncol();

  int EST=0;
  if (est == "mic_e")
  {
    EST = 1;
  }
  
  /* Allocate parameters for computation */
  param.alpha=alpha;
  param.c=C; 
  param.est=EST;
  
  /*Prepare data structure */
  mmat.data = x.begin();
  mmat.m = nr;
  mmat.n = nc;
  
  pstats = mine_compute_pstats(&mmat, &param);
  
  NumericMatrix mres(pstats->n, 4);
  for (i=0; i<pstats->n; i++)
  {
    // Rcout << pstats->mic[i] << std::endl; 
    mres(i,2) = pstats->mic[i];
    mres(i,3) = pstats->tic[i];
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
  
  /* Free used memory */
  free(pstats->mic);
  free(pstats->tic);
  free(pstats);
  
  return(mres);
}


//'Compute statistics (MIC and normalized TIC) between each pair of the two
//' collections of variables (convenience function).
//' If n and m are the number of variables in X and Y respectively, then the
//'   statistic between the (row) i (for X) and j (for Y) is stored in mic[i, j]
//' and tic[i, j].
//'  
//' @param x Numeric Matrix of m-by-n with n variables and m samples.
//' @param y Numeric Matrix of m-by-p with p variables and m samples.
//' @param alpha float (0, 1.0] or >=4 if alpha is in (0,1] then B will be max(n^alpha, 4) where n is the
//' number of samples. If alpha is >=4 then alpha defines directly the B
//' parameter. If alpha is higher than the number of samples (n) it will be
//' limited to be n, so B = min(alpha, n).
//' @param c float (> 0) determines how many more clumps there will be than columns in
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
//' @export
// [[Rcpp::export]]
NumericMatrix mine_compute_cstats(NumericMatrix x, NumericMatrix y, double alpha=0.6, double C=15, String est="mic_approx")
{
  mine_parameter param;
  mine_matrix x_, y_;
  mine_cstats *cres;
  int i, j, l;

  int EST=0;
  if (est == "mic_e")
  {
    EST = 1;
  }
  
  /* Allocate parameters for computation */
  param.alpha=alpha;
  param.c=C;
  param.est=EST;
  
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
  {
    free(cres);
    return(NA_REAL);
  }
  
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
  
  /* Free */
  free(cres->mic);
  free(cres->tic);
  free(cres);
  return(rres);
}