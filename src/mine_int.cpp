#include <Rcpp.h>
#include "mine.h"
// #include <stdlib.h>
// #include <stdio.h>
// #include <math.h>


using namespace Rcpp;

// enum mes_value {mic=1, mas=2, mev=3, mcn=4, tic=5, gmic=6};

// static std::map<std::string, mes_value> measure_map;


//' 

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
    
  // switch(measure){
  // case 1: res = mine_mic(minescore);
  //   break;
  // case 2: res = mine_mas(minescore); 
  //   break;
  // case 3: res = mine_mev(minescore);
  //   break;
  // case 4: 
  //   if (!std::isnan(eps))
  //     res = mine_mcn(minescore, eps);
  //   else
  //     res = mine_mcn_general(minescore);
  //   break;
  // case 5: res = mine_tic(minescore, norm);
  //   break;
  // case 6: res = mine_gmic(minescore, p);
  //   break;
  // default: res = NA_REAL;
  // };
  res = mine_mic(minescore);
  
  /* Free score struct */
  mine_free_score(&minescore);

  return(res);
}

