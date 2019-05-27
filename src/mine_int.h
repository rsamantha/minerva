#ifndef MINE_INT_H
#define MINE_INT_H

const std::map<std::string, int> create_measure_map();
const std::map<std::string, int> create_est_map();


/* Create a constant map for available measures 
 * of the mine statistic.
 * Values are coded as follows:
 * {"mic", 1}, {"mas", 2}, {"mev", 3}, {"mcn", 4}, {"tic", 5}, {"gmic", 6};
 * It works with upper case also.
 */
const std::map<std::string, int> MEASURE=create_measure_map();
//



const std::map<std::string, int> EST=create_est_map();


#include <Rcpp.h>

double mine_stat(Rcpp::NumericVector x, Rcpp::NumericVector y, double alpha, double C, Rcpp::String est, Rcpp::String measure, double eps, double p, bool norm);

#endif

