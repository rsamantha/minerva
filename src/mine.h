#include <stdlib.h>


/* The mine_problem structure describes the problem. */
/* x and y are the two variables of length n. */
typedef struct mine_problem
{
  int n;
  double *x;
  double *y;
} mine_problem;


/* The mine_parameter structure describes the MINE parameters. */
/* alpha is the exponent in B(n) = n^alpha and must be in  */
/* (0,1], and c determines how many more clumps there will  */
/* be than columns in every partition. c = 15 meaning that  */
/* when trying to draw Gx grid lines on the x-axis, the  */
/* algorithm will start with at most 15*Gx clumps. c must  */
/* be > 0. */
typedef struct mine_parameter
{
  double alpha;
  double c;
} mine_parameter;


/* The mine_score structure describes the maximum  */
/* normalized mutual information scores. I[i][j]  */
/* contains the score using a grid partitioning  */
/* x-values into i+2 bins and y-values into j+2 bins.  */
/* p and I are of length m and each I[i] is of length p[i]. */
typedef struct mine_score
{
  int m;
  int *p;
  double **I;
} mine_score;


/* Computes the maximum normalized mutual information scores 
 * and returns a mine_score structure.
 */
mine_score *mine_compute_score(mine_problem *prob, 
			       mine_parameter *param);

/* This function checks the parameters. It should be called 
 * before calling mine_compute_score(). It returns NULL if 
 * the parameters are feasible, otherwise an error message is returned.
 */
char *check_parameter(mine_parameter *param);


/* Returns the Maximal Information Coefficient (MIC). */
double mic(mine_score *score);


/* Returns the Maximum Asymmetry Score (MAS). */
double mas(mine_score *score);


/* Returns the Maximum Edge Value (MEV). */
double mev(mine_score *score);


/* Returns the Minimum Cell Number (MCN). */
double mcn(mine_score *score);


/* This function frees the memory used by a mine_score and 
 *  destroys the score structure.
 */
void mine_free_score(mine_score **score);
