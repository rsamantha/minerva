/*  
    This code is written by Davide Albanese <davide.albanese@gmail.com>.
    (C) 2012 Davide Albanese, (C) 2012 Fondazione Bruno Kessler.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "core.h"
#include "mine.h"


#define MAX(a, b) ((a) > (b) ? (a):(b))
#define MIN(a, b) ((a) < (b) ? (a):(b))

/* Computes the maximum normalized mutual information scores 
 * and returns a mine_score structure.
 */
mine_score *mine_compute_score(mine_problem *prob, mine_parameter *param)
{
  int i, j, q, p;
  double B;
  mine_score *score;
  int gx, gy, Gy_max, k;
  int *Gy;
  int *idx_x, *idx_y;
  double *x_x, *y_y, *x_y, *y_x;
  int *Qm, *Qm_s, *Pm_s;
  double *It;
  
  
  score = (mine_score *) malloc (sizeof(mine_score));
  score->m = 0;
  /* x-1 for each y */
  score->p = NULL;
  score->I = NULL;
  Gy = NULL;

  B = MAX(pow(prob->n, param->alpha), 4);
  Gy_max = MAX((int)floor(B/2.0), 2);

  score->m = Gy_max-1;
  score->p = (int *) malloc (score->m * sizeof(int));
  Gy = (int *) malloc (score->m * sizeof(int));
  for (gy=2; gy<Gy_max+1; gy++)
    {
      gx = (int) floor(B / gy);
      if ((gx * gy) <= B)
	{
	  score->p[gy-2] = gx-1;
	  Gy[gy-2] = gy;
	}
    }
    
  score->I = (double **) malloc (score->m * sizeof(double *));
  
  /* sort and argsort of x and y */ 
  x_x = (double *) malloc (prob->n * sizeof(double));
  y_y = (double *) malloc (prob->n * sizeof(double));
  idx_x = (int *) malloc (prob->n * sizeof(int));
  idx_y = (int *) malloc (prob->n * sizeof(int));
  for (i=0; i<prob->n; i++)
    {
      x_x[i] = prob->x[i];
      y_y[i] = prob->y[i];
      idx_x[i] = i;
      idx_y[i] = i;
    }
  sort(x_x, idx_x, prob->n);
  sort(y_y, idx_y, prob->n);
  
  /* sort x by y-value and sort y by x-value */ 
  x_y = (double *) malloc (prob->n * sizeof(double));
  y_x = (double *) malloc (prob->n * sizeof(double));
  for (j=0; j<prob->n; j++)
    {
      x_y[j] = prob->x[idx_y[j]];
      y_x[j] = prob->y[idx_x[j]];
    }
  
  Qm = (int *) malloc (prob->n * sizeof(int));
  Qm_s = (int *) malloc (prob->n * sizeof(int));
  Pm_s = (int *) malloc (prob->n * sizeof(int));
  
  /* x vs. y */
  for (i=0; i<score->m; i++)
    {
      score->I[i] = (double *) malloc ((score->p[i]) * sizeof(double));
      k = MAX((int) (param->c * (score->p[i]+1)), 1);
      
      q = EquipartitionYAxis(y_y, prob->n, Gy[i], Qm_s);

      for (j=0; j<prob->n; j++)
	Qm[idx_y[j]] = Qm_s[j];
      
      for (j=0; j<prob->n; j++)
	Qm_s[j] = Qm[idx_x[j]];
      
      p = GetSuperclumpsPartition(x_x, prob->n, Qm_s, k, Pm_s);
      
      ApproxOptimizeXAxis(x_x, y_x, prob->n, Qm_s, q, Pm_s, p,
			  score->p[i]+1, score->I[i]);
    }
  
  /* y vs. x */
  for (i=0; i<score->m; i++)
    {
      It = (double *) malloc ((score->p[i]) * sizeof(double));
      k = MAX((int) (param->c * (score->p[i]+1)), 1);
      
      q = EquipartitionYAxis(x_x, prob->n, Gy[i], Qm_s);
      
      for (j=0; j<prob->n; j++)
	Qm[idx_x[j]] = Qm_s[j];
      
      for (j=0; j<prob->n; j++)
	Qm_s[j] = Qm[idx_y[j]];
      
      p = GetSuperclumpsPartition(y_y, prob->n, Qm_s, k, Pm_s);

      ApproxOptimizeXAxis(y_y, x_y, prob->n, Qm_s, q, Pm_s, p,
			  score->p[i]+1, It);

      for (j=0; j<score->p[i]; j++)
	score->I[j][i] = MAX(It[j], score->I[j][i]);
      
      free(It);
    }
  
  free(Qm);
  free(Qm_s);
  free(Pm_s);
  free(x_x);
  free(y_y);
  free(x_y);
  free(y_x);
  free(idx_x);
  free(idx_y);
  free(Gy);

  return score;
}

/* This function checks the parameters. It should be called 
 * before calling mine_compute_score(). It returns NULL if 
 * the parameters are feasible, otherwise an error message is returned.
 */
char *check_parameter(mine_parameter *param)
{
  if ((param->alpha <= 0.0) || (param->alpha > 1.0))
    return "alpha must be in (0, 1.0]";

  if (param->c <= 0.0)
    return "c must be > 0.0";
  
  return NULL;
}

/* Returns the Maximal Information Coefficient (MIC). */
double mic(mine_score *score)
{
  int i, j;
  double score_max;
  
  score_max = 0.0;
  for (i=0; i<score->m; i++)
    for (j=0; j<score->p[i]; j++)
      if (score->I[i][j] > score_max)
	score_max = score->I[i][j];
    
  return score_max;
}

/* Returns the Maximum Asymmetry Score (MAS). */
double mas(mine_score *score)
{
  int i, j;
  double score_max, s;
  
  score_max = 0.0;
  for (i=0; i<score->m; i++)
    for (j=0; j<score->p[i]; j++)
      {
	s = fabs(score->I[i][j] - score->I[j][i]);
	if (s > score_max)
	    score_max = s;
      }
  
  return score_max;
}

/* Returns the Maximum Edge Value (MEV). */
double mev(mine_score *score)
{
  int i, j;
  double score_max;
  
  score_max = 0.0;
  for (i=0; i<score->m; i++)
    for (j=0; j<score->p[i]; j++)
      if (((j==0) || (i==0)) && score->I[i][j] > score_max)
	  score_max = score->I[i][j];
	   
   return score_max;
}

/* Returns the Minimum Cell Number (MCN). */
double mcn(mine_score *score)
{
  int i, j, b, b_max;
  double score_max;
  
  b_max = 4;
  score_max = -1.0;

  for (i=0; i<score->m; i++)
    for (j=0; j<score->p[i]; j++)
      {
	b = (i+2) * (j+2);
	if ((score->I[i][j] > score_max) || 
	    ((score->I[i][j] == score_max) && (b < b_max)))
	  {
	    score_max = score->I[i][j];
	    b_max = b;
	  }
      }
  return log(b_max) / log(2.0);
}

/* This function frees the memory used by a mine_score and 
 *  destroys the score structure.
 */
void mine_free_score(mine_score **score)
{
  int i;
  mine_score *score_ptr = *score;

  if (score_ptr != NULL)
    {
      if (score_ptr->m != 0)
	{
	  free(score_ptr->p);
	  for (i=0; i<score_ptr->m; i++)
	    free(score_ptr->I[i]);
	  free(score_ptr->I);
	}
      
      free(score_ptr);
    }
}
