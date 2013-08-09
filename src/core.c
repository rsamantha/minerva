/*
    This code is written by Davide Albanese <davide.albanese@gmail.com>.
    (C) 2012 Davide Albanese.

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
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "core.h"

#define MAX(a, b) ((a) > (b) ? (a):(b))
#define MIN(a, b) ((a) < (b) ? (a):(b))


/* computes H(Q) */
double HQ(int **N, int q, int p, int n)
{
  int i;
  double prob, logprob, H;

  H = 0.0;

  for (i=q;i--;)//(i=0;i<q; i++)//0 to q
    {
      prob = (double) N[i][p-1] / (double) n;
      if (prob != 0)
	{
	  logprob = log(prob);
	  H += prob * logprob;
	}
    }

  return -H;
}


/* computes H(<c_0, c_s, c_t>) */
double HP3(const int *np, int p, int s, int t)
{
  int sum1, sum2, tot;
  double prob1, prob2, H;

  if (s == t)
    return 0.0;

  sum1 = np[s-1];
  sum2 = np[t-1] - np[s-1];
  
  tot = np[t-1];
  prob1 = (double) sum1 / (double) tot;
  prob2 = (double) sum2 / (double) tot;

  H = 0.0;

  if (prob1 != 0)
    H += prob1 * log(prob1);

  if (prob2 != 0)
    H += prob2 * log(prob2);

  return -H;
}


/* computes H(<c_0, c_s, c_t>, Q) */
double HPQ3(int **N, const int *np, int q, int p, int s, int t)
{
  int i;
  int sum1 ,sum2, tot;
  double prob1, prob2, H;

  tot = np[t-1];

  H = 0.0;
  for (i=q;i--;)//(i = 0; i<q; i++ )
    {
      sum1 = N[i][s-1];
      sum2 = N[i][t-1] - N[i][s-1];
      prob1 = (double) sum1 / (double) tot;
      prob2 = (double) sum2 / (double) tot;

      if (prob1 != 0)
	H += prob1 * log(prob1);

      if (prob2 != 0)
	H += prob2 * log(prob2);
    }
  return -H;
}



/* computes H(<c_s, c_t>, Q) */
double HPQ2(int **N, const int *np, int q, int p, int s, int t)
{
  int i;
  int sum, tot;
  double prob, H;

  if (s == t)
    return 0.0;

  tot = np[t-1] - np[s-1];

  H = 0.0;
  for (i=q; i--; )
    {
      sum = N[i][t-1] - N[i][s-1];
      prob = (double) sum / (double) tot;
      if (prob != 0)
	H += prob * log(prob);
    }

  return -H;
}


/* Dy must be sorted in increasing order.
 * Qm must be a preallocated vector of size n.
 * Returns q.
 */
int EquipartitionYAxis(const double *Dy, int n, int y, int *Qm)
{
  int i, j, s, h, curr;
  double rowsize;

  for (i=0; i<n; i++)
    Qm[i] = -1;

  i = 0;
  curr = 0;
  rowsize =  (double) n / (double) y;
  while (i < n)
    {
      s = 1;
      for (j=i+1; j<n; j++)
	if (Dy[i] == Dy[j])
	  s++;
	else
	  break;

      h = 0;
      for (j=n;j--;)//(j=0; j<n; j++)
	if (Qm[j] == curr)
	  h++;

      if ((fabs(h+s-rowsize) >= fabs(h-rowsize)) && (h != 0))
	{
	  curr++;
	  rowsize = (double) (n-i) / (double) (y-curr);
	}

      for (j=s; j--;)//(j=0; j<s; j++)
	Qm[i+j] = curr;

      i += s;
    }

  return curr + 1;
}


/* Dx and Qm must be sorted in increasing order by Dx-values.
 * Pm must be a preallocated vector of size n.
 * Returns p.
 */
int GetSuperclumpsPartition(const double *Dx, int n, const int *Qm, int k, int *Pm)
{
  int i, j, s, h, b, flag, curr;
  int *Qm_n, *Pm_c;
  double colsize;

  /* clumps */

  Qm_n = (int *) malloc (n * sizeof(int));
  for(i=0; i<n; i++)
    Qm_n[i] = Qm[i];
  
  i = 0;
  b = -1;
  while (i < n)
    {
      s = 1;
      flag = 0;
      for (j=i+1; j<n; j++)
	if (Dx[i] == Dx[j])
	  {
	    s++;
	    if (Qm_n[i] != Qm_n[j])
	      flag = 1;
	  }
	else
	  break;

      if ((s > 1) && (flag == 1))
	{
	  for (j=s; j--;)
	    Qm_n[i+j] = b;
	  b--;
	}

      i += s;
    }

  Pm_c = (int *) malloc (n * sizeof(int));
  curr = 0;
  Pm_c[0] = curr;
  for (i=1; i<n; i++)
    {
      if (Qm_n[i] != Qm_n[i-1])
	curr++;
      Pm_c[i] = curr;
    }

  free(Qm_n);

  if (k >= (curr + 1))
    {
      for (i=0; i<n; i++)
	Pm[i] = Pm_c[i];
      free(Pm_c);
      return curr + 1;
    }

  /* superclumps*/

  for (i=0; i<n; i++)
    Pm[i] = -1;

  i = 0;
  curr = 0;
  colsize =  (double) n / (double) k;
  while (i < n)
    {
      s = 1;
      for (j=i+1; j<n; j++)
	if (Pm_c[i] == Pm_c[j])
	  s++;
	else
	  break;

      h = 0;
      for (j=0; j<n; j++)
	if (Pm[j] == curr)
	  h++;

      if ((fabs(h+s-colsize) >= fabs(h-colsize)) && (h != 0))
	{
	  curr++;
	  colsize = (double) (n-i) / (double) (k-curr);
	}

      for (j = s; j--;)//(j=0; j<s; j++)
	Pm[i+j] = curr;

      i += s;
    }

  free (Pm_c);
  return curr + 1;
}


/*
 * Dx, Dy, Qm and Pm must be sorted in increasing order
 * by Dx-value.
 * I must be a preallocated array of dimension x-1.
 * It will contain I_{k,2}, ..., I_{k, x}.
 */
void ApproxOptimizeXAxis(const double *Dx,const double *Dy, int n,
			 const int *Qm, int q, const int *Pm, int p,
			 int x, double *I)
{
  int i, j, s, t, l;
  int **N;
  int *np; 
  int *c;
  double **IM; // mutual information matrix, I_t,l
  double f, fmax, r1, r2;
  double hq;
  double **hpq2;


  /* if Pk==1 return I=0 */
  if (p <= 1)
    {
      for (i=0; i<x-1; i++)
	I[i] = 0.0;
      return;
    }

  /* N matrix */
  N = (int **) malloc (q * sizeof(int *));
  for (i=0; i<q; i++)
    {
      N[i] = (int *) malloc (p * sizeof(int));
      for (j=0; j<p; j++)
	N[i][j] = 0;
    }

  /* np vector */
  np = (int *) malloc (p * sizeof(int));
  for (j=0; j<p; j++)
    np[j] = 0;
  
  /* /\* c_1, ..., c_k *\/ */
  /* c = (int *) malloc (p * sizeof(int)); */
  /* for (j=0; j<p; j++) */
  /*   c[j] = 0; */
  
  /* fill N and np */
  for (i=0; i<n; i++)
    {
      N[Qm[i]][Pm[i]] += 1;
      np[Pm[i]] += 1;
      //c[Pm[i]] += 1;
    }

  // Update Filosi 12/12/12
  /* Entropy function will be updated too... */
  /* Create cumulative distribution */
  for (i=0; i<q; i++)
    for (j=1; j<p; j++)
      N[i][j] += N[i][j-1];
 
  
  /* c_1, ..., c_k */
  c = (int *) malloc (p * sizeof(int));
  c[0] = np[0];
  for (i=1; i<p; i++)
    {
      //c[i] += c[i-1];
      c[i] = np[i] + c[i-1];
      np[i] += np[i-1];
    }
  
  /* IM matrix */
  IM = (double **) malloc ((p+1) * sizeof(double *));
  for (i=0; i<=p; i++)
    {
      IM[i] = (double *) malloc ((x+1) * sizeof(double));
      for (j=0; j<=x; j++)
        IM[i][j] = 0.0;
    }
  
  hq = HQ(N, q, p, n);
  /* Find the optimal partitions of size 2 */
  for (t=2; t<=p; t++)
    {
      fmax = -DBL_MAX;
      for (s=1; s<=t; s++)
	{
	  f = HP3(np, p, s, t) - HPQ3(N, np, q, p, s, t);
	  if (f > fmax)
	    {
	      IM[t][2] = hq + f;
	      fmax = f;
	    }
	}
    }
  
  /* precomputed H(<c_s, c_t>, Q) matrix */
  hpq2 = (double **) malloc ((p+1) * sizeof(double *));
  for (i=0; i<=p; i++)
    {
      hpq2[i] = (double *) malloc ((p+1) * sizeof(double));
      for (j=0; j<=p; j++)
	hpq2[i][j] = 0.0;
    }
  for (t=3; t<=p; t++)
    for (s=2; s<=t; s++)
      hpq2[s][t] = HPQ2(N, np, q, p, s, t);
  
  /* inductively build the rest of the table of optimal partitions */
  for (l=3; l<=x; l++)
    for (t=l; t<=p; t++)
      {
	fmax = -DBL_MAX;
	for (s=l-1; s<=t; s++)
	  {
	    r1 = (double) c[s-1] / (double) c[t-1];
	    r2 = (double) (c[t-1] - c[s-1]) / (double) c[t-1];
	    f = (r1 * (IM[s][l-1] - hq)) - (r2 * hpq2[s][t]);

	    if (f > fmax)
	      {
		IM[t][l] = hq + f;
		fmax = f;
	      }
	  }
      }

  /* line 19 */
  if (x > p)
    {
      for (i=p+1; i<=x; i++)
	IM[p][i] = IM[p][p];
    }

  /* fill I */
  for (i=2; i<=x; i++)
    I[i-2] = IM[p][i] / MIN(log(i), log(q));

  /* free */
  for (i=q;i--;)//(i=0; i<q; i++)
    free(N[i]);
  free(N);
  free(np);
  free(c);
  for (i=p+1;i--;)//(i=0; i<=p; i++)
    {
      free(IM[i]);
      free(hpq2[i]);
    }
  free(IM);
  free(hpq2);

  return;
}


void swap(double *x, int *idx, int i, int j)
{
  double x_t;
  int idx_t;

  x_t = x[i];
  x[i] = x[j];
  x[j] = x_t;

  idx_t = idx[i];
  idx[i] = idx[j];
  idx[j] = idx_t;
}


void quicksort(double *x, int *idx, int l, int u)
{
  int i, m;

  if (l >= u)
    return;

  m = l;
  for (i=l+1; i<=u; i++)
    if (x[i] < x[l])
      swap(x, idx, ++m, i);
  swap(x, idx, l, m);
  quicksort(x, idx, l, m-1);
  quicksort(x, idx, m+1, u);
}


void sort(double *x, int *idx, int n)
{
  quicksort(x, idx, 0, n-1);
}
