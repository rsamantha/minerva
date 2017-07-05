/* This code is written by Michele Filosi  <michele.filosi@gmail.com>
   Roberto Visintainer <r.visintainer@gmail.com>.
   2012 

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
#include "mine.h"
#include <R.h>
#include <Rinternals.h>


double pearson(mine_problem *myprobl){
  double r=0.0, xmean=0.0, ymean=0.0;
  double sx=0.0, sy=0.0;
  int i;
  
  for (i=0; i<myprobl->n; i++){
    xmean += myprobl->x[i];
    ymean += myprobl->y[i];
  }
  xmean /= myprobl->n;
  ymean /= myprobl->n;
  
  for (i=0; i<myprobl->n; i++){
    sx += (myprobl->x[i] - xmean) * (myprobl->x[i] - xmean);
    sy += (myprobl->y[i] - ymean) * (myprobl->y[i] - ymean);
  }
  sx = sqrt((sx/myprobl->n));
  sy = sqrt((sy/myprobl->n));
  
  for (i=0; i<myprobl->n; i++)
    r += (((myprobl->x[i] - xmean)/sx) * ((myprobl->y[i] - ymean)/sy));
  
  r /= myprobl->n;
  return(r * r);
}

SEXP mineRonevar (SEXP x, SEXP y, SEXP alpha, SEXP C, SEXP eps, SEXP est){
  
  double *restmp;
  double EPS;
  mine_problem *prob;
  mine_parameter *param;
  mine_score *minescore;
  SEXP res;
  
  PROTECT(alpha = coerceVector(alpha,REALSXP));
  PROTECT(C = coerceVector(C,INTSXP));
  PROTECT(res=allocVector(REALSXP,7));
  restmp=REAL(res);
    
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=asReal(alpha);
  param->c=asReal(C);
  param->est=asInteger(est);
  
  prob = (mine_problem *) Calloc(1,mine_problem);
  prob->n=length(x);
  prob->x=REAL(x);
  prob->y=REAL(y);
  
  minescore=mine_compute_score(prob,param);
  restmp[0]=mine_mic(minescore);
  restmp[1]=mine_mas(minescore);
  restmp[2]=mine_mev(minescore);
  if (!isNull(eps)){
    EPS=asReal(eps);
    restmp[3]=mine_mcn(minescore,EPS);  
  }else{
    restmp[3]=mine_mcn_general(minescore);
  }
  
  restmp[4]=restmp[0] - pearson(prob);

  restmp[5]=mine_gmic(minescore, -1); 
  restmp[6]=mine_tic(minescore, FALSE);
  
  /* Free */
  Free(prob);
  Free(param);
  mine_free_score(&minescore);
  UNPROTECT(3);
  return(res);
}

SEXP mineRall (SEXP x, SEXP nrx, SEXP ncx, SEXP alpha, SEXP C, SEXP eps, SEXP est)
{
  R_len_t i, j, rx, cx;
  double score, EPS;
  double **pointers;
  mine_problem *prob;
  mine_parameter *param;
  mine_score *minescore;
  SEXP res, mydim, resmic, resmas, resmev, resmcn, resmicmr, resgmic, restic, names;
  
  param = (mine_parameter *) Calloc(1,mine_parameter);
  param->alpha=asReal(alpha);
  param->c=asReal(C);
  param->est=asInteger(est);
  
  /* Matrix dimension */
  rx=asInteger(nrx);
  cx=asInteger(ncx);
  
  PROTECT(x=coerceVector(x,REALSXP));
  
  /* Initialize data matrix */
  pointers = (double **) Calloc(cx, double *);
  for (i = 0; i<cx;i++){
    pointers[i] = (double *) Calloc(rx, double);
    pointers[i] = &REAL(x)[i * rx];
  }
  
  /* Initialize result matrix */
  PROTECT(resmic=allocVector(REALSXP,cx*cx));
  PROTECT(resmas=allocVector(REALSXP,cx*cx));
  PROTECT(resmev=allocVector(REALSXP,cx*cx));
  PROTECT(resmcn=allocVector(REALSXP,cx*cx));
  PROTECT(resmicmr=allocVector(REALSXP,cx*cx));
  PROTECT(resgmic=allocVector(REALSXP,cx*cx));
  PROTECT(restic=allocVector(REALSXP,cx*cx));
  PROTECT(res=allocVector(VECSXP,7));
  
  /* Allocating result list */
  SET_VECTOR_ELT(res, 0, resmic);
  SET_VECTOR_ELT(res, 1, resmas);
  SET_VECTOR_ELT(res, 2, resmev);
  SET_VECTOR_ELT(res, 3, resmcn);
  SET_VECTOR_ELT(res, 4, resmicmr);
  SET_VECTOR_ELT(res, 5, resgmic);
  SET_VECTOR_ELT(res, 6, restic);
  
  /* Set the mine_problem */
  prob = (mine_problem *) Calloc(1,mine_problem);
  prob->n=rx;
  
  for (i = 0; i<cx; i++){
    prob->x=pointers[i];
    for (j = i; j<cx; j++){
      prob->y=pointers[j];
      /* Computing MINE scores */
      minescore=mine_compute_score(prob,param);
      score=mine_mic(minescore);
      REAL(resmic)[(cx*j) + i] = score;
      REAL(resmic)[(cx*i) + j] = score;
						
      score-=pearson(prob);
      REAL(resmicmr)[(cx*j) + i] = score;
      REAL(resmicmr)[(cx*i) + j] = score;
      
      score=mine_mas(minescore);
      REAL(resmas)[(cx*j) + i] = score;
      REAL(resmas)[(cx*i) + j] = score;

      score=mine_mev(minescore);
      REAL(resmev)[(cx*j) + i] = score;
      REAL(resmev)[(cx*i) + j] = score;
	  
      if (!isNull(eps)){
        EPS=asReal(eps);
        score=mine_mcn(minescore,EPS);
      }else{
        score=mine_mcn_general(minescore);
      }
      REAL(resmcn)[(cx*j) + i] = score;
      REAL(resmcn)[(cx*i) + j] = score;
	  
	  score=mine_gmic(minescore, -1);
      REAL(resgmic)[(cx*j) + i] = score;
      REAL(resgmic)[(cx*i) + j] = score;

      score=mine_tic(minescore, FALSE);
      REAL(restic)[(cx*j) + i] = score;
      REAL(restic)[(cx*i) + j] = score;

	  
      /* Free score */
      mine_free_score(&minescore);
    }
  }
    
  /* Set matrix dimension to be passed to R object */
  PROTECT(mydim=allocVector(INTSXP, 2));
  INTEGER(mydim)[0] = cx;
  INTEGER(mydim)[1] = cx;
  
  setAttrib(resmic, R_DimSymbol, mydim);
  setAttrib(resmas, R_DimSymbol, mydim);
  setAttrib(resmev, R_DimSymbol, mydim);
  setAttrib(resmcn, R_DimSymbol, mydim);
  setAttrib(resmicmr, R_DimSymbol, mydim);
  setAttrib(resgmic, R_DimSymbol, mydim);
  setAttrib(restic, R_DimSymbol, mydim);
  
  PROTECT(names=allocVector(STRSXP,7));
  SET_STRING_ELT(names,0,mkChar("MIC"));
  SET_STRING_ELT(names,1,mkChar("MAS"));  
  SET_STRING_ELT(names,2,mkChar("MEV"));  
  SET_STRING_ELT(names,3,mkChar("MCN"));
  SET_STRING_ELT(names,4,mkChar("MICR2"));
  SET_STRING_ELT(names,5,mkChar("GMIC"));
  SET_STRING_ELT(names,6,mkChar("TIC"));
  setAttrib(res, R_NamesSymbol,names);
  
  /* Free memeory */
  UNPROTECT(11);
  Free(pointers);
  Free(param);
  Free(prob);
  /* Return a named list of 4 matrix */
  return(res);
}

