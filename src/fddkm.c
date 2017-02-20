/*==========================================================================
 Copyright (C) R. Priam 
Date Version 0.1.0 : 2017-01-31
Date last update   : 2017-02-16
==========================================================================
This file is part of Rcoclust.

Rcoclust is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Rcoclust is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Rcoclust.  If not, see <http://www.gnu.org/licenses/>.
========================================================================== */

#include "rcoclust.h"
SEXP fddkm(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_Ai_bycol, SEXP R_Aj_bycol, SEXP R_Ax_bycol, SEXP R_nnzj_, SEXP R_zik_, SEXP R_wjk_, SEXP R_zi_, SEXP R_wj_, SEXP params_)
{
  double *Ax_byrow, *Ax_bycol, *zik, *wjk;
  int *Ai_byrow, *Aj_byrow, *Ai_bycol, *Aj_bycol, *nnzi_, *nnzj_, *zi, *wj;
  PROTECT( R_Ai_byrow = coerceVector(R_Ai_byrow, INTSXP) );  Ai_byrow = INTEGER(R_Ai_byrow);
  PROTECT( R_Aj_byrow = coerceVector(R_Aj_byrow, INTSXP) );  Aj_byrow = INTEGER(R_Aj_byrow);
  PROTECT( R_Ax_byrow = coerceVector(R_Ax_byrow, REALSXP) ); Ax_byrow = REAL(R_Ax_byrow);
  PROTECT( R_Ai_bycol = coerceVector(R_Ai_bycol, INTSXP) );  Ai_bycol = INTEGER(R_Ai_bycol);
  PROTECT( R_Aj_bycol = coerceVector(R_Aj_bycol, INTSXP) );  Aj_bycol = INTEGER(R_Aj_bycol);
  PROTECT( R_Ax_bycol = coerceVector(R_Ax_bycol, REALSXP) ); Ax_bycol = REAL(R_Ax_bycol);
  PROTECT( R_nnzi_    = coerceVector(R_nnzi_, INTSXP) );     nnzi_    = INTEGER(R_nnzi_);
  PROTECT( R_nnzj_    = coerceVector(R_nnzj_, INTSXP) );     nnzj_    = INTEGER(R_nnzj_);
  PROTECT( R_zik_     = coerceVector(R_zik_, REALSXP) );     zik      = REAL(R_zik_);
  PROTECT( R_wjk_     = coerceVector(R_wjk_, REALSXP) );     wjk      = REAL(R_wjk_);
  PROTECT( R_zi_      = coerceVector(R_zi_, INTSXP) );       zi       = INTEGER(R_zi_);
  PROTECT( R_wj_      = coerceVector(R_wj_, INTSXP) );       wj       = INTEGER(R_wj_);
  PROTECT( params_    = coerceVector(params_, REALSXP) );
  int n               = (int)REAL(params_)[0];
  int d               = (int)REAL(params_)[1];
  int g               = (int)REAL(params_)[2];
  int nnz             = (int)REAL(params_)[3];
  double alpha        = REAL(params_)[4];
  double beta         = REAL(params_)[5];
  double delta        = REAL(params_)[6];
  int transfrm        = (int)REAL(params_)[7];
  int maxiter         = (int)REAL(params_)[8];
  int debug           = (int)REAL(params_)[9];
  double tol          = REAL(params_)[10];
  struct mat_byrow *A_byrow = (struct mat_byrow *)malloc(sizeof(struct mat_byrow));
  struct mat_bycol *A_bycol = (struct mat_bycol *)malloc(sizeof(struct mat_bycol));
  read_data_byrow(A_byrow, Ai_byrow, Aj_byrow, Ax_byrow, nnzi_, nnzj_, n, d, nnz, transfrm);
  read_data_bycol(A_bycol, Ai_bycol, Aj_bycol, Ax_bycol, nnzi_, nnzj_, n, d, nnz, transfrm);
  if (delta<0) {
    double maxij=A_bycol->values[0][0];
    for(int j=0; j<d; j++) {
      for(int index_i=0; index_i<A_bycol->nb_row[j]; index_i++) {
        double xij = A_bycol->values[j][index_i];
        if (xij>maxij) maxij=xij;
      }
    }
    delta=maxij;
  }
  double *Dik    = (double *)malloc(sizeof(double)*(g));
  double *Djk    = (double *)malloc(sizeof(double)*(g));
  double *wkbeta = (double *)malloc(sizeof(double)*(g));
  double *zkbeta = (double *)malloc(sizeof(double)*(g));
  double *e  = (double *)malloc(sizeof(double)*(maxiter));
  for(int iter=0;iter<maxiter;iter++) {e[iter]=-1;} //t[iter]=-1;}
  if (debug) Rprintf("F-DDKM\n");
  if (debug) Rprintf("n=%d, d=%d g=%d nnz=%d\n",n,d,g,nnz);
  if (debug) Rprintf("alpha=%f, beta=%f, delta=%f\n",alpha, beta,delta);
  if (debug) Rprintf("transfrm=%d, maxiter=%d\n",transfrm,maxiter);
  if (debug) Rprintf("debug=%d, tol=%f\n",debug,tol);
  if (debug) Rprintf("iter=");
  int iter_end=maxiter-1;
  for(int iter=0;iter<maxiter;iter++) {
    if (debug) Rprintf("%d#",iter);
    compute_sums_FDDKM(zik,wjk,wkbeta,zkbeta,alpha,beta,n,d,g,tol);
    for(int i=0; i<n; i++) {
      for(int k=0;k<g;k++) {
        Dik[k]=0;
        for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
          double xij = A_byrow->values[i][index_j];
          int j      = A_byrow->col_index[i][index_j];
          Dik[k]=Dik[k]+pow(wjk[j*g+k],beta)*xij*(xij-2*delta)/wkbeta[k];
        }
        Dik[k]=Dik[k]+delta*delta*wkbeta[k]/wkbeta[k];
        Dik[k]=Dik[k]+tol;
      }
      for(int k=0; k<g; k++) {
        zik[i*g+k]=0;
        for (int l=0; l<g; l++) {
          zik[i*g+k] = zik[i*g+k] + pow(Dik[k]/Dik[l],1/(alpha-1));
        }
        zik[i*g+k] = 1/(zik[i*g+k]);
      }
    }
    bound_values(zik,n*g,tol);
    compute_sums_FDDKM(zik,wjk,wkbeta,zkbeta,alpha,beta,n,d,g,tol);
    for(int j=0; j<d; j++) {
      for(int k=0;k<g;k++) {
        Djk[k]=0;
        for(int index_i=0; index_i<A_bycol->nb_row[j]; index_i++) {
          double xij = A_bycol->values[j][index_i];
          int i      = A_bycol->row_index[j][index_i];
          Djk[k]=Djk[k]+pow(zik[i*g+k],beta)*xij*(xij-2*delta)/zkbeta[k];
        }
        Djk[k]=Djk[k]+delta*delta*zkbeta[k]/zkbeta[k];
        Djk[k]=Djk[k]+tol;
      }
      for(int k=0; k<g; k++) {
        wjk[j*g+k]=0;
        for (int l=0; l<g; l++) {
          wjk[j*g+k] = wjk[j*g+k] + pow(Djk[k]/Djk[l],1/(beta-1));
        }
        wjk[j*g+k] = 1/(wjk[j*g+k]);
      }
    }
    bound_values(wjk,d*g,tol);
    compute_sums_FDDKM(zik,wjk,wkbeta,zkbeta,alpha,beta,n,d,g,tol);
    e[iter]=0;
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        for(int k=0; k<g; k++) {
          double xij = A_byrow->values[i][index_j];
          int j      = A_byrow->col_index[i][index_j];
          e[iter]    = e[iter]+ (pow(wjk[j*g+k],beta)/wkbeta[k])*(pow(zik[i*g+k],alpha)/zkbeta[k])*xij*(xij-2*delta);
        }
      }
    }
    e[iter]=e[iter]+delta*delta*g;
    if (iter>7) if (abs(e[iter]-e[iter-1])/e[iter]<1E-8) { iter_end=iter; break; }
    R_CheckUserInterrupt();
  }
  if (debug) Rprintf("\nF-DDKM ending...\n");
  //
  compute_zi(zik,zi,n,g,tol);
  compute_zi(wjk,wj,d,g,tol);
  //
  free_data_byrow(A_byrow, n);
  free_data_bycol(A_bycol, d);
  free(Dik);
  free(Djk);
  free(wkbeta); 
  free(zkbeta);
  //
  SEXP Robj; PROTECT(Robj = allocVector(REALSXP, iter_end+1));
  for (int iter=0; iter<=iter_end; iter++) { REAL(Robj)[iter] = e[iter]; }
  free(e);
  UNPROTECT(14);
  return Robj;
}
