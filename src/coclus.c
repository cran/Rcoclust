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
//SEXP coclus(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_Ai_bycol, SEXP R_Aj_bycol, SEXP R_Ax_bycol, SEXP R_nnzj_, SEXP R_zi_, SEXP R_wj_, SEXP params_)
SEXP coclus(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_nnzj_, SEXP R_zi_, SEXP R_wj_, SEXP params_)
{
  double *Ax_byrow, *Ax_bycol;
  int *Ai_byrow, *Aj_byrow, *Ai_bycol, *Aj_bycol, *nnzi_, *nnzj_, *zi, *wj;
  PROTECT( R_Ai_byrow = coerceVector(R_Ai_byrow, INTSXP) );  Ai_byrow = INTEGER(R_Ai_byrow);
  PROTECT( R_Aj_byrow = coerceVector(R_Aj_byrow, INTSXP) );  Aj_byrow = INTEGER(R_Aj_byrow);
  PROTECT( R_Ax_byrow = coerceVector(R_Ax_byrow, REALSXP) ); Ax_byrow = REAL(R_Ax_byrow);
  //PROTECT( R_Ai_bycol = coerceVector(R_Ai_bycol, INTSXP) );  Ai_bycol = INTEGER(R_Ai_bycol);
  //PROTECT( R_Aj_bycol = coerceVector(R_Aj_bycol, INTSXP) );  Aj_bycol = INTEGER(R_Aj_bycol);
  //PROTECT( R_Ax_bycol = coerceVector(R_Ax_bycol, REALSXP) ); Ax_bycol = REAL(R_Ax_bycol);
  PROTECT( R_nnzi_    = coerceVector(R_nnzi_, INTSXP) );     nnzi_    = INTEGER(R_nnzi_);
  PROTECT( R_nnzj_    = coerceVector(R_nnzj_, INTSXP) );     nnzj_    = INTEGER(R_nnzj_);
  PROTECT( R_zi_      = coerceVector(R_zi_, INTSXP) );       zi       = INTEGER(R_zi_);
  PROTECT( R_wj_      = coerceVector(R_wj_, INTSXP) );       wj       = INTEGER(R_wj_);
  PROTECT( params_    = coerceVector(params_, REALSXP) );
  int n               = (int)REAL(params_)[0];
  int d               = (int)REAL(params_)[1];
  int g               = (int)REAL(params_)[2];
  int nnz             = (int)REAL(params_)[3];
  int transfrm        = (int)REAL(params_)[4];
  int maxiter         = (int)REAL(params_)[5];
  int debug           = (int)REAL(params_)[6];
  double tol          = REAL(params_)[7];
  struct mat_byrow *A_byrow = (struct mat_byrow *)malloc(sizeof(struct mat_byrow));
  //struct mat_bycol *A_bycol = (struct mat_bycol *)malloc(sizeof(struct mat_bycol));
  read_data_byrow(A_byrow, Ai_byrow, Aj_byrow, Ax_byrow, nnzi_, nnzj_, n, d, nnz, transfrm);
  //read_data_bycol(A_bycol, Ai_bycol, Aj_bycol, Ax_bycol, nnzi_, nnzj_, n, d, nnz, transfrm);
  double *xik    = (double *)malloc(sizeof(double)*(n*g));
  double *xjk    = (double *)malloc(sizeof(double)*(d*g));
  double *xio    = (double *)malloc(sizeof(double)*(n));
  double *xoj    = (double *)malloc(sizeof(double)*(d));
  double *xok    = (double *)malloc(sizeof(double)*(g));
  double *xko    = (double *)malloc(sizeof(double)*(g));
  double xoo     = 0;
  zero_vect(xio,n);
  zero_vect(xoj,d);
  for(int i=0; i<n; i++) {
    for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
      double xij = A_byrow->values[i][index_j];
      int j      = A_byrow->col_index[i][index_j];
      xio[i]=xio[i]+xij;
      xoj[j]=xoj[j]+xij;
      xoo=xoo+xij;
    }
  }
  zero_vect(xik,g);
  zero_vect(xjk,g);
  double *e  = (double *)malloc(sizeof(double)*(maxiter));
  for(int iter=0;iter<maxiter;iter++) {e[iter]=-1;} //t[iter]=-1;}
  if (debug) Rprintf("COCLUS\n");
  if (debug) Rprintf("n=%d, d=%d g=%d nnz=%d\n",n,d,g,nnz);
  if (debug) Rprintf("transfrm=%d, maxiter=%d\n",transfrm,maxiter);
  if (debug) Rprintf("debug=%d, tol=%f\n",debug,tol);
  if (debug) Rprintf("iter=");
  int iter_end=maxiter-1;
  for(int iter=0;iter<maxiter;iter++) {
    if (debug) Rprintf("%d#",iter);
    zero_vect(xik,n*g);
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        double xij = A_byrow->values[i][index_j];
        int j      = A_byrow->col_index[i][index_j];
        for (int k=0; k<g; k++) { xik[i*g+k]=xik[i*g+k]+xij*(wj[j]==k); }
      }
    }
    zero_vect(xok,g);
    for (int k=0; k<g; k++) for(int j=0; j<d; j++) { xok[k]=xok[k]+(wj[j]==k)*xoj[j]; }
    zero_vect_int(zi,n);
    for (int i=0; i<n; i++) {
      double valmax=xik[i*g+0]-xio[i]*xok[0]/xoo;
      int kmax=0;
      for (int k=1; k<g; k++) {
        if (xik[i*g+k]-xio[i]*xok[k]/xoo > valmax) {
          valmax=xik[i*g+k]-xio[i]*xok[k]/xoo;
          kmax=k;         
        }
      }
      zi[i]=kmax;
    }
    zero_vect(xjk,d*g);
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        double xij = A_byrow->values[i][index_j];
        int j      = A_byrow->col_index[i][index_j];
        for (int k=0; k<g; k++) { xjk[j*g+k]=xjk[j*g+k]+xij*(zi[i]==k); }
      }
    }
    zero_vect(xko,g);
    for (int k=0; k<g; k++) for(int i=0; i<n; i++) { xko[k]=xko[k]+(zi[i]==k)*xio[i]; }
    zero_vect_int(wj,d);
    for (int j=0; j<d; j++) {
      double valmax=xjk[j*g+0]-xoj[j]*xok[0]/xoo;
      int kmax=0;
      for (int k=1; k<g; k++) {
        if (xjk[j*g+k]-xoj[j]*xko[k]/xoo > valmax) {
          valmax=xjk[j*g+k]-xoj[j]*xko[k]/xoo;
          kmax=k;         
        }
      }
      wj[j]=kmax;
    }
    //compute_sums_FDDKM(zik,wjk,wkbeta,zkbeta,alpha,beta,n,d,g,tol);
    e[iter]=0;
    for (int j=0; j<d; j++) for (int k=0; k<g; k++) { 
      e[iter] = e[iter] + (xjk[j*g+k]-xoj[j]*xko[k]/xoo) * (wj[j]==k);
    }
    e[iter] = e[iter]/xoo;
    if (iter>7) if (abs(e[iter]-e[iter-1])/e[iter]<1E-8) { iter_end=iter; break; }
    R_CheckUserInterrupt();
  }
  if (debug) Rprintf("\nCOCLUS ending...\n");
  free_data_byrow(A_byrow, n);
  //free_data_bycol(A_bycol, d);
  free(xik);
  free(xjk);
  free(xio);
  free(xoj);
  //
  SEXP Robj; PROTECT(Robj = allocVector(REALSXP, iter_end+1));
  for (int iter=0; iter<=iter_end; iter++) { REAL(Robj)[iter] = e[iter]; }
  free(e);
  UNPROTECT(9);
  return Robj;
}
