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
SEXP ddkm(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_Ai_bycol, SEXP R_Aj_bycol, SEXP R_Ax_bycol, SEXP R_nnzj_, SEXP R_zi_, SEXP R_wj_, SEXP params_)
{
  double *Ax_byrow, *Ax_bycol;
  int *Ai_byrow, *Aj_byrow, *Ai_bycol, *Aj_bycol, *nnzi_, *nnzj_, *zi, *wj;
  PROTECT( R_Ai_byrow = coerceVector(R_Ai_byrow, INTSXP) );  Ai_byrow = INTEGER(R_Ai_byrow);
  PROTECT( R_Aj_byrow = coerceVector(R_Aj_byrow, INTSXP) );  Aj_byrow = INTEGER(R_Aj_byrow);
  PROTECT( R_Ax_byrow = coerceVector(R_Ax_byrow, REALSXP) ); Ax_byrow = REAL(R_Ax_byrow);
  PROTECT( R_Ai_bycol = coerceVector(R_Ai_bycol, INTSXP) );  Ai_bycol = INTEGER(R_Ai_bycol);
  PROTECT( R_Aj_bycol = coerceVector(R_Aj_bycol, INTSXP) );  Aj_bycol = INTEGER(R_Aj_bycol);
  PROTECT( R_Ax_bycol = coerceVector(R_Ax_bycol, REALSXP) ); Ax_bycol = REAL(R_Ax_bycol);
  PROTECT( R_nnzi_    = coerceVector(R_nnzi_, INTSXP) );     nnzi_    = INTEGER(R_nnzi_);
  PROTECT( R_nnzj_    = coerceVector(R_nnzj_, INTSXP) );     nnzj_    = INTEGER(R_nnzj_);
  PROTECT( R_zi_      = coerceVector(R_zi_, INTSXP) );       zi       = INTEGER(R_zi_);
  PROTECT( R_wj_      = coerceVector(R_wj_, INTSXP) );       wj       = INTEGER(R_wj_);
  PROTECT( params_    = coerceVector(params_, REALSXP) );
  int n               = (int)REAL(params_)[0];
  int d               = (int)REAL(params_)[1];
  int g               = (int)REAL(params_)[2];
  int nnz             = (int)REAL(params_)[3];
  double delta        = REAL(params_)[4];
  int transfrm        = (int)REAL(params_)[5];
  int maxiter         = (int)REAL(params_)[6];
  int debug           = (int)REAL(params_)[7];
  double tol          = REAL(params_)[8];
  struct mat_byrow *A_byrow = (struct mat_byrow *)malloc(sizeof(struct mat_byrow));
  struct mat_bycol *A_bycol = (struct mat_bycol *)malloc(sizeof(struct mat_bycol));
  read_data_byrow(A_byrow, Ai_byrow, Aj_byrow, Ax_byrow, nnzi_, nnzj_, n, d, nnz, transfrm);
  read_data_bycol(A_bycol, Ai_bycol, Aj_bycol, Ax_bycol, nnzi_, nnzj_, n, d, nnz, transfrm);
  if (delta<0) {
    double maxij=A_bycol->values[0][0];
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        double xij = A_byrow->values[i][index_j];
        int j      = A_byrow->col_index[i][index_j];
        if (xij>maxij) { maxij=xij;}
      }
    }
    delta=maxij;
  }
  int *wk = (int *)malloc(sizeof(int)*(g));
  int *zk = (int *)malloc(sizeof(int)*(g));
  double *e  = (double *)malloc(sizeof(double)*(maxiter));
  for(int iter=0;iter<maxiter;iter++) {e[iter]=-1;} //t[iter]=-1;}
  if (debug) Rprintf("DDKM\n");
  if (debug) Rprintf("n=%d, d=%d g=%d nnz=%d\n",n,d,g,nnz);
  if (debug) Rprintf("delta=%f, transfrm=%d, maxiter=%d\n",delta,transfrm,maxiter);
  if (debug) Rprintf("debug=%d, tol=%f\n",debug,tol);
  if (debug) Rprintf("iter=");
  int iter_end=maxiter-1;
  for(int iter=0;iter<maxiter;iter++) {
    if (debug) Rprintf("%d#",iter);
    zero_vect_int(zk,g);
    zero_vect_int(wk,g);
    compute_sums_DDKM(zi,wj,wk,zk,n,d,g);
    zero_vect_int(zi,n);
    for(int i=0; i<n; i++) {
      double valmin=DBL_MAX;
      int kmin=g;
      for(int k=0;k<g;k++) {
        double Dik=0;
        for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
          double xij = A_byrow->values[i][index_j];
          int j      = A_byrow->col_index[i][index_j];
          Dik=Dik+(int)(wj[j]==k)*xij*(xij-2*delta)/(wk[k]+tol);
        }
        Dik=Dik+delta*delta*wk[k]/(wk[k]+tol);
        if (Dik<valmin) {
          valmin=Dik;
          kmin=k;
        }
      }
      zi[i]=kmin;
    }
    zero_vect_int(zk,g);
    zero_vect_int(wk,g);
    compute_sums_DDKM(zi,wj,wk,zk,n,d,g);
    zero_vect_int(wj,d);
    for(int j=0; j<d; j++) {
      double valmin=DBL_MAX;
      int kmin=0;
      for(int k=0;k<g;k++) {
        double Djk=0;
        for(int index_i=0; index_i<A_bycol->nb_row[j]; index_i++) {
          double xij = A_bycol->values[j][index_i];
          int i      = A_bycol->row_index[j][index_i];
          Djk=Djk+(int)(zi[i]==k)*xij*(xij-2*delta)/(zk[k]+tol);
        }
        Djk=Djk+delta*delta*zk[k]/(zk[k]+tol);
        if (Djk<valmin) {
          valmin=Djk;
          kmin=k;
        }
      }
      wj[j]=kmin;
    }
    compute_sums_DDKM(zi,wj,wk,zk,n,d,g);
    e[iter]=0;
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        for(int k=0; k<g; k++) {
          double xij = A_byrow->values[i][index_j];
          int j      = A_byrow->col_index[i][index_j];
          e[iter]    = e[iter]+ ((int)(wj[j]==k)/(wk[k]+tol))*((int)(zi[i]==k)/(zk[k]+tol))*xij*(xij-2*delta);
        }
      }
    }
    e[iter]=e[iter]+delta*delta*g;
    //Rprintf("\niter=%d L=%f|",iter,e[iter]);
    if (iter>7) if (abs(e[iter]-e[iter-1])/e[iter]<1E-8) { iter_end=iter; break; }
    R_CheckUserInterrupt();
  }
  if (debug) Rprintf("\nDDKM ending...\n");
  free_data_byrow(A_byrow, n);
  free_data_bycol(A_bycol, d);
  free(zk);
  free(wk);
  //
  SEXP Robj; PROTECT(Robj = allocVector(REALSXP, iter_end+1));
  for (int iter=0; iter<=iter_end; iter++) { REAL(Robj)[iter] = e[iter]; }
  free(e);
  UNPROTECT(12);
  return Robj;
}
