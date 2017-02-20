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
SEXP randomp(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_nnzj_, SEXP R_Ar, SEXP params_)
{
  double *Ax_byrow, *Ax_bycol, *zi, *wj, *Ar;
  int *Ai_byrow, *Aj_byrow, *Ai_bycol, *Aj_bycol, *nnzi_, *nnzj_;
  PROTECT( R_Ai_byrow = coerceVector(R_Ai_byrow, INTSXP) );  Ai_byrow = INTEGER(R_Ai_byrow);
  PROTECT( R_Aj_byrow = coerceVector(R_Aj_byrow, INTSXP) );  Aj_byrow = INTEGER(R_Aj_byrow);
  PROTECT( R_Ax_byrow = coerceVector(R_Ax_byrow, REALSXP) ); Ax_byrow = REAL(R_Ax_byrow);
  PROTECT( R_nnzi_    = coerceVector(R_nnzi_, INTSXP) );     nnzi_    = INTEGER(R_nnzi_);
  PROTECT( R_nnzj_    = coerceVector(R_nnzj_, INTSXP) );     nnzj_    = INTEGER(R_nnzj_);
  PROTECT( R_Ar       = coerceVector(R_Ar, REALSXP) );       Ar       = REAL(R_Ar);
  PROTECT( params_    = coerceVector(params_, REALSXP) );
  int n               = (int)REAL(params_)[0];
  int d               = (int)REAL(params_)[1];
  int nnz             = (int)REAL(params_)[2];
  int r               = (int)REAL(params_)[3];
  double sgr          = REAL(params_)[4];
  int typer           = (int)REAL(params_)[5];
  int transfrm        = (int)REAL(params_)[6];
  int debug           = (int)REAL(params_)[7];
  if (debug) Rprintf("Random projection\n");
  if (debug) Rprintf("n=%d, d=%d nnz=%d\n",n,d,nnz);
  if (debug) Rprintf("dimr=%d, sgr=%f typer=%d\n",r,sgr,typer);
  if (debug) Rprintf("transfrm=%d, debug=%d\n",transfrm,debug);
  struct mat_byrow *A_byrow = (struct mat_byrow *)malloc(sizeof(struct mat_byrow));
  read_data_byrow(A_byrow, Ai_byrow, Aj_byrow, Ax_byrow, nnzi_, nnzj_, n, d, nnz, transfrm);
  double *R_s = (double *)malloc(sizeof(double)*(d));
  zero_vect(Ar,n*r);
  for (int s=0; s<r; s++) {
    vect_random(typer, sgr, R_s, d);
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        int j      = A_byrow->col_index[i][index_j];
        double xij = A_byrow->values[i][index_j];
        Ar[i*r+s]=Ar[i*r+s]+xij*R_s[j];
      }
    }
  }
  free(R_s);
  UNPROTECT(7);  
  return R_NilValue;
}
