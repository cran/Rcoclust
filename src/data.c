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

double fun_transfrm(double tfij, int n, int dj, int transfrm) {
  if (transfrm==1) return((int)(tfij>0));
  if (transfrm==2 || transfrm==3) return(tfij * (log(1+(1+n)/(1+dj))));
  if (transfrm==4 || transfrm==5) return(tfij / sqrt(dj));
  return (tfij);
}

void read_data_byrow(struct mat_byrow *A_byrow, const int *Ai_byrow, const int *Aj_byrow, const double *Ax_byrow, const int *nnzi_, const int *nnzj_, const int n, const int d, const int nnz, const int transfrm) {
  A_byrow->col_index = (int **)malloc(sizeof(int*)*n);
  A_byrow->values = (double **)malloc(sizeof(double*)*n);
  A_byrow->nb_col = (int *)malloc(sizeof(int)*n);	
  for (int i=0; i<n; i++) {
    A_byrow->col_index[i] = (int *)malloc(sizeof(int)*nnzi_[i]);		
    A_byrow->values[i]    = (double *)malloc(sizeof(double)*nnzi_[i]);
    A_byrow->nb_col[i]    = nnzi_[i];
  }
  int *nuj_=NULL;
  if (transfrm==4 || transfrm==5) {
    nuj_      = (int *)malloc(sizeof(int)*(d)); // column margins
    zero_vect_int(nuj_,d);
    int c_j=0, i=0;
    for (int c=0; c<nnz; c++) {
      if (Ai_byrow[c]>i) {i=Ai_byrow[c]; c_j=0;} // next matrix row (need to be sorted by increasing values)
      int j=Aj_byrow[c];
      int xij=Ax_byrow[c];
      nuj_[j]=nuj_[j]+xij;
      c_j=c_j+1;
    }
  }
  int c_j=0;
  int i=0;
  for (int c=0; c<nnz; c++) {
    if (Ai_byrow[c]>i) {i=Ai_byrow[c]; c_j=0;} // next matrix row (need to be sorted by increasing values)
    A_byrow->col_index[i][c_j]=Aj_byrow[c];
    //A_byrow->values[i][c_j]=Ax_byrow[c];
    if (transfrm==4 || transfrm==5) {
      A_byrow->values[i][c_j] = fun_transfrm(Ax_byrow[c], n, nuj_[(int)Aj_byrow[c]], transfrm);
    } else {A_byrow->values[i][c_j] = fun_transfrm(Ax_byrow[c], n, nnzj_[(int)Aj_byrow[c]], transfrm);}
    c_j=c_j+1;
  }
  A_byrow->n = n;
  A_byrow->d = d;
  if (transfrm==3 || transfrm==5) { // normalization if
    double *xi2 = (double *)malloc(sizeof(double)*(n));
    for (int i=0; i<n; i++) {xi2[i]=0;}
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        double xij = A_byrow->values[i][index_j];
        xi2[i]=xi2[i]+xij*xij;
      }
    }
    for (int i=0; i<n; i++) { xi2[i] = sqrt(xi2[i])+1E-8; }
    for(int i=0; i<n; i++) {
      for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
        double xij = A_byrow->values[i][index_j];
        A_byrow->values[i][index_j] = xij/xi2[i];
      }
    }
  }
}

void read_data_bycol(struct mat_bycol *A_bycol, const int *Ai_bycol, const int *Aj_bycol, const double *Ax_bycol, const int *nnzi_, const int *nnzj_, const int n, const int d, const int nnz, const int transfrm) {
  A_bycol->row_index = (int **)malloc(sizeof(int*)*d);
  A_bycol->values = (double **)malloc(sizeof(double*)*d);
  A_bycol->nb_row = (int *)malloc(sizeof(int)*d);
  for (int j=0; j<d; j++) {
    A_bycol->row_index[j] = (int *)malloc(sizeof(int)*nnzj_[j]);
    A_bycol->values[j]    = (double *)malloc(sizeof(double)*nnzj_[j]);
    A_bycol->nb_row[j]    = nnzj_[j];
  }
  int* nuj_=NULL;
  if (transfrm==4 || transfrm==5) {
    nuj_      = (int *)malloc(sizeof(int)*(d)); // column margins
    zero_vect_int(nuj_,d);
    int c_i=0;
    int j=0;
    for (int c=0; c<nnz; c++) {
      if (Aj_bycol[c]>j) {j=Aj_bycol[c]; c_i=0;} // next matrix col (need to be sorted by increasing values)
      int xij=Ax_bycol[c];
      nuj_[j]=nuj_[j]+xij;
      c_i=c_i+1;
    }
  }
  int c_i=0;
  int j=0;
  for (int c=0; c<nnz; c++) {
    if (Aj_bycol[c]>j) {j=Aj_bycol[c]; c_i=0;} // next matrix col (need to be sorted by increasing values)
    A_bycol->row_index[j][c_i]=Ai_bycol[c];
    //A_bycol->values[j][c_i]=Ax_bycol[c];
    if (transfrm==4 || transfrm==5) {
      A_bycol->values[j][c_i] = fun_transfrm(Ax_bycol[c], n, nuj_[j], transfrm);
    } else {A_bycol->values[j][c_i] = fun_transfrm(Ax_bycol[c], n, nnzj_[j], transfrm);}
    c_i=c_i+1;
  }
  A_bycol->n = n;
  A_bycol->d = d;
  if (transfrm==3 || transfrm==5) { // normalization if
    double *xi2 = (double *)malloc(sizeof(double)*(n));
    for (int i=0; i<n; i++) {xi2[i]=0;}
    for(int j=0; j<d; j++) {
      for(int index_i=0; index_i<A_bycol->nb_row[j]; index_i++) {
        double xij = A_bycol->values[j][index_i];
        int i      = A_bycol->row_index[j][index_i];
        xi2[i]=xi2[i]+xij*xij;
      }
    }
    for (int i=0; i<n; i++) { xi2[i] = sqrt(xi2[i])+1E-8; }
    for(int j=0; j<d; j++) {
      for(int index_i=0; index_i<A_bycol->nb_row[j]; index_i++) {
        double xij = A_bycol->values[j][index_i];
        int i      = A_bycol->row_index[j][index_i];
        A_bycol->values[j][index_i]=xij/xi2[i];
      }
    }
  }
}

void free_data_byrow(struct mat_byrow *A_byrow, const int n) {
  for (int i=0; i<n; i++) {
    free(A_byrow->col_index[i]);
    free(A_byrow->values[i]);
  }
  free(A_byrow->col_index);
  free(A_byrow->values);
  free(A_byrow->nb_col);
  free(A_byrow);
}

void free_data_bycol(struct mat_bycol *A_bycol, const int d) {
  for (int j=0; j<d; j++) {
    free(A_bycol->row_index[j]);
    free(A_bycol->values[j]);
  }
  free(A_bycol->row_index);
  free(A_bycol->values);
  free(A_bycol->nb_row);
  free(A_bycol);
}
