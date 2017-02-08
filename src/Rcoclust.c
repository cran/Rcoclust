/*==========================================================================
 Copyright (C) R. Priam 
 Date : January 2017
 ==========================================================================
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
 ========================================================================== */

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <Rmath.h>
#include <math.h>

#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

#define SWAP(a,b)     { double temp=(a);(a)=(b);(b)=temp; }
#define SWAP_INT(a,b) { int temp=(a);(a)=(b);(b)=temp; }

// ----------------------------------------------------------------------------

void printMatrix(double *M, int nrow, int ncol) {
  for (int i=0; i<nrow; i++) {
    for (int j=0; j<ncol; j++) { 
      Rprintf("%f ", M[i*ncol+j]);
    }
    Rprintf("\n");
  }
}

void printMatrixInt(int *M, int nrow, int ncol) {
  for (int i=0; i<nrow; i++) {
    for (int j=0; j<ncol; j++) { 
      Rprintf("%d ", M[i*ncol+j]);
    }
    Rprintf("\n");
  }
}

double fun_transfrm(double tfij, int n, int dj, int transfrm) {
  if (transfrm==1) return((int)(tfij>0));
  if (transfrm==2 || transfrm==3) return(tfij * (log(1+(1+n)/(1+dj))));
  if (transfrm==4 || transfrm==5) return(tfij / sqrt(dj));
  return (tfij);
}

// --------------------------------------------------------------------------------------------------------

void compute_zi(const double *zik, int *zi, const int n, const int g, const double tol) {
  int kimax;
  double valimax;
  for (int i=0; i<n; i++) {
    kimax=0;
    valimax=zik[i*g+0];
    for (int k=1; k<g; k++) if (zik[i*g+k]>valimax) {valimax=zik[i*g+k]; kimax=k;}
    zi[i]=kimax;
  }
}

// --------------------------------------------------------------------------------------------------------

struct mat_byrow {
  int **col_index;
  double **values;
  int *nb_col;
  int n;
  int d;
};

struct mat_bycol {
  int **row_index;
  double **values;
  int *nb_row;
  int n;
  int d;
};

void zero_vect(double *v, const int dim) { for (int l=0; l<dim; l++) v[l]=0+1E-7; }
void zero_vect_int(int *v, const int dim) { for (int l=0; l<dim; l++) v[l]=0; }
void copy_vect(double *v_out, const double *v_in, const int dim) {for (int l=0; l<dim; l++) v_out[l]=v_in[l];}

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

// --------------------------------------------------------------------------------------------------------

void bound_values(double *v, int n, double tol) {
  for (int i=0; i<n; i++) {
    if (v[i]>1-tol) v[i]=1-tol;
    if (v[i]<tol) v[i]=tol;
  }
}

void prod_vect_mat(double *v_out, double *v_in, double *vMat_in, int dim1, int dim2) {
  zero_vect(v_out,dim2);
  for (int l=0; l<dim2; l++) for (int c=0; c<dim1; c++) v_out[l]=v_out[l]+vMat_in[c*dim2+l]*v_in[c];
}

void sum_colmat2vect(double *v_out, const double *vectMat_in, const int dim1, const int dim2) {
  zero_vect(v_out,dim2);
  for (int l=0; l<dim2; l++) for (int c=0; c<dim1; c++) v_out[l]=v_out[l]+vectMat_in[c*dim2+l];
}

void normalize_posterior(double *pik, int n, int g, int withmax, int withbound, double tol) {
  for (int i=0; i<n; i++) {
    if (withmax==1) { 
      double pikmax=pik[i*g+0];
      for (int k=1; k<g; k++) pikmax=MAX(pik[i*g+k],pikmax);
      for (int k=0; k<g; k++) pik[i*g+k] = exp(pik[i*g+k] - pikmax);
    }
    double sumi=tol;
    for (int k=0; k<g; k++) {sumi=sumi+pik[i*g+k];}
    for (int k=0; k<g; k++) {pik[i*g+k] = pik[i*g+k]/sumi;}
  }
  if (withbound==1) bound_values(pik,n*g,tol);
}

void prod_A_vectB(const struct mat_byrow *A_byrow, 
                  const double *vectB, double *v_out, const int n, const int d, 
                  const int byrow, const int dim) {
  if (byrow==1) {zero_vect(v_out,d*dim);} else {zero_vect(v_out,n*dim);}
  for (int i=0; i<n; i++) {
    for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
      double xij = A_byrow->values[i][index_j];
      int j      = A_byrow->col_index[i][index_j];
      if (byrow==1) {
        for (int k=0; k<dim; k++) v_out[j*dim+k]=v_out[j*dim+k]+vectB[i*dim+k]*xij;
      } else {
        for (int l=0; l<dim; l++) v_out[i*dim+l]=v_out[i*dim+l]+vectB[j*dim+l]*xij;
      }}}
}

double innerprod_vectA_vectB(const double* vectA, const double* vectB, const int n) {
  double inner=0.0f;
  for (int i=0; i<n; i++) inner=inner+vectA[i]*vectB[i];
  return (inner);
}

double euclidist_vectA_vectB(const double* vectA, const double* vectB, const int n) {
  double euclid=0.0f;
  for (int i=0; i<n; i++) euclid=euclid+(vectA[i]-vectB[i])*(vectA[i]-vectB[i]);
  return (euclid);
}

// --------------------------------------------------------------------------------------------------------

// COCLUS
SEXP COCLUS(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_Ai_bycol, SEXP R_Aj_bycol, SEXP R_Ax_bycol, SEXP R_nnzj_, SEXP R_zi_, SEXP R_wj_, SEXP params_)
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
  int transfrm        = (int)REAL(params_)[4];
  int maxiter         = (int)REAL(params_)[5];
  int debug           = (int)REAL(params_)[6];
  double tol          = REAL(params_)[7];
  struct mat_byrow *A_byrow = (struct mat_byrow *)malloc(sizeof(struct mat_byrow));
  struct mat_bycol *A_bycol = (struct mat_bycol *)malloc(sizeof(struct mat_bycol));
  read_data_byrow(A_byrow, Ai_byrow, Aj_byrow, Ax_byrow, nnzi_, nnzj_, n, d, nnz, transfrm);
  read_data_bycol(A_bycol, Ai_bycol, Aj_bycol, Ax_bycol, nnzi_, nnzj_, n, d, nnz, transfrm);
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
    if (iter>15) if (abs(e[iter]-e[iter-1])/e[iter]<1E-6) { iter_end=iter; break; }
  }
  if (debug) Rprintf("\nCOCLUS ending...\n");
  free_data_byrow(A_byrow, n);
  free_data_bycol(A_bycol, d);
  free(xik);
  free(xjk);
  free(xio);
  free(xoj);
  //
  SEXP RlogL; PROTECT(RlogL = allocVector(REALSXP, iter_end+1));
  for (int iter=0; iter<=iter_end; iter++) { REAL(RlogL)[iter] = e[iter]; }
  free(e);
  UNPROTECT(12);
  return RlogL;
}

// --------------------------------------------------------------------------------------------------------

void compute_sums_FDDKM(const double* zik,const double* wjk, double* wkbeta, double* zkbeta, const double alpha, const double beta, const int n, const int d, const int g, const double tol) {
  for(int k=0;k<g;k++) {
    wkbeta[k]=0;
    for(int j=0;j<d;j++) wkbeta[k]=wkbeta[k]+pow(wjk[j*g+k],beta);
    wkbeta[k]=wkbeta[k]+tol;
  }
  for(int k=0;k<g;k++) {
    zkbeta[k]=0;
    for(int i=0;i<n;i++) zkbeta[k]=zkbeta[k]+pow(zik[i*g+k],alpha);
    zkbeta[k]=zkbeta[k]+tol;
  }
}

void compute_sums_DDKM(const int* zi,const int* wj, int* wk, int* zk, const int n, const int d, const int g)
{
  zero_vect_int(zk,g);
  zero_vect_int(wk,g);
  for(int k=0;k<g;k++) for(int j=0;j<d;j++) wk[k]=wk[k]+(int)(wj[j]==k);
  for(int k=0;k<g;k++) for(int i=0;i<n;i++) zk[k]=zk[k]+(int)(zi[i]==k);
}

// DDKM
SEXP DDKM(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_Ai_bycol, SEXP R_Aj_bycol, SEXP R_Ax_bycol, SEXP R_nnzj_, SEXP R_zi_, SEXP R_wj_, SEXP params_)
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
    if (iter>15) if (abs(e[iter]-e[iter-1])/e[iter]<1E-6) { iter_end=iter; break; }
  }
  if (debug) Rprintf("\nDDKM ending...\n");
  free_data_byrow(A_byrow, n);
  free_data_bycol(A_bycol, d);
  free(zk);
  free(wk);
  //
  SEXP RlogL; PROTECT(RlogL = allocVector(REALSXP, iter_end+1));
  for (int iter=0; iter<=iter_end; iter++) { REAL(RlogL)[iter] = e[iter]; }
  free(e);
  UNPROTECT(12);
  return RlogL;
}

// F-DDKM
SEXP FDDKM(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_Ai_bycol, SEXP R_Aj_bycol, SEXP R_Ax_bycol, SEXP R_nnzj_, SEXP R_zik_, SEXP R_wjk_, SEXP R_zi_, SEXP R_wj_, SEXP params_)
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
    if (iter>15) if (abs(e[iter]-e[iter-1])/e[iter]<1E-6) { iter_end=iter; break; }
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
  SEXP RlogL; PROTECT(RlogL = allocVector(REALSXP, iter_end+1));
  for (int iter=0; iter<=iter_end; iter++) { REAL(RlogL)[iter] = e[iter]; }
  free(e);
  UNPROTECT(14);
  return RlogL;
}

// --------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------

// Random projection
SEXP RPMAT(SEXP R_Ai_byrow, SEXP R_Aj_byrow, SEXP R_Ax_byrow, SEXP R_nnzi_, SEXP R_nnzj_, SEXP R_Ar, SEXP params_)
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
  int r               = (int)REAL(params_)[3]; //dim reduc random projection
  double sgr          = REAL(params_)[4];      //var reduc random projection
  int typer           = (int)REAL(params_)[5]; //type reduc random projection
  int transfrm        = (int)REAL(params_)[6];
  int debug           = (int)REAL(params_)[7];
  if (debug) Rprintf("Random projection\n");
  if (debug) Rprintf("n=%d, d=%d nnz=%d\n",n,d,nnz);
  if (debug) Rprintf("dimr=%d, sgr=%f typer=%d\n",r,sgr,typer);
  if (debug) Rprintf("transfrm=%d, debug=%d\n",transfrm,debug);
  struct mat_byrow *A_byrow = (struct mat_byrow *)malloc(sizeof(struct mat_byrow));
  read_data_byrow(A_byrow, Ai_byrow, Aj_byrow, Ax_byrow, nnzi_, nnzj_, n, d, nnz, transfrm);
  double *R      = (double *)malloc(sizeof(double)*(d*r));
  GetRNGstate ();
  for (int j=0; j<d; j++) for (int s=0; s<r; s++) {
    if (typer==0) {
      R[j*r+s]=rnorm (0.0 ,sgr);
    } else if (typer==1) {
      int u=runif(0,1);
      if (u<=0.5) {R[j*r+s]=-1;} else{R[j*r+s]=+1;} 
    } else {
      int u=runif(0,1); 
      if (u<=2/3) R[j*r+s]=0;
      if (u>2/3 & u<=5/6) R[j*r+s]=sqrt(3);
      if (u>5/6) R[j*r+s]=-sqrt(3);
    }
  }
  PutRNGstate ();
  zero_vect(Ar,n*r);
  for(int i=0; i<n; i++) {
    for(int index_j=0; index_j<A_byrow->nb_col[i]; index_j++) {
      int j      = A_byrow->col_index[i][index_j];
      double xij = A_byrow->values[i][index_j];
      for (int s=0; s<r; s++) {
        Ar[i*r+s]=Ar[i*r+s]+xij*R[j*r+s];
      }
    }
  }
  free(R);
  UNPROTECT(7);  
  return R_NilValue;
}

// --------------------------------------------------------------------------------------------------------
