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

void zero_vect(double *v, const int dim) { for (int l=0; l<dim; l++) v[l]=0+1E-8; }
void zero_vect_int(int *v, const int dim) { for (int l=0; l<dim; l++) v[l]=0; }
void copy_vect(double *v_out, const double *v_in, const int dim) {for (int l=0; l<dim; l++) v_out[l]=v_in[l];}

void bound_values(double *v, int n, double tol) {
  for (int i=0; i<n; i++) {
    if (v[i]>1-tol) v[i]=1-tol;
    if (v[i]<tol) v[i]=tol;
  }
}

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

void compute_sums_FDDKM(const double* zik,const double* wjk, double* wkbeta, double* zkbeta, const double alpha, const double beta, const int n, const int d, const int g, const double tol) {
  zero_vect(zkbeta,g);
  zero_vect(wkbeta,g);
  for(int k=0;k<g;k++) { for(int j=0;j<d;j++) {wkbeta[k]=wkbeta[k]+pow(wjk[j*g+k],beta);} wkbeta[k]=wkbeta[k]+tol; }
  for(int k=0;k<g;k++) { for(int i=0;i<n;i++) {zkbeta[k]=zkbeta[k]+pow(zik[i*g+k],alpha);} zkbeta[k]=zkbeta[k]+tol; }
}

void compute_sums_DDKM(const int* zi,const int* wj, int* wk, int* zk, const int n, const int d, const int g)
{
  zero_vect_int(zk,g); for(int k=0;k<g;k++) for(int i=0;i<n;i++) zk[k]=zk[k]+(int)(zi[i]==k);
  zero_vect_int(wk,g); for(int k=0;k<g;k++) for(int j=0;j<d;j++) wk[k]=wk[k]+(int)(wj[j]==k);
}

void vect_random(int typer, double sgr, double* R_s, int d) {
  //zero_vect(R_s,d);
  GetRNGstate ();
  for (int j=0; j<d; j++) {
    if (typer==0) {
      R_s[j]=rnorm (0.0 ,sgr);
    } else if (typer==1) {
      double u=runif(0,1);
      if (u<=0.5) {R_s[j]=-1;} else{R_s[j]=+1;} 
    } else {
      double u=runif(0,1); 
      if (u<=2/3) R_s[j]=0;
      if (u>2/3 & u<=5/6) R_s[j]=1;//\sqrt(3);
      if (u>5/6) R_s[j]=-1;//\sqrt(3);
    }
  }
  PutRNGstate ();
}
