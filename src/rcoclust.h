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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

// ----------------------------------------------------------------------------

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

double fun_transfrm(double , int , int , int );
void read_data_byrow(struct mat_byrow *, const int *, const int *, const double *, const int *, const int *, const int , const int , const int , const int );
void read_data_bycol(struct mat_bycol *, const int *, const int *, const double *, const int *, const int *, const int , const int , const int , const int );
void free_data_byrow(struct mat_byrow *, const int );
void free_data_bycol(struct mat_bycol *, const int );

void zero_vect(double *, const int );
void zero_vect_int(int *, const int );
void copy_vect(double *, const double *, const int );

void bound_values(double *, int , double );
void compute_zi(const double *, int *, const int , const int , const double );
void compute_sums_FDDKM(const double* ,const double* , double* , double* , const double , const double , const int , const int , const int , const double );
void compute_sums_DDKM(const int* ,const int* , int* , int* , const int , const int , const int );
void vect_random(int typer, double , double* , int );

SEXP coclus(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
SEXP ddkm(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
SEXP fddkm(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
SEXP randomp(SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , SEXP );
