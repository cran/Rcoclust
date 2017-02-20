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
#include <R_ext/Rdynload.h>
#include "rcoclust.h"

static const R_CallMethodDef callFuncs[]  = {
  {"C_coclus",  (DL_FUNC) &coclus, 8},
  {"C_ddkm",    (DL_FUNC) &ddkm,   11},
  {"C_randomp", (DL_FUNC) &randomp, 7},
  {NULL, NULL, 0}
};

void R_init_Rcoclust(DllInfo *dll) {
  R_registerRoutines(dll, NULL, callFuncs, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
