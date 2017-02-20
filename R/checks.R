#/*==========================================================================
#  Copyright (C) R. Priam 
#Date Version 0.1.0 : 2017-01-31
#Date last update   : 2017-02-19
#==========================================================================
#  This file is part of Rcoclust.
#
#Rcoclust is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#Rcoclust is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with Rcoclust.  If not, see <http://www.gnu.org/licenses/>.
#========================================================================== */

.check_dataI <- function(envrdata) {
  if (!exists("envrdata")) stop("envrdata does not exist.")
  if (is.null("envrdata")) stop("envrdata is null.")
  if (is.null("envrdata$n")) stop("envrdata$n is null")
  if (is.null("envrdata$d")) stop("envrdata$d is null.")
  if (is.null("envrdata$name")) stop("envrdata$name is null")
  if (is.null("envrdata$Ais_byrow")) stop("envrdata$Ais_byrow is null.")
  if (is.null("envrdata$Ajs_byrow")) stop("envrdata$Ajs_byrow is null.")
  if (is.null("envrdata$Axs_byrow")) stop("envrdata$Axs_byrow is null.")
  if (is.null("envrdata$nnzi")) stop("envrdata$nnzi is null.")
  if (is.null("envrdata$nnzj")) stop("envrdata$nnzj is null.")
  if (is.null("envrdata$nnz")) stop("envrdata$nnz is null.")
  if (length(envrdata$Ais_byrow)!=envrdata$nnz) stop("length(envrdata$Ais_byrow)!=envrdata$nnz")
  if (length(envrdata$Ajs_byrow)!=envrdata$nnz) stop("length(envrdata$Ajs_byrow)!=envrdata$nnz")
  if (length(envrdata$Axs_byrow)!=envrdata$nnz) stop("length(envrdata$Axs_byrow)!=envrdata$nnz")
  if (max(envrdata$Ais_byrow)>=envrdata$n) stop("Ais_byrow elements must be in [0;n-1]")
  if (min(envrdata$Ais_byrow)<0) stop("Ais_byrow elements must be in [0;n-1]")
  if (max(envrdata$Ajs_byrow)>=envrdata$d) stop("Ajs_byrow elements must be in [0;d-1]")
  if (min(envrdata$Ajs_byrow)<0) stop("Ajs_byrow elements must be in [0;d-1]")
  if (min(envrdata$Axs_byrow)<0) stop("Axs_byrow elements must be >=0")
  if (sum(envrdata$nnzi)!=envrdata$nnz) stop("sum(envrdata$nnzi)!=envrdata$nnz")
  if (sum(envrdata$nnzj)!=envrdata$nnz) stop("sum(envrdata$nnzj)!=envrdata$nnz")
}

.check_dataJ <- function(envrdata) {
  if (!exists("envrdata")) stop("envrdata does not exist.")
  if (is.null("envrdata")) stop("envrdata is null.")
  if (is.null("envrdata$n")) stop("envrdata$n is null.")
  if (is.null("envrdata$d")) stop("envrdata$d is null.")
  if (is.null("envrdata$name")) stop("envrdata$name is null.")
  if (is.null("envrdata$Ais_bycol")) stop("envrdata$Ais_bycol is null.")
  if (is.null("envrdata$Ajs_bycol")) stop("envrdata$Ajs_bycol is null.")
  if (is.null("envrdata$Axs_bycol")) stop("envrdata$Axs_bycol is null.")
  if (is.null("envrdata$nnzi")) stop("envrdata$nnzi is null.")
  if (is.null("envrdata$nnzj")) stop("envrdata$nnzj is null.")
  if (is.null("envrdata$nnz")) stop("envrdata$nnz is null.")
  if (length(envrdata$Ais_bycol)!=envrdata$nnz) stop("length(envrdata$Ais_bycol)!=envrdata$nnz")
  if (length(envrdata$Ajs_bycol)!=envrdata$nnz) stop("length(envrdata$Ajs_bycol)!=envrdata$nnz")
  if (length(envrdata$Axs_bycol)!=envrdata$nnz) stop("length(envrdata$Axs_bycol)!=envrdata$nnz")
  if (max(envrdata$Ais_bycol)>=envrdata$n) stop("Ais_bycol elements must be in [0;n-1]")
  if (min(envrdata$Ais_bycol)<0) stop("Ais_bycol elements must be in [0;n-1]")
  if (max(envrdata$Ajs_bycol)>=envrdata$d) stop("Ajs_bycol elements must be in [0;d-1]")
  if (min(envrdata$Ajs_bycol)<0) stop("Ajs_bycol elements must be in [0;d-1]")
  if (min(envrdata$Axs_bycol)<0) stop("Axs_bycol elements must be >=0")
  if (sum(envrdata$nnzi)!=envrdata$nnz) stop("sum(envrdata$nnzi)!=envrdata$nnz")
  if (sum(envrdata$nnzj)!=envrdata$nnz) stop("sum(envrdata$nnzj)!=envrdata$nnz")
}

.check_params_randp <- function(envrdata,dimr,sgr,vect_Ar,transfrm,debug) {
  if (dimr<1) stop("dimr must be >=1");
  if (sgr<=0) stop("sgr must be >0");
  if (!transfrm%in%0:3) stop("transfrm must be in {0,1,2,3}");
  if (!debug%in%0:1) stop("debug must be in {0,1}");
  if (is.null(vect_Ar)) stop("vect_Ar is null");
  if (length(as.numeric(vect_Ar))!=envrdata$n*dimr) stop("length(as.numeric(vect_Ar))!=envrdata$n*dimr");
}

.check_params_coclus <- function (g,envrdata,zi,wj,transfrm,maxiter,debug) {
  if (!g%in%2:ceiling(envrdata$n/10)) stop("g must be integer in {2, 3, ..., [n/10]}");
  if (!transfrm%in%0:3) stop("transfrm must be in {0,1,2,3}");
  if (!maxiter%in%10:80) stop("maxiter must be in {10,11, ..., 80}");
  if (!debug%in%0:1) stop("debug must be in {0,1}");
  if (is.null(zi)) stop("zi is null");
  if (length(as.numeric(zi))!=envrdata$n) stop("length(as.numeric(zi))!=envrdata$n");
  if (is.null(wj)) stop("wj is null");
  if (length(as.numeric(wj))!=envrdata$d) stop("length(as.numeric(wj))!=envrdata$d");
}

.check_params_ddkm <- function (g,envrdata,zi,wj,delta,transfrm,maxiter,debug) {
  .check_params_coclus(g,envrdata,zi,wj,transfrm,maxiter,debug);
  if (!(delta==-1 || (delta>=min(envrdata$Axs_byrow) && delta<=max(envrdata$Axs_byrow))))
    stop("delta must be either =-1 for automatic computation of max({xij}) either in [min({xij});max({xij)}]");
}

.check_params_get_envrdata <- function(A_ijx,lbs,name,datacol) { 
  if (is.null(A_ijx) || is.null(lbs) || is.null(name) || is.null(datacol)) stop("null parameter");
  if (!datacol%in%0:1) stop("datacol must be in {0,1}");
  if (class(A_ijx)!="matrix") stop("A_ijx must be a matrix");
  if (ncol(A_ijx)!=3) stop("A_ijx must have 3 columns");
  if (class(lbs)!="numeric" && class(lbs)!="integer") stop("lbs must be a vector");
  if (max(A_ijx[,1])+1!=length(lbs)) stop("max(A_ijx[,1])+1!=length(lbs)");
}
  
.check_params_fddkm <- function (g,envrdata,zi,wj,vect_cik,vect_wjk,alpha,beta,delta,transfrm,maxiter,debug) {
  .check_params_ddkm(g,envrdata,zi,wj,delta,transfrm,maxiter,debug);
  if (alpha <= 1) stop("alpha must be >1");
  if (beta <= 1) stop("beta must be >1");
  if (is.null(vect_cik)) stop("vect_cik is null");
  if (length(as.numeric(vect_cik))!=envrdata$n*g) stop("length(as.numeric(vect_cik))!=envrdata$n*g");
  if (is.null(vect_wjk)) stop("vect_wjk is null");
  if (length(as.numeric(vect_wjk))!=envrdata$d*g) stop("length(as.numeric(vect_wjk))!=envrdata$d*g");
}
