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

cc_coclus <- function(g,envrdata,zi,wj,transfrm,maxiter,debug) {
  .check_dataI(envrdata);
  #.check_dataJ(envrdata);
  .check_params_coclus(g,envrdata,zi,wj,transfrm,maxiter,debug);

  ans1 <- .Call("C_coclus",
                as.integer(envrdata$Ais_byrow), as.integer(envrdata$Ajs_byrow),
                as.numeric(envrdata$Axs_byrow), as.integer(envrdata$nnzi), 
                as.integer(envrdata$nnzj), as.integer(zi), as.integer(wj), 
                as.numeric(c(envrdata$n, envrdata$d, g, envrdata$nnz, transfrm, maxiter, debug,1E-7)),
                PACKAGE="Rcoclust");
  
  return (list(obj=ans1, zi=zi, wj=wj));
}

cc_ddkm <- function(g,envrdata,zi,wj,delta,transfrm,maxiter,debug) {
  .check_dataI(envrdata);
  .check_dataJ(envrdata);
  .check_params_ddkm(g,envrdata,zi,wj,delta,transfrm,maxiter,debug);
  
  ans1 <- .Call("C_ddkm",
                as.integer(envrdata$Ais_byrow), as.integer(envrdata$Ajs_byrow),
                as.numeric(envrdata$Axs_byrow), as.integer(envrdata$nnzi), 
                as.integer(envrdata$Ais_bycol), as.integer(envrdata$Ajs_bycol),
                as.numeric(envrdata$Axs_bycol), as.integer(envrdata$nnzj), 
                as.integer(zi), as.integer(wj),
                as.numeric(c(envrdata$n, envrdata$d, g, envrdata$nnz, delta, transfrm, maxiter, debug,1E-7)),
                PACKAGE="Rcoclust");
  
  return (list(obj=ans1, zi=zi, wj=wj));
}

randp <- function(envrdata,dimr,sgr,vect_Ar,transfrm,debug) {
  .check_dataI(envrdata);
  .check_params_randp(envrdata,dimr,sgr,vect_Ar,transfrm,debug);
  
  ans1 <- .Call("C_randomp",
                as.integer(envrdata$Ais_byrow), as.integer(envrdata$Ajs_byrow),
                as.numeric(envrdata$Axs_byrow), as.integer(envrdata$nnzi), as.integer(envrdata$nnzj),
                as.numeric(vect_Ar),
                as.numeric(c(envrdata$n,envrdata$d,envrdata$nnz,dimr,sgr,0,transfrm,debug)),
                PACKAGE="Rcoclust");
}

get_envrdata <- function(A_ijx,lbs,name,datacol) { 
  .check_params_get_envrdata(A_ijx,lbs,name,datacol);
  
  n    = max(A_ijx[,1])+1;
  d    = max(A_ijx[,2])+1;
  nnzi = as.numeric(table(sort(A_ijx[,1])));
  nnzj = as.numeric(table(sort(A_ijx[,2])));
  nnz  = nrow(A_ijx);
  
  if (datacol==1) {
    tA_ijx = data.frame(Ais_bycol=rep(0:(n-1),nnzi),
                        Ajs_bycol=A_ijx[,2],
                        Axs_bycol=A_ijx[,3]);
    tA_ijx = tA_ijx[order(tA_ijx$Ajs_bycol),];
  }
  
  envrdata = new.env(parent=emptyenv());
  assign("name",name,envrdata);#
  assign("n",n,envrdata);#
  assign("d",d,envrdata);#
  assign("Ais_byrow", A_ijx[,1],envrdata);#0,...,d-1
  assign("Ajs_byrow", A_ijx[,2],envrdata);#
  assign("Axs_byrow", A_ijx[,3],envrdata);#
  if (datacol==1) assign("Ais_bycol",tA_ijx[,1],envrdata);#0,...,d-1
  if (datacol==1) assign("Ajs_bycol",tA_ijx[,2],envrdata);#
  if (datacol==1) assign("Axs_bycol",tA_ijx[,3],envrdata);#
  assign("nnzi",nnzi,envrdata);#
  assign("nnzj",nnzj,envrdata);#
  assign("nnz",nrow(A_ijx),envrdata);#
  assign("lbs",lbs,envrdata);#
  
  return(envrdata);
}
