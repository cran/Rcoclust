\name{randp}
\alias{randp}
\alias{C_randomp}
\title{Random projection for sparse matrix}
\description{Implementation of the reduction of matrix via a right multiplication with a generated Gaussian random matrix. The data vectors for the columns are not required for the input parameters in envrdata.}
\usage{randp(envrdata,dimr,sgr,vect_Ar,transfrm,debug)}
\arguments{
  \item{envrdata}{environment with data.}
  \item{dimr}{reduced dimension for random projection.}
  \item{sgr}{standard deviation for random projection.}
  \item{vect_Ar}{vector with the transformed matrix after random projection.}
  \item{transfrm}{transformation of the data (0:none, 1:binarization, 2:tf-idf, 3:tf-idf+rows normalization).}
  \item{debug}{flag for debug, if equal to 1 shows some informations to user.}
}
\value{
  The function alters the matrix in vect_Ar with the result of the reduction.
}
\examples{
  library(Rcoclust);
  
  #load data
  data(data_cstr);
  envrdata=get_envrdata(A_ijx,lbs,name,0);
  
  #retrieve matrix size and number of classes
  n=envrdata$n;
  d=envrdata$d;
  g=length(unique(envrdata$lbs));
  
  #random projection
  dimr    = min(500,envrdata$d); #reduced dimension
  sgr     = 1;                   #std for projection
  vect_Ar = rep(0,n*dimr);       #reduced matrix in vector form
  Rcoclust::randp(envrdata,dimr,sgr,vect_Ar,3,0);
  Ar=matrix(vect_Ar,nrow=n,byrow = TRUE);
  
  #kmeans+random projection
  km=kmeans(Ar,g);
  table(km$cluster,envrdata$lbs);
}
