\name{randp}
\alias{randp}
\title{Random projection for sparse matrix}
\description{Reduced matrix after right multiplication with a generated Gaussian random matrix.}
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
  The function alters the matrix in vect_Ar. It returns 0 if the function parameters are detected as erroneous, otherwise 1.
}
\examples{
  library(Rcoclust);
  
  #load data
  data(data_cstr);
  
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
