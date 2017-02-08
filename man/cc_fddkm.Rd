\name{cc_fddkm}
\alias{cc_fddkm}
\title{Fuzzy Diagonal Double K-means}
\description{Fuzzy co-clustering via a diagonal reordering.}
\usage{cc_fddkm(g,envrdata,zi_fddkm,wj_fddkm,vect_cik_fddkm,vect_wjk_fddkm,
             alpha,beta,delta,transfrm,maxiter,debug)}
\arguments{
  \item{g}{number of clusters for rows and columns.}
  \item{envrdata}{environment with data.}
  \item{zi_fddkm}{row clusters.}
  \item{wj_fddkm}{column clusters.}
  \item{vect_cik_fddkm}{vector with the posteriors for the row clusters.}
  \item{vect_wjk_fddkm}{vector with the posteriors for the columns clusters.}
  \item{alpha}{coefficient for fuzzification per rows.}
  \item{beta}{coefficient for fuzzification per columns.}
  \item{delta}{common constant value per block.}
  \item{transfrm}{transformation of the data (0:none, 1:binarization, 2:tf-idf, 3:tf-idf+rows normalization).}
  \item{maxiter}{maximum number of iterations.}
  \item{debug}{flag for debug, if equal to 1 shows some informations to user.}
}
\value{
  The function returns a vector with the value of the objective function per iterations, and zi, wj which are the estimated cluster labels. 
}
\examples{
  library(Rcoclust);
  
  #load data
  data(data_news4);#data(data_ng5);
  
  #retrieve matrix size and number of classes
  n=envrdata$n;
  d=envrdata$d;
  g=length(unique(envrdata$lbs));
  
  #fddkm
  bestresu=NULL;
  for (m in 1:20) {
    zi_fddkm=as.integer( sample(x = 1:g,size = n,replace = TRUE) );
    wj_fddkm=as.integer( sample(x = 1:g,size = d,replace = TRUE) );
    vect_cik_fddkm=rep(0,n*g);
    vect_wjk_fddkm=rep(0,d*g);
    for (i in 1:n) vect_cik_fddkm[(i-1)*g+zi_fddkm[i]]=1;
    for (j in 1:d) vect_wjk_fddkm[(j-1)*g+wj_fddkm[j]]=1;
    resu=Rcoclust::cc_fddkm(g,envrdata,
                         zi_fddkm,wj_fddkm,vect_cik_fddkm,vect_wjk_fddkm,
                         1.002,1.002,-1,3,80,0);
    if (m==1) {bestresu=resu;} else {
      if (resu$obj[length(resu$obj)]<bestresu$obj[length(bestresu$obj)])
        bestresu=resu;
    }
  }
  print(table(envrdata$lbs,bestresu$zi));
}
