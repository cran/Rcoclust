\name{cc_ddkm}
\alias{cc_ddkm}
\alias{C_ddkm}
\title{DDKM}
\description{Implementation with sparsity of the method which aims at a diagonal reordering.}
\usage{cc_ddkm(g,envrdata,zi,wj,delta,transfrm,maxiter,debug)}
\arguments{
  \item{g}{number of clusters for rows and columns.}
  \item{envrdata}{environment with data.}
  \item{zi}{row clusters.}
  \item{wj}{column clusters.}
  \item{delta}{common constant value per block.}
  \item{transfrm}{transformation of the data (0:none, 1:binarization, 2:tf-idf, 3:tf-idf+rows normalization).}
  \item{maxiter}{maximum number of iterations.}
  \item{debug}{flag for debug, if equal to 1 shows some informations to user.}
}
\value{
  The function alters the parameters zi and wj. It returns a vector with the value of the objective function per iterations, and zi, wj which are the estimated cluster labels.
}
\examples{
  library(Rcoclust);
  
  #load data
  data(data_news4);
  envrdata=get_envrdata(A_ijx,lbs,name,1);
  
  #retrieve matrix size and number of classes
  n=envrdata$n;
  d=envrdata$d;
  g=length(unique(envrdata$lbs));
  
  #ddkm
  bestresu=NULL;
  for (m in 1:100) {
    zi_ddkm=as.integer( sample(x = 1:g-1,size = n,replace = TRUE) );
    wj_ddkm=as.integer( sample(x = 1:g-1,size = d,replace = TRUE) );
    resu=Rcoclust::cc_ddkm(g,envrdata,zi_ddkm,wj_ddkm,-1,3,80,0);
    if (m==1) bestresu=resu;
    if (m>1) {
      if (resu$obj[length(resu$obj)]<bestresu$obj[length(bestresu$obj)])
        bestresu=resu;
    }
  }
  print(table(envrdata$lbs,bestresu$zi));
}
