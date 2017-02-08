\name{cc_coclus}
\alias{cc_coclus}
\title{CoClus method}
\description{CoClus method via a graph modularity maximization.}
\usage{cc_coclus(g,envrdata,zi_coclus,wj_coclus,transfrm,maxiter,debug)}
\arguments{
  \item{g}{number of clusters for rows and columns.}
  \item{envrdata}{environment with data.}
  \item{zi_coclus}{row clusters.}
  \item{wj_coclus}{column clusters.}
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
  data(data_news4);
  
  #retrieve matrix size and number of classes
  n=envrdata$n;
  d=envrdata$d;
  g=length(unique(envrdata$lbs));
  
  #coclus
  bestresu=NULL;
  for (m in 1:100) {
    zi_coclus=as.integer( sample(x = 1:g,size = n,replace = TRUE) );
    wj_coclus=as.integer( sample(x = 1:g,size = d,replace = TRUE) );
    resu=Rcoclust::cc_coclus(g,envrdata,zi_coclus,wj_coclus,3,80,0);
    if (m==1) bestresu=resu;
    if (m>1) {
      if (resu$obj[length(resu$obj)]>bestresu$obj[length(bestresu$obj)])
        bestresu=resu;
    }
  }
  print(table(envrdata$lbs,bestresu$zi));
}
