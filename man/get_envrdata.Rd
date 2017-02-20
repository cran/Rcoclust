\name{get_envrdata}
\alias{get_envrdata}
\alias{A_ijx}
\alias{lbs}
\alias{name}
\title{Construction of the data for the co-clustering procedures}
\description{Function which takes as en entry a matrix A_ijx in a sparse format for the rows, and adds in a variable of type environment: the data vectors in format coordinate list for the columns while aggregating the number of rows (n), columns (d), non null elements in total (nnz), per rows (nnzi) and per columns (nnzj), plus the name and the true labels (lbs) of the corresponding dataset.}
\usage{get_envrdata(A_ijx,lbs,name,datacol)}
\arguments{
  \item{A_ijx}{an object of class matrix with three columns corresponding to the format called coordinate list (COO)  (row, column, value).}
  \item{lbs}{integer vector with the true class labels for the rows.}
  \item{name}{char array with the name of the dataset.}
  \item{datacol}{integer with value 0 or 1 if it is added the three vectors for the column-wise data.}
}
\value{
  As a result an object of class environment which contains the data for the co-clustering procedures.
}
\examples{
  library(Rcoclust);
  
  #load data
  data(data_cstr);
  
  #library(Matrix); 
  #A = sparseMatrix(i=A_ijx[,1]+1,j=A_ijx[,2]+1,
  #                  x=A_ijx[,3]);
  envrdata = get_envrdata(A_ijx,lbs,name,1);
  
  #retrieve matrix size and number of classes
  n=envrdata$n;
  d=envrdata$d;
  g=length(unique(envrdata$lbs));
}
