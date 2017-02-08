\name{envrdata}
\docType{data}
\alias{envrdata}
\title{General format for the data matrix with a sparse structure, as an environment or a list}
\description{
  It is compound of the following fields.
  \itemize{
    \item name: name of the dataset.
    \item n: number of rows.
    \item d: number of columns.
    \item nnz: number of non zero (nz) cells.
    \item Ais_byrow: vector of nnz components with their row indexes when the matrix is read by rows.
    \item Ajs_byrow: vector of nnz components with their column indexes when the matrix is read by rows. 
    \item Axs_byrow: vector of nnz components with their values when the matrix is read by rows. 
    \item nnzi: vector of n components with the counts of nz cells per rows.  
    \item Ais_bycol: vector of nnz components with their row indexes when the matrix is read by columns.
    \item Ajs_bycol: vector of nnz components with their column indexes when the matrix is read by columns. 
    \item Axs_bycol: vector of nnz components with their values when the matrix is read by columns. 
    \item nnzj: vector of d components with the counts of nz cells per columns.
    \item lbs: vector of n components of the true (or from some algorithm) class labels when available.
  }
}
\keyword{datasets}
