\name{Rcoclust-package}
\alias{Rcoclust-package}
\alias{Rcoclust}
\docType{package}
\title{
\packageTitle{Rcoclust}
}
\description{
\packageDescription{Rcoclust}
}
\details{
This package includes several functions for co-clustering (see references):

1) CoClus,

2) Diagonal Double K-means, 

3) Fuzzy Diagonal Double K-means.

The main functions are wrappers to C functions in a shared dynamic library. A relevant initialisation may be wanted (see references). A common parameter is named envrdata for the input matrix in a sparse format. A function for reducing a dimension of a sparse matrix via a random projection is included. For any query about the package, thanks for sending an electronic mail.
}
\author{
\packageAuthor{Rcoclust}

Maintainer: \packageMaintainer{Rcoclust}
}

\references{
 M. Ailem, F. Role, M. Nadif, Co-clustering Document-term Matrices by Direct Maximization of Graph Modularity. CIKM, 2015, 1807-1810.

 E. Bingham, H. Mannila, Random projection in dimensionality reduction: applications to image and text data. ACM KDD, NY, 2001, 245-250.

 G. Govaert, M. Nadif, Co-Clustering, December 2013, Wiley-ISTE, 256 pages.
 
 T. Hofmann Probabilistic latent semantic indexing, SIGIR, 1999, 50-57.
 
 C. Laclau, M. Nadif, Hard and fuzzy diagonal co-clustering for document-term partitioning, Neurocomputing 193, 2016, 133-147.

 R. Priam, M. Nadif, G. Govaert, Generalized topographic block model, Neurocomputing 173, 2016, 442-449.

 W.C. Tjhi, L. Chen, A partitioning-based algorithm to fuzzy co-cluster documents and words, Pattern Recognition Lett. 27, 2006, 151-159.
}
