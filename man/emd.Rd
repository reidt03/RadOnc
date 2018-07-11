\name{emd}
\alias{emd}
\alias{emd2d}
\alias{emdw}
\alias{emdr}
\title{
Earth Mover's Distance
}
\description{
\code{emd} computes Earth Mover's Distance (related to 1st Mallows and
Wasserstein distances) between distributions. \code{emd} and \code{emdw}
use (weight,location) notation whereas \code{emd2d} compares two
distributions represented as matrices over a grid.
}
\usage{
emd(A, B, dist="euclidean", ...)
emdw(A, wA, B, wB, dist="euclidean", ...)
emd2d(A, B, xdist = 1, ydist = 1, dist="euclidean", ...)
emdr(A, B, extrapolate=NA, flows=FALSE, dist="euclidean", ...)
}
\arguments{
  \item{A}{matrix A}
  \item{B}{matrix B}
  \item{extrapolate}{if set to 1 or 2 the mass of A or B respectively is
    used to extrapolate the distance by penalization using the mass
    quotient assuming the other signature is truncated and thus more
    unlikely to match.
    It has any effect only if the other specified signature has larger
    mass.}
  \item{flows}{logical indicating whether flows should be returned
    in the \code{"flows"} attribute of the result.}
  \item{wA}{weights for locations specified by A}
  \item{wB}{weights for locations specified by B}
  \item{xdist}{
    distance between columns (scalar) or a vector of positions of the columns
  }
  \item{ydist}{
    distance between rows (scalar) of a vector of positions of the rows
  }
  \item{dist}{distance to be used for the computation of the cost over the
    locations. Must be either \code{"euclidean"}, \code{"manhattan"} or
    a closure taking two vectors and returning a scalar number. The
    latter case is much less efficient because it requires R evaluation
    for every possible combination of flows.}
  \item{...}{additional parameters for future use, currently unused}
}
\details{
  \code{emd2d} interprets the two matrices \code{A} and \code{B} as
  a distibution over a two-dimensional grid. The distance between the
  grid points in each direction is defined by \code{xdist} and
  \code{ydist}. Both matrices must have the same dimensionality.

  \code{emd} uses first column of each matrix as the weigths and the
  remaining columns as location coordinates in a up to four-dimensional
  space. \code{A} and \code{B} must have the same number of columns.
  
  \code{emdw} separates the weights from the location matrices but is
  otherwise identical to \code{emd}.
  
  \code{emdr} uses the original EMD implementation by Yossi Rubner
  from Stanford. In case \code{A} and \code{B} are not densities, the
  weighted sum of flows is normalized by the smaller total mass of the
  two. The version of the \code{emd} package released on CRAN
  contains only this implementation and all other functions are just
  front-ends for the call to \code{emdr}.
}
\value{
  Earth Mover's Distance between of the distributions \code{A} and \code{B}.
  If \code{A} and \code{B} are not distributions then \code{A} is the
  source and \code{B} is the target.
}
\author{
  R code by Simon Urbanek, EMD code by Yossi Rubner.
}
\note{
  This is an open-source version of the package which contains only the
  implementation by Yossi Rubner.
}
\examples{
A <- matrix(1:6 / sum(1:6), 2)
B <- matrix(c(0, 0, 0, 0, 0, 1), 2)
emd2d(A, B)
# if we bring the rows closer, the distance will be reduced
# since mass from the first row has to move down
emd2d(A, B,, 0.1)

# use Manhattan distance instead
emd2d(A, B, dist="manhattan")
# same, but using R-side closure
emd2d(A, B, dist=function(x, y) sum(abs(x - y)))

# the positions can overlap - this is a degenerate case of that
emd2d(A, B, rep(0, 3), rep(0, 2))
# and just a sanity check
emd2d(A, A) + emd2d(B, B)

# and the weight/location code should, hopefully have the same results
A. <- matrix(c(1:6 / sum(1:6), 1,2,1,2,1,2, 1,1,2,2,3,3), 6)
B. <- matrix(c(1, 2, 3), 1)
stopifnot(emd(A., B.) == emd2d(A, B))
stopifnot(emd(A., B.) == emdw(A.[,-1], A.[,1], B.[,-1,drop=FALSE], B.[,1]))
}
\keyword{multivariate}
