\name{contr.Tukey}
\alias{contr.Tukey}
\title{
Pairwise contrast matrix
}
\description{
Returns a matrix of Tukey-type contrasts (all-pairwise)
}
\usage{
contr.Tukey(n, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{an integer: the number of levels, corresponding to the number of columns in the output matrix}
  \item{names}{a vector of names with same length as n}
}
\details{
This functions is used for creating contrast matrices for use with glht or other methods for fitting linear contrasts. The rows of the resulting matrices contain the coefficients of contrasts relating to a factor with n levels. The names of the levels can be optionally given with the argument 'names', other wise they default to the numbers from 1 to n.
}
\value{
A matrix with n columns and n(n-1)/2 rows
}
\author{
Andrea Onofri
}
\examples{
contr.Tukey(8, LETTERS[1:8])
}
