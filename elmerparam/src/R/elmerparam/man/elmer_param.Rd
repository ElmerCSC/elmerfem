\name{elmer_param}
\alias{elmer_param}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to provide an interface to ElmerParam}
\description{
Uses ElmerParam to define a function of a bunch of parameters, using
external computation programs such as ElmerSolver.
}
\usage{
elmer_param(xr = NULL, xi = NULL, tag = "", nfun = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xr}{Array of real parameters (R1\ldots{}R$n$)}
  \item{xi}{Array of integer parameters (I1\ldots{}I$n$)}
  \item{tag}{Unique tag, e.g. for creation of unique file names. (T)}
  \item{nfun}{Length of result array}
}
\details{
}
\value{
A vector of output values (O1\ldots{}O\code{nfun}).
}
\references{The ElmerParam Manual}
\author{Erik Edelmann <ErikEdelmann@csc.fi>}
\note{}
\seealso{}
\examples{
}
\keyword{Elmer}% at least one, from doc/KEYWORDS
\keyword{ElmerParam}% __ONLY ONE__ keyword per line
