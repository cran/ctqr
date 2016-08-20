\name{ctqr.control}
\alias{ctqr.control}
\title{
Auxiliary Function for Root Search
}
\description{
This functions can be used within a call to \code{\link{ctqr}}, 
to control the operational parameters of the root search algorithm.
}
\usage{
ctqr.control(tol = 1e-06, maxit = 1000, a = 0.5, b = 1.25)
}
\arguments{
  \item{tol}{
	positive convergence tolerance: the algorithm stops when the maximum 
  absolute change between two consecutive estimates is smaller than \kbd{tol}.
  }
  \item{maxit}{maximum number of iterations.}
  \item{a,b}{
    numeric scalar with \kbd{0 < a < 1} and \kbd{b > 1}. See \sQuote{Details}.
  }
}
\details{
  For a current estimate \kbd{beta}, a new estimate is computed as 
  \kbd{beta_new = beta + delta*s(beta)}, where \kbd{s(beta)} is the current 
  value of the objective function and \kbd{delta} is a positive multiplier. 
  If \kbd{sum(s(beta_new)^2) < sum(s(beta)^2)}, the iteration is accepted and 
 \kbd{delta} is multiplied by \kbd{b}. Otherwise, \kbd{beta_new} is rejected and 
  \kbd{delta} is multiplied by \kbd{a}. By default, \kbd{a = 0.5} and \kbd{b = 1.25}. 
  Choosing \kbd{a,b} closer to \kbd{1} may result in a more accurate estimate, 
  but will require a larger number of iterations.
}
\value{
The function returns its arguments. If some was not correctly specified, it is set to its default and a warning message is returned.
}
\seealso{
\code{\link{ctqr}}
}
\keyword{utilities}
