\name{ctqr-package}
\alias{ctqr-package}
\docType{package}
\title{Censored and Truncated Quantile Regression}
\description{
This package can be used to fit quantile regression models to survival data. 
The true conditional outcome distribution represents a nuisance parameter which is estimated preliminarly. 
Asymptotic theory of two-steps estimators is used to estimate the asymptotic covariance matrix.
}
\details{
\tabular{ll}{
Package: \tab ctqr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-02-18\cr
License: \tab GPL-2\cr
}
The main function \code{\link{ctqr}} is used for model fitting. 
Other documented functions are \code{\link{predict.ctqr}}, to obtain prediction from a \code{ctqr} object, 
\code{\link{plot.ctqr}}, to plot quantile regression coefficients, and \code{\link{ctqr.control}}, 
that can be used to set the operational parameters for the estimation algorithm.
}
\author{
Paolo Frumento

Maintainer: Paolo Frumento <paolo.frumento@ki.se>
}
\references{
Frumento, P., and Bottai, M. (2016). An estimating equation for censored and truncated quantile regression.
}
\keyword{ package }
\seealso{
\code{\link{pchreg}}, that is used to compute a preliminary estimate of the conditional outcome distribution.
}