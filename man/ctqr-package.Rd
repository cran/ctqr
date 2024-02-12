\name{ctqr-package}
\alias{ctqr-package}
\docType{package}
\title{Censored and Truncated Quantile Regression}
\description{
Fit quantile regression models to survival data, allowing for right censoring, left truncation, and interval censoring.
}
\details{
\tabular{ll}{
Package: \tab ctqr\cr
Type: \tab Package\cr
Version: \tab 2.1\cr
Date: \tab 2024-02-12\cr
License: \tab GPL-2\cr
}
The main function \code{\link{ctqr}} is used for model fitting. Other documented functions are \code{\link{predict.ctqr}}, to obtain prediction from a \code{ctqr} object, \code{\link{plot.ctqr}}, to plot quantile regression coefficients, and \code{\link{ctqr.control}}, that can be used to set the operational parameters for the estimation algorithm.
}
\author{
Paolo Frumento

Maintainer: Paolo Frumento <paolo.frumento@unipi.it>
}
\references{
Frumento, P., and Bottai, M. (2017). An estimating equation for censored and truncated quantile regression. \emph{Computational Statistics and Data Analysis}, Vol.113, pp.53-63. ISSN: 0167-9473.

Frumento, P. (2022). A quantile regression estimator for interval-censored data. \emph{The International Journal of Biostatistics}, 19 (1), pp. 81-96.
}
\keyword{ package }
\seealso{
\code{\link{pchreg}}, that is used to compute a preliminary estimate of the conditional outcome distribution.
}
