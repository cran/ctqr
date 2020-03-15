\name{plot.ctqr}
\alias{plot.ctqr}
\title{
Plot Quantile Regression Coefficients
}
\description{
Plots quantile regression coefficients 
\eqn{\beta(p)} as a function of \eqn{p},
based on a fitted model of class \dQuote{\code{ctqr}}.
}
\usage{
\method{plot}{ctqr}(x, which = NULL, ask = TRUE, \ldots)
}
\arguments{
  \item{x}{
an object of class \dQuote{\code{ctqr}}.
}
  \item{which}{
  an optional numerical vector indicating which coefficient(s) to plot. If \kbd{which = NULL}, 
  all coefficients are plotted.
}
\item{ask}{logical. If \kbd{which = NULL} and \kbd{ask = TRUE} (the default), you will be asked
interactively which coefficients to plot.}
\item{\ldots}{additional graphical parameters, that can include \kbd{xlim, ylim, xlab, ylab, col, lwd}. 
See \code{\link{par}}.}
}
\details{
With this command, a plot of \eqn{\beta(p)} versus \eqn{p} is created, provided that at least
two quantiles have been estimated. Dashed lines represent 95\% confidence intervals, while the horizontal dotted line indicates the zero.
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\seealso{
\code{\link{ctqr}}
}
\examples{
  
  # using simulated data
  
  n <- 1000
  x <- runif(n)
  t <- 1 + x + rexp(n)
  c <- runif(n, 1,10)
  y <- pmin(c,t)
  d <- (t <= c)
  
  par(mfrow = c(1,2))
  plot(ctqr(Surv(y,d) ~ x, p = seq(0.05,0.95,0.05)), ask = FALSE) 
}
\keyword{methods}

