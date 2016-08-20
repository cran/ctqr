\name{internals}
\alias{internals}
\title{Internal Functions}
\description{
Functions for internal use only.
}
\usage{
findagoodestimator(dat,w)
safe.pch(obj, eps = 1e-6)
quickpred(obj,y)
check.in.ctqr(z,y,d,x,off,weights, breaks)
start(CDF, p, x, y, off, weights) 
ee.u(beta, tau, CDF, V)
ee.cens(beta, tau, CDF, V)
ee.cens.trunc(beta, tau, CDF, V)
qr.gs(beta0, check, p, CDF, a = 0.5, b = 1.5, maxit = 1000, tol = 1e-6, type)

\method{nobs}{ctqr}(object, \ldots)
\method{vcov}{ctqr}(object, \ldots)
\method{coef}{ctqr}(object, \ldots)
\method{print}{ctqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{summary}{ctqr}(object, \ldots)
\method{print.summary}{ctqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)

firstep(obj, z,y,d,w)
asy.qr(z,y,d,x, weights, p, CDF, fit, fit.ok)
}
\keyword{internal}


