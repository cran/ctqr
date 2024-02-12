\name{internals}
\alias{internals}

\alias{findagoodestimator}
\alias{safe.pch}
\alias{quickpred}
\alias{check.in.ctqr}
\alias{start}
\alias{ee.u}
\alias{ee.cens}
\alias{ee.cens.trunc}
\alias{ee.icens}
\alias{qr_gs}

\alias{nobs.ctqr}
\alias{vcov.ctqr}
\alias{coef.ctqr}
\alias{print.ctqr}
\alias{summary.ctqr}
\alias{print.summary.ctqr}

\alias{firstep.ct}
\alias{firstep.ic}
\alias{asy.qr.ct}
\alias{asy.qr.ic}

\alias{safesolve}
\alias{regularize_corr}

\title{Internal Functions}
\description{
Functions for internal use only.
}
\usage{
findagoodestimator(dat,w,type)
safe.pch(obj, eps = 1e-6)
quickpred(obj,y)
check.in.ctqr(z,y,d,x,off,weights, breaks)
start(CDF, p, x, y, off, weights) 
ee.u(beta, tau, CDF, V)
ee.cens(beta, tau, CDF, V)
ee.cens.trunc(beta, tau, CDF, V)
ee.icens(beta, tau, CDF, V)
qr_gs(beta0, check, p, CDF, a = 0.5, b = 1.5, maxit = 1000, tol = 1e-6, type)

\method{nobs}{ctqr}(object, \ldots)
\method{vcov}{ctqr}(object, \ldots)
\method{coef}{ctqr}(object, \ldots)
\method{print}{ctqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)
\method{summary}{ctqr}(object, \ldots)
\method{print.summary}{ctqr}(x, digits = max(3L, getOption("digits") - 3L), \ldots)

firstep.ct(obj, z,y,d,w)
firstep.ic(obj, y)
asy.qr.ct(z,y,d,x, weights, p, fit, fit.ok)
asy.qr.ic(z,y,d,x, weights, p, fit, fit.ok)

safesolve(A)
regularize_corr(A, c = 0.9)
}
\keyword{internal}


