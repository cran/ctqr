\name{ctqr}
\alias{ctqr}
\title{
Censored and Truncated Quantile Regression
}
\description{
Fits a quantile regression model to possibly censored and truncated data, e.g., survival data.
}
\usage{
ctqr(formula, data, weights, p = 0.5, CDF, control = ctqr.control(), \ldots)
}
\arguments{
  \item{formula}{
an object of class \dQuote{\code{\link{formula}}}: a symbolic description of the regression model. 
The response must be a \kbd{Surv} object as returned by \code{\link{Surv}} (see \sQuote{Details}).
}
  \item{data}{
an optional data frame containing the variables in the model.
}
  \item{weights}{
an optional vector of weights to be used in the fitting process. The weights will always be normalized to sum to the sample size. This implies that, for example, using double weights will not halve the standard errors.
}
  \item{p}{
numerical vector indicating the order of the quantile(s) to be fitted.
}
  \item{CDF}{
an object of class \dQuote{\code{pch}}, i.e., the result of a call to \code{\link{pchreg}}. If missing, 
it will be computed internally with default settings. See \sQuote{Details}.
}
  \item{control}{
a list of operational parameters for the optimization algorithm, usually passed via \code{\link{ctqr.control}}.
}
\item{\ldots}{for future arguments.}
}
\details{
This function implements the method described in Frumento and Bottai (2017) for censored, truncated quantile regression, and Frumento (2022) for interval-censored quantile regression.

The left side of \kbd{formula} must be of the form \kbd{Surv(time, event)} if the data are right-censored, \kbd{Surv(time0, time, event)} if the data are right-censored and left-truncated (\kbd{time0 < time}, \code{time0} can be \kbd{-Inf}), and \kbd{Surv(time1, time2, type = "interval2")} if the data are interval-censored (use \kbd{time1 = time2} for exact observations, \kbd{time1 = -Inf} or \kbd{NA} for left-censored, and \kbd{time2 = Inf} or \kbd{NA} for right-censored). Using \kbd{Surv(time)} is also allowed and indicates that the data are neither censored nor truncated.

The conditional distribution function (\code{CDF}) of the response variable represents a nuisance parameter 
and is estimated preliminarly via \code{\link{pchreg}}. If missing, \code{CDF = pchreg(formula)} is used as default. See the \dQuote{Note} and the documentation of \command{pchreg}.

Estimation is carried out using an algorithm for gradient-based optimization. 
To estimate the asymptotic covariance matrix, standard two-step procedures are used
(e.g., Ackerberg et al., 2012).
}
\value{
An object of class \dQuote{\code{ctqr}}, which is a list with the following items:
\item{p}{the quantile(s) being estimated.}
\item{coefficients}{a named vector or matrix of quantile regression coefficients.}
\item{call}{the matched call.}
\item{n.it}{the number of iterations.}
\item{converged}{logical. The convergence status.}
\item{fitted}{the fitted values.}
\item{terms}{the \code{terms} object used.}
\item{mf}{the model frame used.}
\item{covar}{the estimated asymptotic covariance matrix.}
\item{CDF}{the used \code{CDF}.}

Note that the dimension of all items, except \code{call}, \code{terms}, \code{mf}, and \code{CDF},
is the same as the dimension of \kbd{p}. For example, if \code{p = c(0.25,0.5,0.75)}, \code{coefficients}
and \code{fitted} will be 3-columns matrices; \code{n.it} and \code{converged} will be vectors of 3 elements; 
and \code{covar} will be a list of three covariance matrices.

The generic accessor functions \code{summary}, \code{plot}, \code{predict}, \code{coef}, \code{terms}, \code{nobs}, 
can be used to extract information from the model. The functions
\code{waldtest} (from the package \pkg{lmtest}), and \code{linearHypothesis} (from the package \pkg{car}) can be 
used to perform Wald test, and to test linear restrictions. These functions, however, 
will only work if \code{p} is scalar.
}

\references{
Ackerberg, D., Chen, X., and Hahn, J. (2012). A practical asymptotic variance estimator 
for two-step semiparametric estimators. The Review of Economics and Statistics, 94 (2), 481-498.

Frumento, P., and Bottai, M. (2017). An estimating equation for censored and truncated quantile regression. \emph{Computational Statistics and Data Analysis}, Vol.113, pp.53-63. ISSN: 0167-9473.

Frumento, P. (2022). A quantile regression estimator for interval-censored data. \emph{The International Journal of Biostatistics}, 19 (1), pp. 81-96.
}

\note{
NOTE 1. The first-step estimator (the \kbd{CDF} argument) is computed using the \code{pchreg} function of the
\pkg{pch} package. To be correctly embedded in \kbd{ctqr}, a \kbd{pch} object must be constructed using  
the same observations, in the same order. 

If the first-step estimator is biased, and there is censoring or truncation,
the estimates of the quantile regression coefficients and their standard errors will also be biased. 

If the data are neither censored nor truncated, the \code{CDF} does not enter the estimating equation of the model. However, since the first-step estimator is used to compute the starting points, 
the final estimates may be sensitive to the supplied \code{CDF}.

NOTE 2. Right-censoring is a special case of interval censoring, in which exact events are identified by \code{time2 = time1}, while censored observations have \code{time2 = Inf}. Note, however, that \code{ctqr(Surv(time1, time2, type = "interval2") ~ x)} will not be identical to \code{ctqr(Surv(time = time1, event = (time2 < Inf)) ~ x)}. The estimating equation used for interval-censored data is that described in Frumento (2022), while that used for right-censored data is that of Frumento and Bottai (2017). The two estimating equations are only asymptotically equivalent.


}

\author{
Paolo Frumento <paolo.frumento@unipi.it>
}

\seealso{
\code{\link{plot.ctqr}}, \code{\link{predict.ctqr}}, \code{\link{pchreg}}
}
\examples{
# Using simulated data


# Example 1 - censored data ####################################################

n <- 1000
x1 <- runif(n); x2 <- runif(n)   # covariates
t <- runif(n, 0, 1 + x1 + x2)    # time variable (e.g., time to death)
c <- runif(n,0,5)                # censoring variable (e.g., end of follow-up)
y <- pmin(t,c)                   # observed variable = min(t,c)
d <- (t <= c)                    # 1 = event (e.g., death), 0 = censored

CDF1 <- pchreg(Surv(y,d) ~ x1 + x2)
model1 <- ctqr(Surv(y,d) ~ x1 + x2, p = 0.5, CDF = CDF1)
model2 <- ctqr(Surv(y,d) ~ x1, p = 0.5, CDF = CDF1)

# model1 is identical to ctqr(Surv(y,d) ~ x1 + x2, p = 0.5)
# model2 is NOT identical to ctqr(Surv(y,d) ~ x1, p = 0.5), 
  # which would have default CDF = pchreg(Surv(y,d) ~ x1)


# Example 2 - censored and truncated data ######################################

n <- 1000
x1 <- runif(n); x2 <- runif(n)   # covariates
t <- runif(n, 0, 1 + x1 + x2)    # time variable
c <- runif(n,0,5)                # censoring variable
y <- pmin(t,c)                   # observed variable = min(t,c)
d <- (t <= c)                    # 1 = event, 0 = censored

z <- rnorm(n) # truncation variable (e.g., time at enrollment)
w <- which(y > z) # data are only observed when y > z
z <- z[w]; y <- y[w]; d <- d[w]; x1 <- x1[w]; x2 <- x2[w]

# implement various CDFs and choose the model with smallest AIC

CDFs <- list(
  pchreg(Surv(z,y,d) ~ x1 + x2, breaks = 5),
  pchreg(Surv(z,y,d) ~ x1 + x2, breaks = 10),
  pchreg(Surv(z,y,d) ~ x1 + x2 + x1:x2, breaks = 5),
  pchreg(Surv(z,y,d) ~ x1 + x2 + x1^2 + x2^2, breaks = 10)
)

CDF <- CDFs[[which.min(sapply(CDFs, function(obj) AIC(obj)))]]
summary(ctqr(Surv(z,y,d) ~ x1 + x2, p = 0.5, CDF = CDF))


# Example 3 - interval-censored data ###########################################
# t is only known to be in the interval (t1,t2) ################################

n <- 1000
x1 <- runif(n); x2 <- runif(n)     # covariates
t <- runif(n, 0, 10*(1 + x1 + x2)) # time variable
t1 <- floor(t)                     # lower extreme of the interval
t2 <- ceiling(t)                   # upper extreme of the interval

model <- ctqr(Surv(t1,t2, type = "interval2") ~ x1 + x2, p = 0.5)
}

\keyword{survival}
\keyword{models}
\keyword{regression}