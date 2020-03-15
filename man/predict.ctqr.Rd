\name{predict.ctqr}
\alias{predict.ctqr}
\title{
Prediction After Quantile Regression
}
\description{
This function returns predictions for an object of class \dQuote{\code{ctqr}}.
}
\usage{
\method{predict}{ctqr}(object, newdata, se.fit = FALSE, ...)
}
\arguments{
  \item{object}{
a \code{ctqr} object.
}
  \item{newdata}{
optional data frame in which to look for variables with which to predict. It must include all the covariates that enter the quantile regression model. If omitted, the fitted values are used.
}
  \item{se.fit}{
logical. If \kbd{TRUE}, standard errors of the predictions are also computed.
}
  \item{\dots}{
 for future methods.
}
}
\details{
This function produces predicted values obtained by evaluating the regression function at \code{newdata} (which defaults to \kbd{model.frame(object)}).

}
\value{
If \kbd{se = FALSE}, a matrix of fitted values, with rows corresponding to different observations, and one column for each value of \kbd{object$p}. If \kbd{se = TRUE}, a list with two items:
\item{fit}{a matrix of fitted values, as described above.}
\item{se.fit}{a matrix of estimated standard errors.}
}
\author{
Paolo Frumento \email{paolo.frumento@unipi.it}
}
\seealso{
\code{\link{ctqr}}
}
\examples{
# Using simulated data

n <- 1000
x1 <- runif(n)
x2 <- runif(n)
t <- 1 + x1 + x2 + runif(n, -1,1)
c <- rnorm(n,3,1)
y <- pmin(t,c)
d <- (t <= c)

model <- ctqr(Surv(y,d) ~ x1 + x2, p = c(0.25,0.5))
pred <- predict(model) # the same as fitted(model)
predict(model, newdata = data.frame(x1 = c(0.2,0.6), x2 = c(0.1,0.9)), se.fit = TRUE)

}
\keyword{regression}

