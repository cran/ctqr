

# automatically finds a "good" estimator of a CDF.
# IMPORTANT REMARK: the CDF can either be computed internally, or be supplied by the user.
# This is why the CDF is ALWAYS computed on the unscaled response and covariates.

# With right-censored and left-truncated data, this is quite irrelevant: the value of F(y)
# and F(z) are computed in advance and are scale-invariant. Instead, with interval-censored
# data, I have to compute F(x*beta) at a current beta: in ee.icens, this requires to de-scale
# the value of x*beta.

findagoodestimator <- function(dat,w, type){
  
  if(type == "interval"){
    CDF <- suppressWarnings(pch::pchreg(
      Surv(y1,y2, type = "interval2") ~ ., data = dat, weights = w, splinex = NULL))
  }
  else{
    CDF <- suppressWarnings(pch::pchreg(
      Surv(z,y,d) ~ ., data = dat, weights = w, splinex = NULL))
  }
  
  fit.ok <- (CDF$conv.status == 0)

  splx <- pch::splinex()
  br <- length(CDF$breaks)
  delta.br <- 1
  count <- 0
  while(!fit.ok){
    
    if(count == 32){warning("could not find a CDF with conv.status = 0"); break}
    br <- br + delta.br
  
    CDF <- suppressWarnings(pch::pchreg(
    Surv(z,y,d) ~ ., data = dat, weights = w, splinex = splx, breaks = br))
    count <- count + 1
  
    if(count == 5){br <- br - 5; delta.br <- -1}
    if(count == 10){br <- br + 5; delta.br <- 0; splx$df <- 3; splx$v <- 0.95}
    if(count == 11){delta.br <- 1}
    if(count == 16){br <- br - 5; delta.br <- -1}
    if(count == 21){br <- br + 5; delta.br <- 0; splx$df <- 1; splx$v <- 1}
    if(count == 22){delta.br <- 1}
    if(count == 27){br <- br - 5; delta.br <- -1}
    fit.ok <- (CDF$conv.status == 0)
  }
  CDF
}




# Prepare a safe version of the CDF to ensure that F(t | x) > 0.
# The output is ready for "quickpred".
safe.pch <- function(obj, eps = 1e-6){

	br <- obj$breaks
	c0 <- 2*br[1] - br[length(br)] # new lower break

	H0 <- -log(1 - eps) # cum.haz at br[1], such that F(br[1]) = eps
	h0 <- H0/(br[1] - c0) # hazard between c0 and br[1]

	br <- c(c0, br); k <- length(br) - 1
	attr(br, "k") <- k; attr(br, "h") <- br[2:(k + 1)] - br[1:k]

	list(lambda = cbind(h0, obj$lambda), Lambda = cbind(0,H0,obj$Lambda + H0), breaks = br, 
		u = attr(obj$mf, "u"), x = obj$x, mf = obj$mf,
		beta = obj$beta, covar = obj$covar)
}

# quick predictor for modified pch objects
# Note: safe.pch and quickpred are only really "needed" in gal (not implemented yet) 
  # and with interval-censored data, where at each iteration I need to evaluate F(xbeta).
quickpred <- function(obj,y){
	end.y <- obj$u(y)
	n <- length(y)
	t <- y - obj$breaks[end.y]
	ind <- cbind(1:n,end.y)
	lambda <- obj$lambda[ind]
	Lambda <- obj$Lambda[ind] + lambda*t
	cbind(lambda = lambda, Lambda = Lambda)	
}

# Computes all relevant first-step quantities to be used for asymptotics
# (to be used when the data are right-censored or left-truncated)
firstep.ct <- function(obj, z,y,d,w){
  
  n <- length(y)
  q <- ncol(x <- obj$x)
  lambda <- obj$lambda
  br <- obj$breaks
  k <- attr(br, "k")
  h <- attr(br, "h")
  u <- obj$u

  # handling obs before the smallest break
  
  out.l <- which(z < br[1])
  z[out.l] <- br[1]
  
  ## relevant quantities
  
  L1 <- t(matrix(br[1:k], k,n))
  L2 <- t(matrix(h, k,n))  
  ty <- pmax(y - L1,0)
  tz <- pmax(z - L1,0)
  ty <- pmin(ty, L2)
  tz <- pmin(tz, L2)
  DD <- matrix(0,n,k)
  DD[cbind(1:n, u(y))] <- d

  ## derivative of Hy and Hz (each column to be multiplied by x)

  dHy <- lambda*ty
  dHz <- lambda*tz

  ## hessian (to be used in t(x)%*%(x*hi[,j]))

  hi <- dHz - dHy

  ## score (each column to be multiplied by x)
  
  si <- DD + hi

  ###### final quantities
  # in the safe.pch object there was an additional break with no associated parameters
  # so k is 1 less, the first column of lambda is not estimated, etc.

  Si <- DHy <- DHz <- NULL
  H <- matrix(0, q*(k - 1), q*(k - 1))

  for(j in 2:k){
    jj <- j - 1
    ind <- (jj*q - q + 1):(jj*q)
    H[ind,ind] <- t(x)%*%(x*(w*hi[,j]))
    Si <- cbind(Si, x*si[,j])
    DHy <- cbind(DHy, x*dHy[,j])
    DHz <- cbind(DHz, x*dHz[,j])
  }
  sel <- which(apply(H,2, function(a) any(a != 0)))
  H <- H[sel,sel, drop = FALSE]
  Si <- Si[,sel, drop = FALSE]
  DHy <- DHy[,sel, drop = FALSE]
  DHz <- DHz[,sel, drop = FALSE]
  V <- safesolve(-H) # Note: obj$covar is a sandwitch, and is not equal to the inverse of H.
  
  # Output: inverse of hessian, hessian, score.i, and derivatives of Hy and Hz.

  list(V = V, H = H, Si = Si, dHy = DHy, dHz = DHz)
}


# Computes the derivative of \hat H(y) with respect to theta.
# (to be used when the data are interval-censored).
# Remark: the score and hessian are already provided as attributed of CDF$mf,
# and they already exclude the parameters that could not be estimated.

firstep.ic <- function(obj, y){
  
  n <- length(y)
  q <- ncol(x <- obj$x[,apply(obj$beta,1, function(a) any(a != 0)), drop = FALSE])
  lambda <- obj$lambda
  br <- obj$breaks
  k <- attr(br, "k")
  h <- attr(br, "h")
  
  ## relevant quantities
  
  L1 <- t(matrix(br[1:k], k,n))
  L2 <- t(matrix(h, k,n))
  
  t.y <- pmax(y - L1,0)
  t.y <- pmin(t.y, L2)
  
  ## derivative of H(y) (each column to be multiplied by x)
  
  dH.y <- lambda*t.y
  
  ###### final quantities
  # in the safe.pch object there was an additional break with no associated parameters
  # so k is 1 less, the first column of lambda is not estimated, etc.
  
  DH.y <- NULL
  
  for(j in 2:k){
    jj <- j - 1
    ind <- (jj*q - q + 1):(jj*q)
    DH.y <- cbind(DH.y, x*dH.y[,j])
  }
  DH.y
}




