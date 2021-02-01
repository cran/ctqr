#' @importFrom stats sd prcomp model.matrix model.response delete.response model.frame terms lm.wfit lm.fit
#' @importFrom stats model.weights model.offset printCoefmat coef nobs vcov .getXlevels pnorm
#' @importFrom survival Surv is.Surv
#' @importFrom graphics plot points abline text
#' @importFrom utils menu getFromNamespace
#' @import pch



check.in.ctqr <-
  function(z,y,d,x,off,weights, breaks){
    
    # x ################################################################
    
    n <- nrow(x)
    v <- qr(x)
    sel <- v$pivot[1:v$rank]
    x <- x[,sel, drop = FALSE]
    q <- ncol(x)
    MX <- apply(x,2,mean)
    SdX <- apply(x,2,sd)
    int <- any(const <- which(SdX == 0))
    vars <- (if(int) (1:q)[-const] else (1:q))
    
    if(int){
      SdX[const] <- x[1,const]
      MX[const] <- 0
    }
    else{MX <- rep.int(0,q)}
    x <- scale(x, center = MX, scale = SdX)
    
    attX <- list(sel = sel, int = int, const = const, 
                 vars = vars, centerX = MX, scaleX = SdX)
    if(length(vars) > 1){
      PC <- prcomp(x[,vars,drop = FALSE], center = FALSE, scale = FALSE)
      newX <- scale(PC$x, center = int)
      x[,vars] <- newX
      attX$rot <- PC$rotation
      attX$centerPC <- (if(int) attr(newX, "scaled:center") else 0)
      attX$scalePC <- attr(newX, "scaled:scale")
    }
    
    attributes(x) <- c(attributes(x), attX)
    
    # scaling y and z ###################################################
    
    r <- range(y[is.finite(y)])
    m <- (if(int) r[1] else 0)
    M <- r[2]
    
    z <- pmax(z, breaks[1] + (breaks[2] - breaks[1])/n)
    z0 <- z; y0 <- y
    y <- (y - m)/(M - m)*10
    z <- (z - m)/(M - m)*10
    attr(y, "m") <- m; attr(y, "M") <- M
    
    # offset and weights #################################################
    
    if(!is.null(off)){off <- off*10/(M - m)}
    else{off <- rep(0,n)}
    if(!is.null(weights)){weights <- weights/mean(weights)}
    else{weights <- 1}
    
    # finish #############################################################
    
    list(z0 = z0, z = z, y0 = y0, y = y, d = c(1 - d), x = x, off = c(off), weights = c(weights))
}



# Compute starting points
start <- function(CDF, p, x, y, off, weights){

  predQ.pch <- getFromNamespace("predQ.pch", ns = "pch")
	lambda <- CDF$lambda
	predQ <- predQ.pch(CDF, p = p)
	m <- attr(y, "m"); M <- attr(y, "M")
	if(length(weights) == 1){weights <- rep.int(1, nrow(x))}
	if(length(off) == 1){off <- rep.int(0, nrow(x))}
	out <- NULL
	for(j in 1:length(p)){
		qp <- (predQ[,j] - m)/(M - m)*10
		ww <- weights
		ww[qp > 10] <- 0
		fit <- lm.wfit(x,qp, ww, offset = off)
		coe <- fit$coef
		out <- cbind(out, coe)
	}
	out[is.na(out)] <- 0
	out
}

# The estimating equations

ee.u <- function(beta, tau, CDF, V){

	eta <- c(V$x%*%cbind(beta) + V$off)
	oy <- (V$y <= eta)
	si <- tau - oy
	s <- colSums(V$x*(si*V$weights))
	list(s = s, eta = eta, omega.y = oy, omega.z = 1, si = si)
}

ee.cens <- function(beta, tau, CDF, V){
	eta <- c(V$x%*%cbind(beta) + V$off)
	oy <- (V$y <= eta)
	si <- tau - oy - V$d*oy*(tau - 1)/CDF$Sy
	s <- colSums(V$x*(si*V$weights), na.rm = TRUE) # some Sy could be 0
	list(s = s, eta = eta, omega.y = oy, omega.z = 1, si = si)
}

ee.cens.trunc <- function(beta, tau, CDF, V){
	eta <- c(V$x%*%cbind(beta) + V$off)
	oy <- (V$y <= eta)
	oz <- (V$z <= eta)
	si <- -oy - V$d*oy*(tau - 1)/CDF$Sy + oz + oz*(tau - 1)/CDF$Sz
	s <- colSums(V$x*(si*V$weights), na.rm = TRUE) # some Sy or Sz could be 0
	if(!any(oy)){s[1] <- Inf} # all omega.y = 0 IS almost a solution.
	list(s = s, eta = eta, omega.y = oy, omega.z = oz, si = si)
}

ee.icens <- function(beta, tau, CDF, V){

  m <- attr(V$y, "m"); M <- attr(V$y, "M")
  eta <- c(V$x%*%cbind(beta) + V$off)
  o1 <- (V$y[,2] <= eta) # know omega = 1
  o0 <- (V$y[,1] >= eta) # know omega = 0
  ou <- which(ohat <- !(o1 | o0)) # don't know omega
  Heta <- quickpred(CDF$CDF, eta*((M - m)/10) + m)[,2]
  Seta <- exp(-Heta)
  
  omega.y <- o1
  omega.y[ou] <- ((CDF$S1 - Seta)/CDF$deltaS)[ou]
  attr(omega.y, "ohat") <- ohat # = TRUE if omega is estimated

  si <- tau - omega.y
  s <- colSums(V$x*(si*V$weights))
  list(s = s, eta = eta, omega.y = omega.y, omega.z = 1, si = si)
}



# Estimation algorithm
qr.gs <- function(beta0, check, p, CDF,
a = 0.5, b = 1.5, maxit = 1000, tol = 1e-6, type){

	if(type == "trunc"){ee <- ee.cens.trunc}
	else if(type == "cens"){ee <- ee.cens}
  else if(type == "interval"){ee <- ee.icens}
	else{ee <- ee.u}
  
  if(type != "interval"){
	  Hy <- (if(type != "u") quickpred(CDF, check$y0)[,2] else 0)
	  Hz <- (if(type == "trunc") quickpred(CDF, check$z0)[,2] else 0)
	  CDF <- list(Sy = pmax(exp(-Hy), 1e-2), Sz = pmax(exp(-Hz), 1e-2), CDF = CDF)
  }
  else{
    H1 <- quickpred(CDF, check$y0[,1])[,2]
    H2 <- quickpred(CDF, check$y0[,2])[,2]
    S1 <- exp(-H1); S1[check$y0[,1] == -Inf] <- 1
    S2 <- exp(-H2); S2[check$y0[,2] == Inf] <- 0
    CDF <- list(S1 = S1, S2 = S2, deltaS = pmax(S1 - S2, 1e-6), CDF = CDF)
  }

	Beta <- n.it <- converged <- NULL
	fitted <- omega.z <- omega.y <- si <- ohat <- NULL
	for(u in 1:length(p)){

		tau <- p[u]
		beta <- beta0[,u]
		S <- ee(beta, tau, CDF, check)
		s <- S$s; obj <- sum(s^2)
		delta <- 0.25/max(abs(s) + 0.001)

		####

		for(i in 1:maxit){

			beta.step <- delta*s; beta.step[is.na(beta.step)] <- 0
			if(max(abs(beta.step)) < tol){break}
			new.beta <- beta + beta.step
			new.S <- ee(new.beta, tau, CDF, check)
			new.s <- new.S$s
			new.obj <- sum(new.s^2)
			if(new.obj < obj){
				beta <- new.beta
				obj <- new.obj
				S <- new.S
				s <- new.s
				delta <- delta*b
			}
			else{delta <- delta*a}
		}    
    
		Beta <- cbind(Beta, beta)
		n.it[u] <- i
		converged[u] <- (i < maxit)
		
		fitted <- cbind(fitted, S$eta)
		omega.y <- cbind(omega.y, S$omega.y)
		omega.z <- cbind(omega.z, S$omega.z)
		si <- cbind(si, S$si)
		ohat <- cbind(ohat, attr(S$omega.y, "ohat"))
	}

	# de-scaling fitted (but not beta)

	ay <- attributes(check$y)
	fitted <- fitted*(ay$M - ay$m)/10 + ay$m
	list(beta = Beta, n.it = n.it, converged = converged, CDF = CDF,
		fitted = fitted, omega.z = omega.z, omega.y = omega.y, si = si, ohat = ohat,
		logLik = rep(NA, length(p)) # for compatibility with gal
	)
}


# main fitting function. There will be a future argument "method", "qr" or "gal".
#' @export
ctqr <- function(formula, data, weights, p = 0.5, CDF, control = ctqr.control(), ...){

	cl <- match.call()

	if((CDF.in <- !missing(CDF))){
		if(!inherits(CDF, "pch"))
		  {stop("'CDF' must be an object of class 'pch'", call. = FALSE)}
		if(any(is.na(match(all.vars(formula), attr(CDF$call, "all.vars")))))
		  {stop("some of the variables in 'formula' is not included in 'CDF'")}
	}

	mf <- match.call(expand.dots = FALSE)
	mf$formula <- formula
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	if(CDF.in){mf$subset <- rownames(CDF$mf)}
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
  if(CDF.in){attr(mf, "na.action") <- attr(CDF$mf, "na.action")}
	mt <- attr(mf, "terms")
	
	# Handling null weights
	
	if(any((w <- model.weights(mf)) < 0)){stop("negative 'weights'")}
	if(is.null(w)){w <- rep.int(1, nrow(mf)); alarm <- FALSE}
	else{
	  alarm <- (w == 0)
	  sel <- which(!alarm)
	  mf <- mf[sel,]
	}
	if(any(alarm)){warning("observations with null weight will be dropped", call. = FALSE)}
	if((n <- nrow(mf)) == 0){stop("zero non-NA cases", call. = FALSE)}
	
	####
	
	convert.Surv <- getFromNamespace("convert.Surv", ns = "pch")
	if(!is.Surv(zyd <- model.response(mf)))
		{stop("the model response must be created with Surv()")}
	if((n <- nrow(zyd)) == 0){stop("zero non-NA cases")}
	type <- attributes(zyd)$type

	zyd <- cbind(zyd)
	if(type == "right"){
		y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]
		type <- (if(any(d == 0)) "cens" else "u")
	}
	else if(type == "counting"){
		z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]; type <- "trunc"
		if(!any(z > min(y))){type <- (if(any(d == 0)) "cens" else "u")}
	}
	else if(type == "interval"){y <- convert.Surv(zyd); z <- rep.int(-Inf,n); d <- zyd[,3]}
	else{stop("only 'right', 'counting', and 'interval2' data are supported")}
	if(!(any(d %in% c(1,3)))){stop("all observation are censored")}
	
	x <- model.matrix(mt,mf)
	weights <- model.weights(mf)
	off <- model.offset(mf)

	# plug-in ###################################################################

	if(missing(CDF)){
		if(type == "interval"){dat <- data.frame(y1 = y[,1], y2 = y[,2], x = x)}
	  else{dat <- data.frame(z = z, y = y, d = d, x = x)}
		CDF <- findagoodestimator(dat, weights, type)
	}

	# fit #######################################################################

	if(any(p <= 0 | p >= 1)){stop("0 < p < 1 is required")}
	p <- sort(p)
	method <- "qr"; fitfun <- qr.gs # or "gal", gal.gs
	check <- check.in.ctqr(z,y,d,x,off,weights, CDF$breaks)
	beta0 <- start(CDF, p, check$x, check$y, check$off, check$weights)
	CDF0 <- CDF; CDF <- safe.pch(CDF, min(1e-6, 1/n/10))
	fit <- fitfun(beta0, check, p, CDF, a = control$a, b = control$b, 
		maxit = control$maxit, tol = control$tol, type = type)

	# de-scaling beta ###########################################################

	ay <- attributes(check$y)
	ax <- attributes(check$x)

	beta <- fit$beta
	if(length(ax$vars) > 1){
		beta[ax$vars,] <- beta[ax$vars,]/ax$scalePC
		beta[ax$const,] <- beta[ax$const,] - colSums(ax$centerPC*beta[ax$vars,, drop = FALSE])
		beta[ax$vars,] <- ax$rot%*%cbind(beta[ax$vars,])
	}
	beta[ax$vars,] <- beta[ax$vars,]/ax$scaleX[ax$vars]
	beta[ax$const,] <- beta[ax$const,] - colSums(ax$centerX*beta)
	beta <- beta*(ay$M - ay$m)/10
	beta[ax$const,] <- beta[ax$const,] + ay$m
	beta[ax$const,] <- beta[ax$const,]/ax$scaleX[ax$const]

	# asymptotic ################################################################

	h <- length(p)
	fit.ok <- colMeans(fit$omega.y)
	fit.ok <- (fit.ok > pmin(0.001,p) & fit.ok < pmax(0.999,p))
	xfullrank <- x[, ax$sel, drop = FALSE]
	rank <- ncol(xfullrank)

	covar <- list()
	for(j in 1:h){covar[[j]] <- matrix(, rank, rank)}
	asy <- (if(type != "interval") asy.qr.ct else asy.qr.ic) # or asy.gal
	if(any(fit.ok)){
		A <- asy(z,y,d,xfullrank, check$weights, p, fit, fit.ok = which(fit.ok))
		covar[which(fit.ok)] <- A
	}
	for(j in 1:h){colnames(covar[[j]]) <- rownames(covar[[j]]) <- colnames(xfullrank)}

	# Finish ####################################################################

	Beta <- matrix(NA,ncol(x),h)
	rownames(Beta) <- colnames(x)
	Beta[ax$sel,] <- beta
	fit$beta <- Beta

	npar <- c(q1 = sum(CDF$beta != 0), q2 = rank + (method == "gal"))
	attr(fit$logLik, "df") <- npar
	l <- fit$logLik
	
	if(any(failed <- !fit.ok)){
		fit$beta[,failed] <- NA
		fit$logLik[failed] <- NA
		fit$converged[failed] <- NA
		fit$n.it[failed] <- NA
		fit$fitted[,failed] <- NA		
	}

	
	attr(mf, "xin") <- ax$sel
	attr(mf, "fit.ok") <- fit.ok
	attr(mf, "npar") <- npar
	attr(mf, "contrasts") <- attr(x, "contrasts")
	colnames(fit$beta) <-  names(fit$logLik) <- names(covar) <- names(fit$n.it) <- 
		names(fit$converged) <- colnames(fit$fitted) <- paste("p =", p)

	fit <- list(p = p, coefficients = fit$beta, call = cl, 
		n.it = fit$n.it, converged = fit$converged,
		fitted = fit$fitted, terms = mt, mf = mf, 
		covar = covar, CDF = CDF0
	)
	if(method == "gal"){fit$logLik <- l}

	if(any(fit$n.it[fit.ok] == control$maxit)){warning("convergence has not been achieved")}
	if(any(failed)){warning("estimation failed at some quantiles")}
	class(fit) <- "ctqr"
	fit
}


# GENERAL COMMENT ON ASYMPTOTICS
# Solve S1(theta) = 0, S2(beta, theta) = 0.

# The formula is the following:
# \hat_cov(\hat_\beta) = D^-1 * Omega * D^-1
# where Omega is the outer product of (S2, S1),
# and D = rbind(c(dS2/dbeta, dS2/dtheta), c(0,dS1/dtheta))
# To get V = D^-1, I use block inversion.

# for censored and truncated data
asy.qr.ct <- function(z,y,d,x, weights, p, fit, fit.ok){
  
	est.H2 <- function(fitted, y,d,x,weights, a, CDF){
		eps <- 1e-6
		ok <- FALSE
		while(!ok){
			eps <- eps*2
			fitted.l <- fitted - eps
			fitted.r <- fitted + eps
			oy.l <- (y <= fitted.l)
			oy.r <- (y <= fitted.r)
			check <- mean(oy.l != oy.r)
			ok <- (check > a)
		}
		u.l <- d*(1 - pnorm(y - fitted.l,0,eps)) # a smooth version of d*oy.l
		u.r <- d*(1 - pnorm(y - fitted.r,0,eps)) # a smooth version of d*oy.r
		u <- (u.r - u.l)/(2*eps)
		-t(x)%*%(x*c(weights*u))
	}

  CDF <- fit$CDF
	n <- length(y)
	eta <- fit$fitted
	oy <- fit$omega.y
	oz <- fit$omega.z
	Sy <- pmax(CDF$Sy, 1e-10)
	Sz <- pmax(CDF$Sz, 1e-10)
	if(is.null(w1 <- model.weights(CDF$CDF$mf))){w1 <- 1}
	w2 <- (if(!is.null(weights)) weights else 1)
	w1 <- w1/mean(w1); w2 <- w2/mean(w2)
	si <- fit$si; si[is.na(si)] <- 0
  
  
	# Notation: 
	# S = "first derivatives" for outer product
	# H = "second derivative"/jacobian/hessian
	# V = inverse of H (in particular, V1 is an estimate of the covariance matrix of the first step).
	# There is a single S1 and H1, but a different S2, H2, V2, H12 for each quantile
  
	#### First step

	A <- firstep.ct(CDF$CDF, z,y,d,w1)
	V1 <- A$V
	S1 <- A$Si

	#### Second step, and mixed

	S2 <- V2 <- H12 <- list()
	for(j in fit.ok){
		tau <- p[j]
		S2[[j]] <- x*si[,j]
		H12[[j]] <- t(x*w2)%*%((1 - d)*oy[,j]*(1 - tau)*A$dHy/Sy - oz[,j]*(1 - tau)*A$dHz/Sz)
		V2.ok <- FALSE
		a <- 0.05
		while(!V2.ok){
			H2 <- est.H2(eta[,j], y,d,x,w2, a = a, CDF = CDF)
			try.V2 <- try(chol2inv(chol(-H2)), silent = TRUE)
			if(!inherits(try.V2, "try-error")){V2[[j]] <- try.V2; V2.ok <- TRUE}
			else{a <- a + 0.01}
		}
	}

	#### Finish

	COV <- list()
	q <- ncol(x); ind2 <- 1:q
	for(j in fit.ok){
	  S <- cbind(S2[[j]], S1)
		Sw <- cbind(S2[[j]]*w2, S1*w1)
		U <- t(S)%*%Sw

		Omega <- rbind(cbind(V2[[j]], V2[[j]]%*%H12[[j]]%*%V1), cbind(t(H12[[j]])*0, V1))		
		Omega <- Omega%*%(U)%*%t(Omega)
		Omega <- Omega[1:q, 1:q, drop = FALSE]
		
		# If Omega is too much larger than the naif (which happens in a very small proportion of cases)
		# I assume that "something went wrong" and just return the naif.
		
		Omega_naif <- V2[[j]]%*%U[ind2,ind2]%*%t(V2[[j]]) # Naif covariance matrix
		test <- sqrt(diag(Omega)/diag(Omega_naif))
		if(any(test > 5)){Omega <- Omega_naif} 

		# no special motivation for the multiplicative factor. Empirically, test is typically less than 2.
	
		###################
		
		COV[[length(COV) + 1]] <- Omega
	}
	names(COV) <- p[fit.ok]
	COV
}


# for interval-censored data
# Note that E[\hat\omega] = F(x*beta), implying that the first derivative
# of S = x*(p - \omega) w.r.t. beta is -x'x*f(x*beta). I could use the first-step
# estimator to compute this quantity. However, taking numerical derivatived of \hat\omega
# appears more stable and even more accurate.

asy.qr.ic <- function(z,y,d,x, weights, p, fit, fit.ok){
  
  omega <- function(y,eta,CDF){

    o1 <- (y[,2] <= eta) # know omega = 1
    o0 <- (y[,1] >= eta) # know omega = 0
    ou <- which(ohat <- !(o1 | o0)) # don't know omega
    Heta <- quickpred(CDF$CDF, eta)[,2]
    Seta <- exp(-Heta)
    omega.y <- o1
    omega.y[ou] <- ((CDF$S1 - Seta)/CDF$deltaS)[ou]
    omega.y
  }
  
  est.H2 <- function(fitted, y,d,x,weights, a, CDF){
    eps <- 1e-6
    ok <- FALSE
    while(!ok){
      eps <- eps*2
      fitted.l <- fitted - eps
      fitted.r <- fitted + eps
      oy.l <- omega(y,fitted.l,CDF)
      oy.r <- omega(y,fitted.r,CDF)
      check <- mean(oy.l != oy.r)
      ok <- (check > a)
    }

    u <- (oy.r - oy.l)/(2*eps)
    -t(x)%*%(x*c(weights*u))
  }
  
  CDF <- fit$CDF
  if(is.null(w1 <- model.weights(CDF$CDF$mf))){w1 <- 1}
  w2 <- (if(!is.null(weights)) weights else 1)
  w1 <- w1/mean(w1); w2 <- w2/mean(w2)
  si <- fit$si; si[is.na(si)] <- 0

  # Notation: 
  # S = "first derivatives" for outer product
  # H = "second derivative"/jacobian/hessian
  # V = inverse of H (in particular, V1 is an estimate of the covariance matrix of the first step)
  # There is a single S1 and H1, but a different S2, H2, V2, H12 for each quantile
  
  #### First step #######################################################
  # Note: obj$covar is a sandwitch, and is not equal to the inverse of H1.
  
  H1 <- attr(CDF$CDF$mf, "h") # hessian
  S1 <- attr(CDF$CDF$mf, "s.i") # score.i
  V1 <- safesolve(H1) # inverse of hessian, a covariance matrix
  # note: first step is ML, but I actually minimize, that's why H1 and not -H1

  # CDF and its derivatives
  
  F.y1 <- 1 - CDF$S1
  F.y2 <- 1 - CDF$S2
  deltaF <- pmax(CDF$deltaS, 1e-4)
  deltaF_square <- deltaF^2
  dF.y1 <- (1 - F.y1)*firstep.ic(CDF$CDF, y[,1])
  dF.y2 <- (1 - F.y2)*firstep.ic(CDF$CDF, y[,2])
  deltadF <- dF.y2 - dF.y1
  deltamix <- F.y1*dF.y2 - F.y2*dF.y1

  #### Second step, and mixed ###########################################

  eta <- fit$fitted
  oy <- fit$omega.y
  ohat <- fit$ohat
  n <- nrow(y)
  S2 <- V2 <- H12 <- list()
  for(j in fit.ok){
    tau <- p[j]

    # H12
    pred.eta <- quickpred(CDF$CDF, eta[,j])
    S.eta <- exp(-pred.eta[,2])
    F.eta <- 1 - S.eta
    dF.eta <- S.eta*firstep.ic(CDF$CDF, eta[,j])
    h12 <- (dF.eta*deltaF - F.eta*deltadF + deltamix)/deltaF_square*ohat[,j]
    H12[[j]] <- t(x*w2)%*%(h12)

    # S2, V2
    S2[[j]] <- x*si[,j]
    V2.ok <- FALSE
    a <- 0.05 + mean(ohat[,j])
    while(!V2.ok){
      H2 <- est.H2(eta[,j], y,d,x,w2, a = a, CDF = CDF)
      try.V2 <- try(chol2inv(chol(-H2)), silent = FALSE)
      if(!inherits(try.V2, "try-error")){V2[[j]] <- try.V2; V2.ok <- TRUE}
      else{a <- a + 0.02}
    }
    
    # Another way of computing V2, that was dismissed
    # f.eta <- S.eta*pred.eta[,1]
    # H2 <- t(x*w2)%*%(x*f.eta)
    # V2[[j]] <- safesolve(H2)
  }

  #### Finish
  
  COV <- list()
  q <- ncol(x); ind2 <- 1:q
  for(j in fit.ok){
    
    # Outer product
    S <- cbind(S2[[j]], S1)
    Sw <- cbind(S2[[j]]*w2, S1*w1)
    U <- t(S*ohat[,j])%*%Sw # note the ohat!!!
    U[ind2,ind2] <- t(S[,ind2])%*%Sw[,ind2] # no ohat here
    U[-ind2,-ind2] <- t(S[,-ind2])%*%Sw[,-ind2] # no ohat here

    # Sandwitch
    Omega <- rbind(cbind(V2[[j]], -V2[[j]]%*%H12[[j]]%*%V1), cbind(t(H12[[j]])*0, V1))		
    Omega <- Omega%*%(U)%*%t(Omega)
    Omega <- Omega[1:q, 1:q, drop = FALSE]
    
    # If Omega is too much larger than the naif (which happens in a very small proportion of cases)
      # I assume that "something went wrong" and just return the naif.
    
    Omega_naif <- V2[[j]]%*%U[ind2,ind2]%*%t(V2[[j]]) # Naif covariance matrix
    test <- sqrt(diag(Omega)/diag(Omega_naif))
    if(any(test > 5)){Omega <- Omega_naif} 
    # no special motivation for the multiplicative factor. Empirically, test is typically less than 2.

    ###################
    
    COV[[length(COV) + 1]] <- Omega
  }
  names(COV) <- p[fit.ok]
  COV
}





