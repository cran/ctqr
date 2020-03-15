#' @importFrom stats sd prcomp model.matrix model.response delete.response model.frame terms lm.wfit lm.fit
#' @importFrom stats model.weights model.offset printCoefmat coef nobs vcov .getXlevels pnorm
#' @importFrom survival Surv is.Surv
#' @importFrom graphics plot points abline text
#' @importFrom utils menu getFromNamespace
#' @import pch


# automatically finds a "good" estimator of a CDF
findagoodestimator <- function(dat,w){
  
  splx <- pch::splinex() 
  CDF <- suppressWarnings(pch::pchreg(
    Surv(z,y,d) ~ ., data = dat, weights = w, splinex = splx)) 
  fit.ok <- (CDF$conv.status == 0)

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
# Note: safe.pch and quickpred are only really "needed" in gal.
quickpred <- function(obj,y){
	end.y <- obj$u(y)
	n <- length(y)
	t <- y - obj$breaks[end.y]
	ind <- cbind(1:n,end.y)
	lambda <- obj$lambda[ind]
	Lambda <- obj$Lambda[ind] + lambda*t
	cbind(lambda = lambda, Lambda = Lambda)	
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


# Estimation algorithm
qr.gs <- function(beta0, check, p, CDF,
a = 0.5, b = 1.5, maxit = 1000, tol = 1e-6, type){

	if(type == "trunc"){ee <- ee.cens.trunc}
	else if(type == "cens"){ee <- ee.cens}
	else{ee <- ee.u}
	Hy <- (if(type != "u") quickpred(CDF, check$y0)[,2] else 0)
	Hz <- (if(type == "trunc") quickpred(CDF, check$z0)[,2] else 0)
	CDF <- list(Sy = pmax(exp(-Hy), 1e-2), Sz = pmax(exp(-Hz), 1e-2))

	Beta <- n.it <- converged <- NULL
	fitted <- omega.z <- omega.y <- si <- NULL
	for(u in 1:length(p)){

		tau <- p[u]
		beta <- beta0[,u]
		S <- ee(beta, tau, CDF, check)
		s <- S$s; obj <- sum(s^2)
		delta <- 0.25/max(abs(s) + 0.001)

		####

		for(i in 1:maxit){
			beta.step <- delta*s; beta.step[is.na(beta.step)] <- Inf
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
	}

	# de-scaling fitted (but not beta)

	ay <- attributes(check$y)
	fitted <- fitted*(ay$M - ay$m)/10 + ay$m
	list(beta = Beta, n.it = n.it, converged = converged, CDF = CDF,
		fitted = fitted, omega.z = omega.z, omega.y = omega.y, si = si,
		logLik = rep(NA, length(p)) # for compatibility with gal
	)
}


# main fitting function. There will be a future argument "method", "qr" or "gal".
#' @export
ctqr <- function(formula, data, weights, p = 0.5, CDF, control = ctqr.control(), ...){

	cl <- match.call()

	if((CDF.in <- !missing(CDF))){
		if(class(CDF) != "pch")
		{stop("'CDF' must be an object of class 'pch'", call. = FALSE)}
		if(any(is.na(match(all.vars(formula(cl)), all.vars(CDF$call)))))
			{stop("some of the variables in 'formula' is not included in 'CDF'")}
	}
	
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "weights", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	if(CDF.in){mf$subset <- rownames(CDF$mf)}
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	if(CDF.in){attr(mf, "na.action") <- 
		unique(c(attr(mf, "na.action"), attr(CDF$mf, "na.action")))}

	mt <- attr(mf, "terms")
	if(!is.Surv(zyd <- model.response(mf)))
		{stop("the model response must be created with Surv()")}
	if((n <- nrow(zyd)) == 0){stop("zero non-NA cases")}
	type <- attributes(zyd)$type
	if(type == "right"){
		y <- zyd[,1]; z <- rep.int(-Inf,n); d <- zyd[,2]
		type <- (if(any(d == 0)) "cens" else "u")
	}
	else if(type == "counting"){
		z <- zyd[,1]; y <- zyd[,2]; d <- zyd[,3]; type <- "trunc"
		if(!any(z > min(y))){type <- (if(any(d == 0)) "cens" else "u")}
	}

	x <- model.matrix(mt,mf)
	weights <- model.weights(mf)
	off <- model.offset(mf)
 
	# plug-in ###################################################################

	if(missing(CDF)){
		dat <- data.frame(z = z, y = y, d = d, x = x)
		CDF <- findagoodestimator(dat,weights)
	}

	# fit #######################################################################

	method <- "qr"; fitfun <- qr.gs # or "gal", gal.gs
	check <- check.in.ctqr(z,y,d,x,off,weights, CDF$breaks)
	beta0 <- start(CDF, p, check$x, check$y, check$off, check$weights)
	CDF <- safe.pch(CDF, min(1e-6, 1/n/10))
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
	asy <- asy.qr # or asy.gal
	if(any(fit.ok)){
		A <- asy(z,y,d,xfullrank, weights, p, CDF, fit, fit.ok = which(fit.ok))
		covar[which(fit.ok)] <- A
	}
	for(j in 1:h){colnames(covar[[j]]) <- rownames(covar[[j]]) <- colnames(xfullrank)}

	# Finish ####################################################################

	Beta <- matrix(NA,ncol(x),h)
	rownames(Beta) <- colnames(x)
	Beta[ax$sel,] <- beta
	fit$beta <- Beta

	npar <- rank + (method == "gal") + sum(CDF$beta != 0)
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
	colnames(fit$beta) <-  names(fit$logLik) <- names(covar) <- names(fit$n.it) <- 
		names(fit$converged) <- colnames(fit$fitted) <- paste("p =",p)

	fit <- list(p = p, coefficients = fit$beta, call = cl, 
		n.it = fit$n.it, converged = fit$converged,
		fitted = fit$fitted, terms = mt, mf = mf, 
		covar = covar
	)
	if(method == "gal"){fit$logLik <- l}

	if(any(fit$n.it[fit.ok] == control$maxit)){warning("convergence has not been achieved")}
	if(any(failed)){warning("estimation failed at some quantiles")}
	class(fit) <- "ctqr"
	fit
}

#' @export
print.ctqr <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	"\n\n", sep = "")

	cat("Coefficients:\n")
	print.default(format(x$coef, digits = digits), print.gap = 2L, quote = FALSE)

	cat("\n")
	cat("Degrees of freedom:", nrow(x$mf), "total;", nrow(x$mf) - attr(x$mf, "npar"), "residuals")

	cat("\n")
	invisible(x)
}


#' @export
summary.ctqr <- function(object, ...){
	
	p <- object$p
	covar <- object$covar
	b <- object$coef
	non_na <- attr(object$mf, "xin")
	q <- length(non_na)

	out <- list()
	for(j in 1:length(p)){
		coe <- cbind(b[,j],NA,NA,NA)
		colnames(coe) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
		rownames(coe) <- rownames(b)
		coe[non_na,2] <- sqrt(diag(covar[[j]]))
		coe[,3] <- coe[,1]/coe[,2]
		coe[,4] <- 2*(1 - pnorm(abs(coe[,3])))
		out[[j]]<- coe
	}

	names(out) <- names(covar)
	out <- list(call = object$call, coefficients = out)
	class(out) <- "summary.ctqr"
	out
}

#' @export
print.summary.ctqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

	cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	coe <- x$coefficients
	cat("Coefficients:\n\n")
	for(j in 1:length(coe)){
		cat(paste(names(coe)[j], "\n"))
		printCoefmat(coe[[j]], digits = digits, signif.stars = TRUE, cs.ind = 1:2, tst.ind = 3, 
			P.values = TRUE, has.Pvalue = TRUE)
		cat("\n")
	}
	invisible(x)
}



#' @export
coef.ctqr <- function(object, ...){
	out <- object$coef
	if(length(object$p) == 1){out <- out[,1]}
	out
}
#' @export
nobs.ctqr <- function(object, ...){nrow(object$model)}
#' @export
vcov.ctqr <- function(object, ...){
	out <- object$covar
	if(length(object$p) == 1){out <- out[[1]]}
	out
}



#' @export
predict.ctqr <- function(object, newdata, se.fit = FALSE, ...){
	
	p <- object$p
	beta <- object$coef
	beta[is.na(beta)] <- 0
	mf <- object$mf
	mt <- terms(mf)
	miss <- attr(mf, "na.action")
	xin <- attr(mf, "xin")
	fit.ok <- attr(mf, "fit.ok")
	nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:(nrow(mf) + length(miss)))[-miss])
	xlev <- .getXlevels(mt, mf)

	if(!missing(newdata)){
		mt <- delete.response(mt)
		if(any(is.na(match(all.vars(mt), colnames(newdata)))))
			{stop("'newdata' must contain all x-variables")}
		mf <- model.frame(mt, data = newdata, xlev = xlev)

		if(nrow(mf) == 0){
			nr <- nrow(newdata)
			out <- data.frame(matrix(NA,nr,length(p)))
			colnames(out) <- paste("p",p, sep = "")
			return(out)
		}
		miss <- attr(mf, "na.action")
		nomiss <- (if(is.null(miss)) 1:nrow(mf) else (1:nrow(newdata))[-miss])
	}

	x <- model.matrix(mt, mf)
	if(is.null(off <- model.offset(mf))){off <- 0}
	n <- length(c(miss, nomiss))
	eta <- matrix(NA, n, length(p))
	eta[nomiss,] <- x%*%beta + off

	colnames(eta) <- paste("p", p, sep = "")
	rownames(eta)[nomiss] <- rownames(mf)
	if(!is.null(miss)){rownames(eta)[miss] <- names(miss)}
	eta[,!fit.ok] <- NA
	out <- as.data.frame(eta)

	if(se.fit){
		x <- x[,xin, drop = FALSE]
		se <- out
		for(j in 1:length(p)){
			se[nomiss,j] <- sqrt(diag(x%*%object$covar[[j]]%*%t(x)))
		}
		out <- list(fit = out, se.fit = se)
	}
	out
}


#' @export
plot.ctqr <- function(x, which = NULL, ask = TRUE, ...){
  
  
  if(sum(fit.ok <- attr(x$mf, "fit.ok")) <= 1)
  {stop("at least 2 quantiles must be fitted to use plot.ctqr")}
  p <- x$p
  beta <- x$coefficients
  se <- sapply(x$covar, function(x) sqrt(diag(x)))
  o <- order(p)
  p <- p[o]
  beta <- beta[,o]
  se <- se[,o]
  
  plot.ctqr.int <- function(p, beta,se,j,L){
    beta <- beta[j,]
    if(all(is.na(beta))){
      plot(0,0, col = "white", xlim = c(0,1), ylim = c(-1,1),
           xlab = "p", ylab = "beta(p)", main = L$labels[j])
      text(0.5,0, "no coefficients to show")
    }
    else{
      se <- se[j,]
      low <- beta - 1.96*se
      up <- beta + 1.96*se
      if(is.null(L$ylim)){L$ylim <- range(c(low,up), na.rm = TRUE)}
      
      plot(p, beta, xlab = L$xlab, ylab = L$ylab, main = L$labels[j], 
           type = "l", lwd = L$lwd, xlim = L$xlim, ylim = L$ylim, col = L$col)
      
      points(p, low, lty = 2, lwd = L$lwd, type = "l", col = L$col)
      points(p, up, lty = 2, lwd = L$lwd, type = "l", col = L$col)
      abline(h = 0, lty = 3)
    }
  }
  
  L <- list(...)
  if(is.null(L$xlim)){L$xlim = c(0.01,0.99)}
  if(is.null(L$lwd)){L$lwd <- 2}
  if(is.null(L$col)){L$col <- "black"}
  if(is.null(L$xlab)){L$xlab <- "p"}
  if(is.null(L$ylab)){L$ylab <- "beta(p)"}
  L$labels <- rownames(x$coefficients)
  q <- length(L$labels)
  
  
  if(!is.null(which) | !ask){
    if(is.null(which)){which <- 1:q}
    for(j in which){plot.ctqr.int(p,beta,se,j,L)}
  }
  else{
    pick <- 1
    while(pick > 0 && pick <= q){
      pick <- menu(L$labels, title = "Make a plot selection (or 0 to exit):\n")
      if(pick > 0 && pick <= q){plot.ctqr.int(p,beta,se,pick,L)}
    }
  }
}








#' @export
ctqr.control <-
function(tol = 1e-6, maxit = 1000, a = 0.5, b = 1.25){

	if(tol <= 0){warning(
		"the value of tol supplied is zero or negative; the default value of 1e-6 was used instead",
		call. = FALSE)
		tol <- 1e-06
	}

	if((maxit <- round(maxit)) <= 0){warning(
		"the value of maxit supplied is not valid; the default value of 1000 was used instead",
		call. = FALSE)
		maxit <- 1000
	}

	if(a <= 0 | a >= 1){warning(
		"the value of a supplied is not between 0 and 1; the default value of 0.5 was used instead",
		call. = FALSE)
		a <- 0.5
	}

	if(b <= 1){warning(
		"the value of b supplied is not > 1; the default value of 1.25 was used instead",
		call. = FALSE)
		b <- 1.25
	}

	list(tol = tol, maxit = maxit, a = a, b = b)
}



check.in.ctqr <-
  function(z,y,d,x,off,weights, breaks){
    
    # x ################################################################

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
    
    r <- range(y)
    m <- (if(int) r[1] else 0)
    M <- r[2]

    z <- pmax(z, breaks[1] + (breaks[2] - breaks[1])/length(y))
    z0 <- z; y0 <- y
    y <- c(y - m)/(M - m)*10
    z <- c(z - m)/(M - m)*10
    attributes(y) <- list(m = m, M = M)
    
    # offset and weights #################################################
  
    if(!is.null(off)){off <- off*10/(M - m)}
    else{off <- 0}
    if(!is.null(weights)){weights <- weights/sum(weights)*length(y)}
    else{weights <- 1}

    # finish #############################################################
    
    list(z0 = z0, z = z, y0 = y0, y = y, d = c(1 - d), x = x, off = c(off), weights = c(weights))
  }



firstep <- function(obj, z,y,d,w){
  
  n <- length(y)
  q <- ncol(x <- obj$x)
  beta <- obj$beta
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
  # Not needed in practice, since I already have CDF$covar

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
  V <- obj$covar[sel, sel, drop = FALSE]
  beta <- c(beta)[sel]

  # output: beta without 0/-Inf, hessian, score.i, and derivatives of Hy and Hz
  # all quantities have size length(beta).

  list(beta = beta, H = H, V = V, Si = Si, dHy = DHy, dHz = DHz)
}






asy.qr <- function(z,y,d,x, weights, p, CDF, fit, fit.ok){
  
	est.H2 <- function(fitted, y,d,x,weights, a){
		eps <- 1e-7
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
		u.l <- d*(1 - pnorm(y - fitted.l,0,eps)) # d*oy.l
		u.r <- d*(1 - pnorm(y - fitted.r,0,eps)) # d*oy.r
		u <- (u.r - u.l)/(2*eps)
		-t(x)%*%(x*c(weights*u))
	}


	n <- length(y)
	eta <- fit$fitted
	oy <- fit$omega.y
	oz <- fit$omega.z
	Sy <- pmax(fit$CDF$Sy, 1e-10)
	Sz <- pmax(fit$CDF$Sz, 1e-10)
	if(is.null(w1 <- model.weights(CDF$mf))){w1 <- 1}
	w2 <- (if(!is.null(weights)) weights else 1)
	si <- fit$si; si[is.na(si)] <- 0
  
  
	# Notation: 
	# S = "first derivatives" for outer product
	# H = "second derivative"/jacobian/hessian
	# V = inverse of H (in particular, V1 is a covariance matrix)
	# There is a single S1 and H1, but a different S2, H2, V2, H12 for each quantile
  
	#### First step

	A <- firstep(CDF, z,y,d,w1)
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
			H2 <- est.H2(eta[,j], y,d,x,w2, a = a)
			try.V2 <- try(chol2inv(chol(-H2)), silent = TRUE)
			if(!inherits(try.V2, "try-error")){V2[[j]] <- try.V2; V2.ok <- TRUE}
			else{a <- a + 0.01}
		}
	}

	#### Finish

	COV <- list()
	for(j in fit.ok){
	  S <- cbind(S2[[j]], S1)
		Sw <- cbind(S2[[j]]*w2, S1*w1)
		Omega <- rbind(cbind(V2[[j]], V2[[j]]%*%H12[[j]]%*%V1), cbind(t(H12[[j]])*0, V1))		
		Omega <- Omega%*%(t(S)%*%Sw)%*%t(Omega)
		COV[[length(COV) + 1]] <- Omega[1:ncol(x), 1:ncol(x), drop = FALSE]
	}
	names(COV) <- p[fit.ok]
	COV
}







