
#' @export
print.ctqr <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	"\n\n", sep = "")

	cat("Coefficients:\n")
	print.default(format(x$coef, digits = digits), print.gap = 2L, quote = FALSE)

	cat("\n")
	cat("Number of observations:", nrow(x$mf), "\n")

	N <- attr(x$CDF$mf, "n.events")
	
	if(inherits(x$CDF, "ct")){
	  cat("N. of events: ", paste(deparse(round(N)), sep = " ", collapse = " "), "\n", sep = "")
	}
	else{
	  cat("Non-censored: ", paste(deparse(round(N[[1]])), sep = " ", collapse = " "), "\n", sep = "")
	  cat("Left-censored: ", paste(deparse(round(N[[2]])), sep = " ", collapse = " "), "\n", sep = "")
	  cat("Right-censored: ", paste(deparse(round(N[[3]])), sep = " ", collapse = " "), "\n", sep = "")
	  cat("Interval-censored: ", paste(deparse(round(N[[4]])), sep = " ", collapse = " "), "\n", sep = "")
	}
	cat("Number of free parameters:", attr(x$mf, "npar")[1], "(CDF) +", attr(x$mf, "npar")[2], "(ctqr)")

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
	
	out <- list(call = object$call, coefficients = out, 
	   nobs = nrow(object$CDF$mf), n.events = attr(object$CDF$mf, "n.events"), 
	   npar = attr(object$mf, "npar"), class = class(object$CDF)[2])
	class(out) <- "summary.ctqr"
	out
}

#' @export
print.summary.ctqr <- function(x, digits = max(3L, getOption("digits") - 3L), ...){

	cat("\nCall: ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  ###################
  
  cat("\n")
  cat("Number of observations:", x$nobs, "\n")
  
  N <- x$n.events
  
  if(x$class == "ct"){
    cat("N. of events: ", paste(deparse(round(N)), sep = " ", collapse = " "), "\n", sep = "")
  }
  else{
    cat("Non-censored: ", paste(deparse(round(N[[1]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("Left-censored: ", paste(deparse(round(N[[2]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("Right-censored: ", paste(deparse(round(N[[3]])), sep = " ", collapse = " "), "\n", sep = "")
    cat("Interval-censored: ", paste(deparse(round(N[[4]])), sep = " ", collapse = " "), "\n", sep = "")
  }
  cat("Number of free parameters:", x$npar[1], "(CDF) +", x$npar[2], "(ctqr)")
  
  ###################
  
  cat("\n\n")
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
nobs.ctqr <- function(object, ...){nrow(object$mf)}
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
	contr <- attr(mf, "contrasts")
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

	x <- model.matrix(mt, mf, contrasts.arg = contr)
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
  
  plot.ctqr.int <- function(p, beta,se,j,L){
    beta <- beta[j,]
    if(all(is.na(beta))){
      plot(0,0, col = "white", xlim = c(0,1), ylim = c(-1,1),
           xlab = "p", ylab = "beta(p)", main = L$labels[j])
      text(0.5,0, "no coefficients to show")
    }
    else{
      se <- rbind(se)[j,]
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



safesolve <- function(A){
  
  eps <- 1e-6
  ok <- FALSE
  count <- 0
  while(!ok){
    B <- try(chol2inv(chol(A + diag(eps - 1e-6, nrow(A)))), silent = TRUE)
    ok <- (!inherits(B, "try-error"))
    eps <- eps*5
    count <- count + 1
  }
  warn <- (count > 1)
  attr(B, "warn") <- warn
  B
}


# This function is actually not used.
# Its use is to receive a covariance matrix, and modify the off-diagonal elements
 # to ensure that cor(A) <= c.
regularize_corr <- function(A, c = 0.9){
  
  q <- ncol(A)
  if(q == 1){return(A)}
  s <- sign(A)
  A <- abs(A)
  se <- sqrt(diag(A))
  for(i1 in 1:(q - 1)){
    for(i2 in (i1 + 1):q){
      cc <- s[i1,i2]*min(A[i1,i2], c*se[i1]*se[i2])
      A[i1,i2] <- A[i2,i1] <- cc
    }
  }
  A
}



