expand <- function(x, k){
  xx <- x
  for(i in seq.int(2L, k)) xx <- rbind(x, xx)
  xx
}



gqz <- function(numnodes = 20L, minweight = 0.000001){
# Calculate Gaussian Quadrature points for the Normal distribution
# using the abscissas and weights for Hermite integration. The
# conversion of the locations and weights is given in Lindsey (1992,
# page 169:3) and Skrondal & Rabe-Hesketh (2004, page 165:1).
# The argument numnodes is the theoretical number of quadrature points,
# locations with weights that are less than the argument minweight will
# be omitted.
# The default vale of minweight = 0.000001 returns 14 masspoints for the
# default numnodes = 20 as in Aitkin, Francis & Hinde (2005).
  out <- gauss.quad(numnodes, "hermite")
  h   <- rbind(out$nodes * sqrt(2.0), out$weights / sum(out$weights))
# Sort the locations and weights into columns in decending order of the
# location vector.
  ord <- order(h[1L, ], decreasing = TRUE)
  h <- h[, ord]
  h <- cbind(h[1L, ], h[2L, ])
  h <- subset(as.data.frame(h), (h[, 2L] >= minweight))
  names(h) <- c("location", "weight")
  h
}



gauss.quad <- function(n, kind = "legendre", alpha = 0.0, beta = 0.0){
  
  n <- as.integer(n)
  if(n < 0L) stop("need non-negative number of nodes")
  if(n == 0L) return(list(nodes = numeric(0L), weights = numeric(0L)))
  
  kind <- match.arg(kind, c("legendre", "chebyshev1", "chebyshev2", "hermite", "jacobi", "laguerre"))
  i    <- seq_len(n)
  i1   <- i[-n]
  
  switch(kind,
    legendre = {
      muzero <- 2.0
      a      <- rep(0.0, n)
      b      <- i1 / sqrt(4.0 * i1^2.0 - 1.0)
    },
    chebyshev1 = {
      muzero <- pi
      a      <- rep(0.0, n)
      b      <- rep(0.5, n - 1.0)
      b[1L]  <- sqrt(0.5)
    },
    chebyshev2 = {
      muzero <- pi / 2.0
      a      <- rep(0.0, n)
      b      <- rep(0.5, n - 1.0)
    },
    hermite = {
      muzero <- sqrt(pi)
      a      <- rep(0.0, n)
      b      <- sqrt(i1 / 2.0)
    },
    jacobi = {
      ab     <- alpha + beta
      muzero <- 2.0^(ab + 1.0) * gamma(alpha + 1.0) * gamma(beta + 1.0) / gamma(ab + 2.0)
      a      <- i
      a[1L]  <- (beta - alpha) / (ab + 2.0)
      i2     <- seq.int(2L, n)
      abi    <- ab + 2.0 * i2
      a[i2]  <- (beta^2.0 - alpha^2.0) / (abi - 2.0) / abi
      b      <- i1
      b[1L]  <- sqrt(4.0 * (alpha + 1.0) * (beta + 1.0) / (ab + 2.0)^2.0 / (ab + 3.0))
      i2     <- i1[-1L]
      abi    <- ab + 2.0 * i2
      b[i2]  <- sqrt(4.0 * i2 * (i2 + alpha) * (i2 + beta) * (i2 + ab) / (abi^2.0 - 1.0) / abi^2.0)
    },
    laguerre = {
      a      <- 2.0 * i - 1.0 + alpha
      b      <- sqrt(i1 * (i1 + alpha))
      muzero <- gamma(alpha + 1.0)
    }
  )
  
  A <- rep(0.0, n * n)
  A[(n + 1L) * (i - 1L) + 1L] <- a
  A[(n + 1L) * (i1 - 1L) + 2L] <- b
  A[(n + 1L) * i1] <- b
  dim(A) <- c(n, n)
  vd <- eigen(A, symmetric = TRUE)
  w <- rev(as.vector(vd$vectors[1L, ]))
  w <- muzero * w^2.0
  x <- rev(vd$values)
  
  list(nodes = x, weights = w)

}



weightslogl.calc.w <- function(p, fjk, weights){
# amended to rescale posterior probabilities
# R.E.D. 15th Feb 2006
# p is a vector of length K containing the mixture proportions
# fjk is a JXK matrix of log densities
  logpf  <- t(apply(fjk, 1L, "+", log(p)))
  Mi     <- apply(logpf, 1L, max)
  Sik    <- logpf - Mi
  ifelse(Sik < -760.0, 0.0, Sik)
  eSik   <- exp(Sik)
  SeSik  <- as.vector(apply(eSik, 1L, sum))
  w      <- eSik / SeSik
  ML.dev <- -2.0 * sum(weights * log(apply(eSik, 1L, sum))) - 2.0 * sum(weights * Mi)
  
  list(w = w, ML.dev = ML.dev)

}



summary.pattNPML <- function(object, digits = max(3L, getOption("digits") - 3L), ...){
  
  np <- length(object$coefficients)
  if(np > 0L){
    m <- seq_len(np)[substr(attr(object$coefficients, "names"), 1L, 4L) == "MASS"]
    mass.points <- object$coefficients[m]
    cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
    cat("Coefficients")
    cat(":\n")
    df.r <- object$lastglm$df.residual

    glm.dispersion <- if(any(object$family$family == c("poisson", "binomial"))) 1.0 else
      if(df.r > 0){
        sum(object$lastglm$weights * object$lastglm$residuals^2.0, na.rm = TRUE) / df.r
      } else Inf
    lastglmsumm <- summary.glm(object$lastglm, dispersion = glm.dispersion)
    fitcoef     <- matrix(lastglmsumm$coeff[, seq_len(3L)], ncol = 3L, dimnames = list(dimnames(lastglmsumm$coef)[[1L]], c(dimnames(lastglmsumm$coeff)[[2L]][seq_len(2L)], "t value") ))  #03-08-06
    print(fitcoef)
  } else {
    cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
    cat("No coefficients. \n")
  }
  
  p <- object$masses
  #names(p) <- paste("MASS", seq_len(ncol(object$post.prob)), sep = "") # omitted from 0.42 on
  dispersion <- 1.0   # now calculate dispersion in the sense of "dispersion parameter":
  
  cat("\nMixture proportions:\n")
  print.default(format(p, digits), print.gap = 2L, quote = FALSE)
  
  if(object$family$family == "gaussian"){
    dispersion <- (object$sdev$sdev)^2.0
    if(object$Misc$lambda <= 1.0 / length(object$masses)){
      cat("\nComponent distribution - MLE of sigma:\t  ", format(signif(object$sdev$sdev, digits)), "\n")
    } else {
      cat("\nMLE of component standard deviations:\n")
      s <- object$sdev$sdevk
      names(s) <- names(p)
      print.default(format(s, digits), print.gap = 2L, quote = FALSE)
      cat("\n")
    }
  } else if(object$family$family == "Gamma"){
    dispersion <- 1 / object$shape$shape
    if(object$Misc$lambda <= 1 / length(object$masses)){
      cat("\nComponent distribution - MLE of shape parameter:\t  ", format(signif(object$shape$shape, digits)), "\n")
    } else {
      cat("\n MLE of component shape parameters:\n")
      s <- object$shape$shapek
      names(s) <- names(p)
      print.default(format(s, digits), print.gap = 2L, quote = FALSE)
      cat("\n")
    }
  } else cat("\n")

  ###  cat("Random effect distribution - standard deviation:\t  ", format(object$rsdev), "\n\n") # 03/09/07
  ##cat("\n")

  cat("Deviance:\t   ", format(round(object$deviance, digits = 2L)), "\n")
  cat("-2 log L:\t   ", format(round(object$disparity, digits = 2L)))
  conv <- object$EMconverged
  if(conv != "none"){
    cnv <- ifelse(conv, "", "no")
    if(!is.null(object$post.prob)) cat("   ", cnv, "EM convergence at iteration ", round(object$EMiter, 0L))
  }
  cat("\n")

  if(np > 0L){
    invisible(c("coefficients" = list(fitcoef), object[-1L], list(dispersion = dispersion), list(lastglmsumm = lastglmsumm)))
  } else {
    invisible(c("coefficients" = list(fitcoef), object[-1L], list(dispersion = dispersion)) )
  }
}



print.pattNPML <- function(x, digits = max(3L, getOption("digits") - 3L), ...){
  np <- length(x$coefficients)
  # print(names(x))
  if(np > 0L){
    m <- seq_len(np)[substr(attr(x$coefficients, "names"), 1L, 4L) == "MASS"]
    mass.points <- x$coefficients[m]
    cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
    cat("Coefficients")
    cat(":\n")
    print.default(format(x$coefficients[seq_len(np)], digits = digits), print.gap = 2L, quote = FALSE) #;cat("\n")
  } else {
    cat("\nCall: ", deparse(x$call), "\n", fill = TRUE)
    cat("No coefficients. \n")
  }

  if(x$family$family == "gaussian" && x$Misc$lambda <= 1.0 / length(x$masses)){ # print sigma only if it is constant over components
    cat("Component distribution - MLE of sigma:\t  ", format(signif(x$sdev$sdev, digits)), "\n")   # 03/09/07
  }
  
  if(x$family$family == "Gamma" && x$Misc$lambda <= 1.0 / length(x$masses)){ # print shape only if it is constant over components
    cat("Component distribution - MLE of shape parameter:\t  ", format(signif(x$shape$shape, digits)), "\n")# 03/09/07
  }

  ###  cat("Random effect distribution - standard deviation:\t  ", format(x$rsdev), "\n\n") # 03/09/07
  cat("\n")

  if(!is.null(x$post.prob)){
    p <- x$masses
    # names(p) <- paste("MASS", seq_len(ncol(x$post.prob)), sep = "") # omitted from 0.42
    cat("Mixture proportions")
    cat(":\n")
    print.default(format(p, digits), print.gap = 2L, quote = FALSE)
    cat("\n")
  }
  cat("Deviance:\t   ", format(round(x$deviance, digits = 2L)), "\n")
  cat("-2 log L:\t   ", format(round(x$disparity, digits = 2L)))
  conv <- x$EMconverged
  if(conv != "none"){
    cnv <- ifelse(conv, "", "no")
    if(!is.null(x$post.prob)) cat("   ", cnv, "EM convergence at iteration ", round(x$EMiter, 0L))
  }
  cat("\n")
  invisible(x)
}



# obsolete
#"BIC.pattNPML" <-
#function(object, ...){
#object$disparity + (length(c(na.omit(object$coefficients), object$masses))-1) * log(sum(object$data$y))
#}
# Achim's proposal - thanks Achim
# allows for using stats AIC(), BIC()
logLik.pattNPML <- function(object, ...){
  structure(-0.5 * object$disparity,
    df = length(c(na.omit(object$coefficients), object$masses)) - 1L,
    class = "logLik")
}



nobs.pattNPML <- function(object, ...){
  sum(object$data$y)
}
