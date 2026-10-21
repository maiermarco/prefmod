# print method for pattMod objects
print.pattMod <- function(x, ...){                     #cmat, nobj, elim, ilabels, ENV)
  
  envList <- x$envList
  obj     <- x$result
  
  #chi2 <- ((cnts - sum(cnts) * p.patt)^2) / (sum(cnts) * p.patt)
  #chi2 <- sum(chi2)
  #lm <- sum(log(p.patt) * cnts)
  #lf <- sum(log(p.cnts[p.cnts > 0]) * cnts[cnts > 0])
  #deviance <- -2 * (lm - lf)
  
  type <- switch(envList$resptype,
    "paircomp" = "paired comparison",
    "rating"   = "ratings",
    "ratingT"  = "repeated ratings",
    "ranking"  = "rankings"
  )
  
  cat("\nResults of pattern model for", type, "\n")
  cat("\nCall:\n")
  wcut <- getOption("width")
  if(wcut > 60L) wcut <- 60L
  str <- deparse(x$call, width.cutoff = wcut)
  cat(str, "\n", fill = TRUE)
  cat("Deviance: ", 2 * abs(envList$ll - envList$fl), "\n")
  cat("log likelihood: ", envList$ll, "\n")
  if(envList$elim != "~1") cat("eliminated term(s): ", paste(envList$elim, sep = "", collapse = ""), "\n")
  cat("\nno of iterations: ", obj$iterations, " (Code:", obj$code, ")\n")

  cov.names <- c("", colnames(envList$covdesmat)[-1L])
  if(regexpr("T", envList$resptype) > 0){                # for time models
    idxvec <- -envList$nitems * seq_len(envList$tpoints)     # we need a different definition of nobj
  } else {
    idxvec <- seq_len(envList$nobj - 1L)
  }
  obj.names <- envList$obj.names[idxvec]
  
  all.obj.names <- rep(obj.names, length(cov.names))
  #obj.names <- rep(LETTERS[seq_len(envList$nobj - 1L)], length(cov.names))
  
  #cov.names <- rep(cov.names, rep(envList$nobj - 1L, length(cov.names)))
  cov.names <- rep(cov.names, rep(length(obj.names), length(cov.names))) # covnames repeated for group of all lambdas
  
#browser()
  #colon <- c(rep("", envList$nobj - 1L), rep(":", length(obj.names) - envList$nobj + 1L))
  colon <- c(rep("", length(obj.names)), rep(":", length(cov.names) - length(obj.names)))
  par.names <- paste(obj.names, colon, cov.names, sep = "")
  
  se   <- round(sqrt(diag(solve(obj$hessian))), digits = 5L)
  zval <- round(obj$estimate / se, digits = 3L)
  pval <- round(2.0 * (1.0 - pnorm(abs(zval))), digits = 4L)
  
  coef.table <- cbind(round(obj$estimate, digits = 5L), se, zval, pval)
  
  if(envList$NItest) par.names <- c(par.names, paste("mis.", par.names, sep = ""))
  
  if(envList$resptype == "paircomp"){ # currently only for PC models
  
    if(envList$MISalpha || envList$MISbeta){
      lenallest  <- length(obj$estimate)
      lenmisspar <- lenallest - envList$nobj - length(envList$ilabels) + 1L
      if(envList$MISalpha) par.names <- c(par.names, paste("mis.alpha", seq_len(envList$nobj)[envList$Malph], sep = ""))
      if(envList$MISbeta)  par.names <- c(par.names, paste("mis.beta",  seq_len(envList$nobj)[envList$Mbeta], sep = ""))
    }

    if(envList$MIScommon) par.names <- c(par.names, "mis.alpha.common")
  
  }
  
  if(!is.null(envList$tpoints)){       # rep measurements
    if(envList$undec) {
      U.names   <- paste("UT", seq_len(envList$tpoints), sep = "")
      par.names <- c(par.names, U.names)
    }
  } else {                              # only one time point
    if(envList$undec) {
      par.names <- c(par.names, "U")
    }
  }
  
  if(!is.null(envList$ilabels))  par.names <- c(par.names, envList$ilabels)
  if(!is.null(envList$iTlabels)) par.names <- c(par.names, envList$iTlabels)
  rownames(coef.table) <- par.names
  colnames(coef.table) <- c("estimate", "se", "z", "p-value")
  cat("\n")
  print(coef.table)

  invisible(x)

}
