# matrix of all patterns
patternmat2 <- function(nobj){   # pattStr <- names(countpattern(dat[1L, ]))
  
  nvar <- nobj * (nobj - 1L) / 2
  Y    <- matrix(0, 2L^nvar, nvar)
  for(i in seq_len(nvar)) Y[, nvar + 1L - i] <- rep(rep(c(0, 1), c(2L^(i - 1L), 2L^(i - 1L))), 2L^(nvar - i))
  
  Y

}
