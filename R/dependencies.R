dependencies <- function(nobj, Y){
  
  m <- matrix(, nrow = nobj - 1L, ncol = nobj)
  
  label.comps <- NULL
  for(j in 2L:nobj){
    for(i in 1L:(j - 1L)){
      label.comps <- c(label.comps, paste(i, j, sep = ""))
    }
  }
  
  k <- 1
  for(j in 2L:nobj){
    for(i in 1L:(j - 1L)){
      m[j - 1L, i] <- k     # lower triangle
      m[i, j]      <- k     # upper triangle
      k            <- k + 1
    }
  }
  
  label.intpars <- NULL
  dncol <- (nobj) * (nobj - 1L) * (nobj - 2L) / 2
  d <- matrix(, ncol = dncol, nrow = nrow(Y))
  l <- 1
  for(k in 1L:(nobj)){
    for(i in 1L:(nobj - 2L)){
      for(j in (i + 1L):(nobj - 1L)){
        d[, l] <- Y[, m[i, k]] * Y[, m[j, k]]
        label.intpars <- c(label.intpars, paste("I", label.comps[m[i, k]], ".", label.comps[m[j, k]], sep = ""))
        l <- l + 1
      }
    }
  }
  
  list(d = d, label.intpars = label.intpars)

}
