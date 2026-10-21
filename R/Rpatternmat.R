# matrix of all rank patterns
Rpatternmat <- function(nobj){
  
  Rpatt <- perm(nobj, nobj) # permutations from gtools
  # RpattStr <- apply(Rpatt, 1, function(x) paste(x, sep = "", collapse = ""))
  # ENV$RpattStr <- sort(RpattStr)

  ncomp <- nobj * (nobj - 1L) / 2
  diffs <- matrix(, ncol = ncomp, nrow = nrow(Rpatt))
  c     <- 0
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      c <- c + 1
      diffs[, c] <- Rpatt[, i] - Rpatt[, j]
    }
  }
  
  diffs

}
