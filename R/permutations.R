permutations <- function(n){
################################################################
# permutations from package e1071
################################################################
  
  if(n == 1L) return(matrix(1)) else if(n < 2L) stop("n must be a positive integer")
  
  z <- matrix(1)
  
  for(i in seq.int(2L, n)){
    x <- cbind(z, i)
    a <- c(seq_len(i), seq_len(i - 1L))
    z <- matrix(0, ncol = ncol(x), nrow = i * nrow(x))
    z[seq_len(nrow(x)), ] <- x
    for(j in seq.int(2L, i) - 1L){
      z[j * nrow(x) + seq_len(nrow(x)), ] <- x[, a[seq_len(i) + j]]
    }
  }
  dimnames(z) <- NULL
  
  z

}
