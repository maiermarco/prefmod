# aus e1071
countpattern2 <- function(x){
  
  intpatt <- apply(as.matrix(x), 1L, bin2int)
  intpatt <- intpatt + 1
  pat     <- tabulate(intpatt, nbins = 2L^ncol(as.matrix(x)))
  
  pat

}
