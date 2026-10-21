# PC design matrix    (12), (13), (23), (14), (24), ...
pcdesign <- function(nobj){
  
  ncomp <- nobj * (nobj - 1L) / 2          # number of comparisons
  obj   <- matrix(0, ncomp, nobj)
  row <- 1
  
  for(j in 2L:nobj){
    for(i in 1L:(j - 1L)){
      obj[row, i] <-  1
      obj[row, j] <- -1
      row <- row + 1
    }
  }
  
  obj

}
