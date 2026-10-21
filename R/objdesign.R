objdesign <- function(nrows, nobj, nrespcat){
# design matrix for objects #
#############################
  
  obj <- matrix(c(0:0), nrows, nobj)
  row <- 1
  
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      d <- nrespcat
      for(c in seq.int(0L, nrespcat)){
        obj[row + c, i] <-  d
        obj[row + c, j] <- -d
        d <- d - 2
      }
      row <- row + nrespcat + 1
    }
  }
  
  obj

}
