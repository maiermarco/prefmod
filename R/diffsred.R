diffsred <- function(x, nobj){
##############################################################
# calculate differences from non-paired comparison responses #
# (i.e. ratings/Likert or rankings)                          #
# Sinclair ordering: (1-2)(1-3)(2-3)(1-4)                    #
# reduces diffs to -1/0/1                                    #
##############################################################
  
  ncomp <- nobj * (nobj - 1L) / 2
  diffs <- matrix(c(0:0), nrow(x), ncomp)
  col <- 1
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      diffs[, col] <- ifelse(x[, i] > x[, j], 1, ifelse(x[, i] < x[, j], -1, 0))
      col <- col + 1
    }
  }
  
  diffs

}
