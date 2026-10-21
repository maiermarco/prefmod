differences <- function(in_mat){
##############################################################
# calculate differences from non-paired comparison responses #
# (i.e. ratings/Likert or rankings)                          #
# Sinclair ordering: (1-2)(1-3)(2-3)(1-4)                    #
##############################################################
  
  ncomp <- get("ncomp", get("ENV", environment(patt.design)))
  nobj  <- get("nobj", get("ENV", environment(patt.design)))
  
  diffs <- matrix(c(0:0), nrow(in_mat), ncomp)
  col <- 1
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      diffs[, col] <- in_mat[, i] - in_mat[, j]
      col <- col + 1
    }
  }
  
  diffs

}
