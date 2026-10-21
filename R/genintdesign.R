genintdesign <- function(pattern, sinclair = TRUE){
################################################################
# interactions design matrix
################################################################
#
# generates design matrix for interaction effects between two comparisons
# with one common object
#

  nobj    <- get("nobj", get("ENV", environment(patt.design)))
  ncomp   <- nobj * (nobj - 1L) / 2      # number of comparisons
  compvec <- matrix(c(0:0), ncomp, 2) # matrix to describe which objects in which comparison
  row     <- 1
  
  if(sinclair){                 # (12), (13), (23), (14), (24), ...
    for(j in seq.int(2L, nobj)) {
      for(i in seq_len(j - 1L)){
        compvec[row, 1L] <- i
        compvec[row, 2L] <- j
        row <- row + 1
      }
    }
  } else {                        # (12), (13), (14), ..., (23), (24), ...
    for(i in seq_len(nobj - 1L)) {
      for(j in (i + 1L):nobj){
        compvec[row, 1L] <- i
        compvec[row, 2L] <- j
        row <- row + 1
      }
    }
  }
  
  vec1 <- NULL
  vec2 <- NULL
  vec3 <- NULL
  desc.names <- NULL
  for(i in seq_len(ncomp - 1L)){
    for(j in (i + 1L):ncomp){
      if(compvec[i, 1L] == compvec[j, 1L]){
        vec1 <- c(vec1, i)
        vec2 <- c(vec2, j)
        vec3 <- c(vec3, +1)  # sign + because common object preferred in both comps
        desc.names <- c(desc.names, paste("I", compvec[i, 1L], ".", compvec[i, 2L], compvec[j, 2L], sep = "", collapse = ""))
      }
      if(compvec[i, 2L] == compvec[j, 2L]){
        vec1 <- c(vec1, i)
        vec2 <- c(vec2, j)
        vec3 <- c(vec3, +1)  # sign + because common object not preferred in both comps
        desc.names <- c(desc.names, paste("I", compvec[i, 2L], ".", compvec[i, 1L], compvec[j, 1L], sep = "", collapse = ""))
      }
      if(compvec[i, 1L] == compvec[j, 2L]){
        vec1 <- c(vec1, i)
        vec2 <- c(vec2, j)
        vec3 <- c(vec3, -1)  # sign - because common object preferred/not preferred
        desc.names <- c(desc.names, paste("I", compvec[i, 1L], ".", compvec[i, 2L], compvec[j, 1L], sep = "", collapse = ""))
      }
      if(compvec[i, 2L] == compvec[j, 1L]){
        vec1 <- c(vec1, i)
        vec2 <- c(vec2, j)
        vec3 <- c(vec3, -1)  # sign - because common object preferred/not preferred
        desc.names <- c(desc.names, paste("I", compvec[i, 2L], ".", compvec[i, 1L], compvec[j, 2L], sep = "", collapse = ""))
      }
    }
  }
  multvec <- cbind(vec1, vec2)   # each row contains column numbers of comparison patterns
                                 # to be multiplied for interactioneffect


  nintpars <- length(desc.names) # number of interaction parameters
  assign("nintpars", nintpars, get("ENV", environment(patt.design)))

  # pattern <- pattern - median(range(pattern))
  # pattern <- pattern * 2
  # print(pattern)

  designmatrix <- NULL
  for(i in seq_len(nintpars)){
    work         <- pattern[, multvec[i, 1L]] * pattern[, multvec[i, 2L]] * vec3[i]
    designmatrix <- cbind(designmatrix, work)
  }
  
  colnames(designmatrix) <- desc.names    # names for interaction parameters,
                                          #  e.g., I3.45, comparisons 3-4 and 3-5
  designmatrix

}
