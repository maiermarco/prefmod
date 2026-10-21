# simulates paired comparison data with given worth
# 2 response categories, no undecided category
#
simPC <- function(nobj, nobs, worth = NULL, seed = NULL, pr = FALSE){
  
  if(is.null(worth)){
    worth <- runif(nobj)
  } else if(length(worth) != nobj){
    stop("length of worth is not equal to nobj\n")
  }
  worth <- worth / sum(worth)
  if(pr) cat("used worth parameters are: ", worth, "\n")
  
  probs <- NULL
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      probs <- c(probs, worth[i] / (worth[i] + worth[j]))
    }
  }
  
  ncomp <- choose(nobj, 2L)
  if(!is.null(seed)) set.seed(seed)
  
  data <- NULL
  for(i in seq_len(ncomp)) data <- cbind(data, sample(c(1, -1), nobs, replace = TRUE, prob = c(probs[i], 1.0 - probs[i])))
  
  data.frame(data)

}

