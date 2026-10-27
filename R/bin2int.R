# bin vector of 0 / 1
# int integer
bin2int <- function(bin){
  
  if(length(bin) > 31L) stop("\"bin\" must not contain more than 31 elements (overflow).")
  # TODO packBits --> double (64)
  
  storage.mode(bin) <- "integer"
  
  n   <- length(bin)
  int <- 0L
  
  for(k in seq.int(0L, n - 1L)) int <- int + bin[n - k] * as.integer(2L^k)
  
  int

}
