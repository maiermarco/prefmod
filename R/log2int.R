log2int <- function(bin)
{
  n   <- length(bin)
  int <- 0
  for(k in 0L:(n - 1L))
    int <- int + as.numeric(bin[n - k]) * 2^k
  int
}
