# tern vector of 0/1/2
# int integer
tern2int <- function(tern)
{
  n   <- length(tern)
  dec <- 0
  for(k in 0L:(n - 1L))
    dec <- dec + tern[n - k] * 3^k
  dec

}
