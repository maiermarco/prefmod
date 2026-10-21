# x integer
# r vector of 0/1/2
int2tern <- function(x){
  
  r <- vector(mode = "integer", length = 0L)
  
  if(x <= 2L) r <- x else{
    while(x > 2L){
      r <- c(x %% 3, r)
      x <- floor(x / 3)
    }
    r <- c(1, r)
  }
  
  r

}
