all_patterns <- function(a, b, nobj = NULL){
##########################################################
# corresponds to: rev( expand.grid(a:b, a:b, a:b, ...) ) #
#   like GLIM's $gfac                                    #
##########################################################s
  
  if(is.null(nobj)) nobj <- get("nobj", get("ENV", environment(patt.design)))
  
  grid <- vector(mode = "list", length = nobj)
  for(i in seq_len(nobj)) grid[[i]] <- a:b
  
  patterns <- rev(expand.grid(grid))
  
  patterns

}
