compnames <- function(nobj){
#################################################################
#  function to generate a vector of labels
#  for comparisons in sinclair ordering
#################################################################
  nobj <- get("nobj", get("ENV", environment(patt.design)))
  str  <- NULL
  
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      str <- c(str, paste("u", i, j, sep = ""))
    }
  }
  assign("undecnames", str, get("ENV", environment(patt.design)))
  
  str

}
