convert2strings <- function(pattmat){
################################
# initializes resulting string #
# uses p2string                #
################################
    
  strvec <- NULL
  for(i in seq_len(nrow(pattmat))) strvec[i] <- p2string(pattmat[i, ])
  
  strvec

}
