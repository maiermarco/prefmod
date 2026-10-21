checkMIS <- function(obj, nitems, MISmodel = "obj", obj.names = NULL, verbose = FALSE){
  
  nobj <- nitems
  ## prepare data
  ncomp <- nobj * (nobj - 1L) / 2
  if(is.character(obj)){                    ## datafilename supplied
    datafile <- obj
    if(file.access(datafile, mode = 0) == 0){
      dat <- as.matrix(read.table(datafile, header = TRUE))  # datafile
    } else {
      stop("\ninput data file does not exist!\n")
    }
  } else if(is.data.frame(obj)){            ## data frame supplied
    dat <- as.matrix(obj)                                  # dataframe
    dat <- apply(dat, 2L, as.numeric)
  } else {
    stop("first argument must be either datafilename or dataframe")
  }
  
  varnames <- colnames(dat)
  if (ncol(dat) > ncomp) {
    covnames <- varnames[(ncomp + 1L):ncol(dat)]
    covs     <- as.data.frame(dat[, (ncomp + 1L):ncol(dat)])
  } else {
    covs <- NULL
  }
  
  ## vector of counts of NAs, columnwise
  misv   <- apply(dat, 2L, function(x) sum(is.na(x)) )
  mism   <- matrix(0, nobj, nobj)
  mism[upper.tri(mism, diag = FALSE)] <- misv
  mism   <- mism + t(mism)
  objMis <- colSums(mism)
  
  if (is.null(obj.names)){
    rownames(mism) <- colnames(mism) <- paste("o", seq_len(nobj), sep = "")
    names(objMis)  <- paste("o", seq_len(nobj), sep = "")
  } else {
    rownames(mism) <- colnames(mism) <- obj.names
    names(objMis)  <- obj.names
  }
  objMis <- colSums(mism)
  if(verbose){
    cat("number of missing comparisons:\n")
    print(mism)
    cat("number of missing comparisons for objects:\n", objMis, "\n")
  }
  if(MISmodel == "obj"){
    RET <- ifelse(objMis == 0, FALSE, TRUE)
  } else if(MISmodel == "comp"){
    RET <- ifelse(misv == 0, FALSE, TRUE)
  } else {
    stop("\nNo output: MISmodel not correctly specified. Use \"obj\" or \"comp\"\n")
  }
  
  invisible(RET)

}
