pattR.fit <- function(
  obj,
  nitems,
  formel    = ~ 1,
  elim      = ~ 1,
  resptype  = "ranking",
  obj.names = NULL,
  ia        = FALSE,
  NItest    = FALSE,
  pr.it     = FALSE
){
  
  call <- match.call()
  if(ia) cat("Warning:\n\tDependencies do not make sense for rankings! \n")

  ENV       <- new.env()
  ENV$pr.it <- pr.it

  ENV$resptype <- "ranking"

  nobj      <- nitems
  ENV$Rnobj <- nobj   # number of ranked objects
  opt       <- options()
  options("warn" = -1L)
  
#######################################
  if(is.character(obj)){                                     # datafilename supplied
    datafile <- obj
    if(file.access(datafile, mode = 0) == 0){
      dat <- as.matrix(read.table(datafile, header = TRUE))  # datafile
    } else {
      stop("\ninput data file does not exist!\n")
    }
  } else if(is.data.frame(obj)){                             # data frame supplied
    dat <- as.matrix(obj)                                    # dataframe
    dat <- apply(dat, 2L, as.numeric)
  } else {
    stop("first argument must be either datafilename or dataframe")
  }
  
##### formel.names <- attr(terms(as.formula(formel)), "term.labels")
##### elim.names   <- attr(terms(as.formula(elim)), "term.labels")
##### covnames     <- unique(c(formel.names, elim.names))
  
  varnames <- colnames(dat)
  if(ncol(dat) > nobj){
    # covnames <- varnames[(nobj + 1L):ncol(dat)]
    # covs     <- as.data.frame(dat[, (nobj + 1L):ncol(dat)])
    ## instead of the above rh 2011-05-13
    formel.names <- attr(terms(as.formula(formel)), "term.labels")
    formel.names <- unique(unlist(strsplit(formel.names, ":", fixed = TRUE))) # 2011-08-31 remove interaction terms
    elim.names   <- attr(terms(as.formula(elim)), "term.labels")
    elim.names   <- unique(unlist(strsplit(elim.names, ":", fixed = TRUE)))     # 2011-08-31 remove interaction terms
    covnames     <- unique(c(formel.names, elim.names))
    covs         <- as.data.frame(dat[, covnames])
  } else {
    covs <- NULL
  }
  
  # for rankings: at least one item not NA
  idx <- apply(dat[, seq_len(nobj)], 1L, function(x){ sum(!is.na(x)) > 0L })
  dat <- dat[idx, ]
  
  # check for proper ranks
  nodups <- apply(dat[, seq_len(nobj)], 1L, function(x){ max(table(x)) == 1L })
  dups   <- !nodups
  if(sum(dups) > 0L){
    norankslines <- seq_len(nrow(dat))[dups]
    cat("Warning:\n\timproper ranks in lines", norankslines, " - removed from data\n")
    dat <- dat[nodups, ]
    # if(!is.null(covs)) covs <- covs[nodups, ] ## replaced 20-08-09
    if(!is.null(covs)) covs <- covs[nodups, , drop = FALSE]
  }
  
  # transform into PCs
  dat <- ifelse(is.na(dat), 99999L, dat) # if removed only comparisons between chosem then reverse results in -lambda
  
  pc.dat <- NULL
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      pc.dat <- cbind(pc.dat, as.numeric(dat[, i] > dat[, j]))
    }
  }
  pc.notnaidx <- NULL
  for(j in seq.int(2L, nobj)){
    for(i in seq_len(j - 1L)){
      pc.notnaidx <- cbind(pc.notnaidx, dat[, i] != dat[, j])   # which values are not missing -> TRUE
    }
  }
  pc.na <- ifelse(pc.notnaidx, 1, NA)   # set FALSE to NA
  dat   <- pc.dat * pc.na
  
  # dat must be data frame
  dat <- data.frame(dat)
  
  # check for NAs in subject covariates
  if(!is.null(covs)){
    covs <- as.data.frame(covs)
    colnames(covs) <- covnames
    NAs <- which(!complete.cases(covs))                  # check for NA
    if(length(NAs) > 0L){
      cat("Warning:\n\tsubject covariates: NAs in lines", NAs, " - removed from data\n")
      notNAs <- which(complete.cases(covs))
      dat    <- dat[notNAs, ]
      # covs <- covs[notNAs, ] ## replaced 20-08-09
      covs <- covs[notNAs, , drop = FALSE]
    }
  }
  
######################################
  
  ## option for obj.names added
  if(is.null(obj.names)) ENV$obj.names <- varnames[seq_len(nobj)] else ENV$obj.names <- obj.names[seq_len(nobj)]
  
  if(NItest)
    if(!anyNA(dat))
      stop("Test for ignorable missing cannot be performed - no NA values!")
  
  ENV$Y <- ifelse(Rpatternmat(nobj) > 0, 1, -1) # pattern matrix
  
  # no undecided with rankings
  ENV$undec <- FALSE
  
  
  ENV$NItest <- NItest
  if(ENV$NItest){
    if(formel != "~1" || elim != "~1"){
      covs   <- NULL
      formel <- ~ 1
      elim   <- ~ 1
      cat("\ncurrently no covariates fitted if NItest == TRUE !!\n")
    }
  }
  
  ENV$ia <- ia
  if(ia){
    depL     <- dependencies(nobj, ENV$Y)
    ENV$XI   <- depL$d
    ilabels  <- depL$label.intpars
    npars.ia <- nobj * (nobj - 1L) * (nobj - 2L) / 2
  } else {
    ilabels  <- NULL
    npars.ia <- 0
  }
  
  X <- ENV$Y %*% pcdesign(nobj)
  X <- -X[, -nobj]                          # basic design matrix

  cList     <- splitCovs(dat, covs, formel, elim, ENV)     # split data according to subject covariates
  partsList <- gen.partsListR(nobj, cList, ENV)     # generate list for all subj covariate x miss values groups
  rm(cList)
  
  
  npar <- (nobj - 1L) * ENV$ncovpar + ENV$undec + npars.ia
  if(ENV$NItest) npar <- (nobj - 1L) * 2L + ENV$undec + npars.ia
  
  lambda   <- rep(0, npar)
  ENV$iter <- 0
  
  ## MAIN FITTING ROUTINE
  result <- nlm(loglik, lambda, X, nobj, partsList, ENV, hessian = TRUE, iterlim = 1000L)
  
  options(opt)
  
  ENV$nobj    <- nobj
  ENV$ilabels <- ilabels
  
  envList   <- mget(ls(ENV), envir = ENV)
  outputobj <- list(
      coefficients = result$estimate,
      ll           = ENV$ll,
      fl           = ENV$fl,
      call         = call,
      result       = result,
      envList      = envList,
      partsList    = partsList
    )
  class(outputobj) <- "pattMod"                         #class: pattern model
  
  outputobj

}
