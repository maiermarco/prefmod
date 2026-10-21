pattRrep.fit <- function(
  obj,
  nitems,
  tpoints   = 1,
  formel    = ~ 1,
  elim      = ~ 1,
  resptype  = "rankingT",
  obj.names = NULL,
  ia        = FALSE,
  iaT       = FALSE,
  NItest    = FALSE,
  pr.it     = FALSE
){
  
  call <- match.call()
  if(ia) cat("Warning:\n\tDependencies do not make sense for rankings! \n")
  if(tpoints < 2) stop("no of timepoints incorrectly specified! if tpoints == 1 use pattR.fit")
  
  ENV <- new.env()
  ENV$pr.it <- pr.it
  
  ENV$resptype <- "ratingT"
  
  nobj      <- nitems * tpoints
  ENV$Rnobj <- nobj   # number of ranked objects
  opt       <- options()
  options("warn" = -1L)

  ENV$nitems  <- nitems
  ENV$tpoints <- tpoints
  
  ### check data file
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
  
  ### formula and variables
  formel.names <- attr(terms(as.formula(formel)), "term.labels")
  elim.names   <- attr(terms(as.formula(elim)), "term.labels")
  covnames     <- unique(c(formel.names, elim.names))
  
  varnames     <- colnames(dat)
  
### variables rightmost columns in data
  if(ncol(dat) > nobj){
    #covnames <- varnames[(nobj + 1L):ncol(dat)]
    #covs <- as.data.frame(dat[, (nobj + 1L):ncol(dat)])
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
  
  # check for proper ranks
  dups   <- apply(dat[, seq_len(tpoints * nitems)], 1L, function(x){ max(table(x)) > tpoints })
  nodups <- !dups
  
  if(sum(dups) > 0){
    norankslines <- seq_len(nrow(dat))[dups]
    cat("Warning:\n\timproper ranks in lines", norankslines, " - removed from data\n")
    dat <- dat[nodups, ]
    if(!is.null(covs)) covs <- covs[nodups, , drop = FALSE]
  }
  
  
  # transform into PCs
  dat <- ifelse(is.na(dat), 99999L, dat) # if removed with partial rankings: only comparisons between chosen then reverse results in -lambda
  ## reduce datamatrix
  ## for each timepoint
  #dat.t <- NULL
  dat.t <- as.data.frame(diffsred(dat[, seq_len(nitems)], nitems))
  if(tpoints > 1){
    for(t in seq.int(2L, tpoints)){
      from  <- nitems * (t - 1L) + 1L
      to    <- from + nitems - 1L
      dat.t <- cbind(dat.t, as.data.frame(diffsred(dat[, from:to], nitems)))
    }
  }
  ## replace 0 with NA
###browser()
  dat.t <- lapply(dat.t, function(x){ ifelse(x == 0, NA, x) })
  ###dat <- dat.t
  pc.dat <- data.frame(dat.t)
  rm(dat.t)
  
  
  # check for NAs in subject covariates
  if(!is.null(covs)){
    covs <- as.data.frame(covs)
    colnames(covs) <- covnames
    NAs <- which(!complete.cases(covs))                  # check for NA
    if(length(NAs) > 0){
      cat("Warning:\n\tsubject covariates: NAs in lines", NAs, " - removed from data\n")
      notNAs <- which(complete.cases(covs))
      dat    <- dat[notNAs, ]
      #covs <- covs[notNAs, ] ## replaced 20-08-09
      covs <- covs[notNAs, , drop = FALSE]
    }
  }
  
  ## option for obj.names added
  if(is.null(obj.names))
    ENV$obj.names <- varnames[seq_len(nobj)]
  else
    ENV$obj.names <- obj.names[seq_len(nobj)]
  ####################################  end data, inits

  ##design matrix
  ENV$Y <- ifelse(Rpatternmat(nitems) > 0, 1, -1) # pattern matrix
  
  # "recursively" expand patternmatrix for timepoints
  np  <- nrow(ENV$Y)                     # no patterns in Y
  npp <- np
  YL  <- ENV$Y
  for(t in seq_len(tpoints - 1L)){
    YL  <- do.call("rbind", lapply(seq_len(np), function(i) YL)) # stacks YL np times
    YR  <- expand.mat(ENV$Y, rep(npp, np))                       # repeats each line of Y np^t times
    YL  <- cbind(YL, YR)
    npp <- npp * np
  }
  ENV$Y <- YL
  rm(YR, YL) # tidy up
  if(NItest)
    if(!anyNA(pc.dat))
      stop("Test for ignorable missing cannot be performed - no NA values!")
  
  
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
  
  ncomp <- choose(nitems, 2L)
  X     <- -(ENV$Y[, seq_len(ncomp)] %*% pcdesign(nitems))[, -nitems]
  for(t in seq.int(2L, tpoints)){
    from <- ncomp * (t - 1L) + 1L
    to   <- from + ncomp - 1L
    X    <- cbind(X, -(ENV$Y[, from:to] %*% pcdesign(nitems))[, -nitems])
  }
  
  
  # dependence parameters between timepoints (AR(1))
  ENV$iaT <- iaT
  if(iaT){
    npars.iaT <- ncomp * (tpoints - 1L)
    ENV$XIT   <- do.call("cbind", lapply(seq_len(npars.iaT), function(i){ ENV$Y[, i] * ENV$Y[, i + ncomp] }))
    ENV$iTlabels <- paste(
        paste("Comp", seq_len(ncomp), sep = ""),
        paste("IT", rep(seq_len(tpoints - 1L), rep(ncomp, tpoints - 1L)), rep(seq.int(2L, tpoints), rep(ncomp, tpoints - 1L)), sep = ""),
      sep = ":")
  } else {
    ENV$iTlabels <- NULL
    npars.iaT    <- 0
  }
  
  
  
  #X <- ENV$Y %*% pcdesign(nobj)
  #X <- -X[, -nobj]                          # basic design matrix
  
  cList     <- splitCovs(pc.dat, covs, formel, elim, ENV)     # split data according to subject covariates
  partsList <- gen.partsListR(nitems, cList, ENV)     # generate list for all subj covariate x miss values groups
  rm(cList)
  
  npar <- tpoints * (nitems - 1L) * ENV$ncovpar + ENV$undec * tpoints + npars.ia + npars.iaT
  if(ENV$NItest) npar <- tpoints * (nitems - 1L) * 2L + ENV$undec * tpoints + npars.ia + npars.iaT
  
  #npar <- (nobj - 1L) * ENV$ncovpar + ENV$undec + npars.ia
  #if(ENV$NItest) npar <- (nobj - 1L) * 2L + ENV$undec + npars.ia
  
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
