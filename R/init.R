init <- function(dfr){
  ##############################################################################
  # initialisations, read data                                                 #
  ##############################################################################
  
  datafile  <- get("datafile",  get("ENV", environment(patt.design)))
  nitems    <- get("nitems",    get("ENV", environment(patt.design)))
  objnames  <- get("objnames",  get("ENV", environment(patt.design)))
  cov.sel   <- get("cov.sel",   get("ENV", environment(patt.design)))
  blnRevert <- get("blnRevert", get("ENV", environment(patt.design)))
  resptype  <- get("resptype",  get("ENV", environment(patt.design)))

  if(!resptype %in% c("paircomp", "rating", "ranking"))
    stop("resptype must be one of \"paircomp\", \"rating\", \"ranking\"")

  if(is.null(dfr)){
    if(file.access(datafile, mode = 0) == 0){
      dat <- as.matrix(read.table(datafile, header = TRUE))  # datafile
    } else {
      stop("\ninput data file does not exist!\n")
    }
  } else {
     dat <- as.matrix(dfr)                                  # dataframe
     dat <- apply(dat, 2L, as.numeric)
  }
  NAs <- which(!complete.cases(dat))                  # check for NA
  if(length(NAs) > 0L){
    cat("\tNAs in lines", NAs, " - removed from data\n")
    notNAs <- which(complete.cases(dat))
    dat    <- dat[notNAs, ]
  }
  
  nsubj    <- nrow(dat)
  ncolumns <- ncol(dat)
  nobj     <- nitems
  ncov     <- ncolumns - nitems
  ncomp    <- nobj * (nobj - 1L) / 2         #    number of comparisons
  
  ncatL <- diff(range(dat[, seq_len(nobj)])) + 1   # number of rating(Likert) categories
  
  if(resptype == "paircomp"){                 # in case of PC
    ncatPC <- ncatL                         # number of response categories
    # ncomp <- nitems
    npatt <- NULL
    # nobj <- 0.5 + sqrt(0.25 + 2 * nitems) #    inverse of choose(nobj, 2L)
    objnames <- get("objnames", get("ENV", environment(patt.design)))
    if(length(objnames) != nobj)            ### default objnames
      objnames <- paste("o", seq_len(nobj), sep = "")
  } else {                                  # in case of  ratings/likert or rankings
    if(length(objnames) != nobj)            ### obj names from data if not defined in call rh 2011-08-24
      objnames <- colnames(dat)[seq_len(nobj)]
    ncomp  <- nobj * (nobj - 1) / 2         #    number of comparisons
    npatt  <- ncatL**nobj                   #    for ratings/likert
    ncatPC <- ncatL * 2 - 1                 #    number of categories for differences
  }

  blnSubjcov <- cov.sel[1L] != ""           #    FALSE if ""

  if(blnSubjcov){
    if(resptype == "paircomp"){         # in case of PC
      inpcovnames <- colnames(dat)[(ncomp + 1L):ncolumns]
    } else {
      inpcovnames <- colnames(dat)[(nobj + 1L):ncolumns]
    }
    if(toupper(cov.sel[1L]) == "ALL"){   # all covariates included
      cov.sel <- inpcovnames
      if(cov.sel[1L] == "") stop("\nno subject covariates in data\n") #rh 2001-08-12
    } else if(length(setdiff(cov.sel, inpcovnames)) > 0L){
      stop("\nsubject covariate name(s) incorrectly specified\n")
    }
    cov.case <- as.matrix(dat[, c(cov.sel)])
    colnames(cov.case) <- cov.sel
    covlevels <- apply(cov.case, 2L, max)
    covlevmin <- apply(cov.case, 2L, min)         # check for subj covs which have no proper values
    covwrong  <- !(covlevels > 1 & covlevmin == 1) # don't have at least 2 values (highest lev > 1) and don't start with 1
    if(any(covwrong)){
      wrongcov <- colnames(cov.case)[covwrong]
      txt      <- paste0("\ncovariates with improper values: ", wrongcov, "\n")
      stop(txt)
    }
    ncov <- length(cov.sel)
  } else {
    cov.case  <- NULL
    covlevels <- NULL
    covnames  <- NULL
    ncov      <- 0
  }
  
  reverse <- ifelse(blnRevert, -1, 1)
  
  ret <- list(
      dat          = dat,
      nsubj        = nsubj,
      nobj         = nobj,
      ncatL        = ncatL,
      ncatPC       = ncatPC,
      resptype     = resptype,
      reverse      = reverse,
      objnames     = objnames,
      npatt        = npatt,
      ncomp        = ncomp,
      blnSubjcov   = blnSubjcov,
      ncov         = ncov,
      covlevels    = covlevels,
      covnames     = cov.sel,
      cov.case     = cov.case,
      nintcovs.out = 40
    )
  # writes all control structures to environment ENV
  for(i in seq_along(ret))
    do.call("assign", list(names(ret)[i], ret[[i]], get("ENV", environment(patt.design)) ))

}
