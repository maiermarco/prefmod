# matrix of all likert patterns
Lpatternmat <- function(datrng, nobj){
  
  Lpatt <- all_patterns(datrng[1L], datrng[2L], nobj) # for ordinal models only
  diffs <- diffsred(Lpatt, nobj)
  diffs <- unique(diffs)
  diffs

}
