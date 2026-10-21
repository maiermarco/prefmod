# aus gtools
perm <- function(n, r, v = seq_len(n)) {
  if(r == 1L) matrix(v, n, 1L) else if(n == 1L) matrix(v, 1L, r) else {
    X <- NULL
    for(i in seq_len(n)) X <- rbind(X, cbind(v[i], Recall(n - 1L, r - 1L, v[-i])))
    X
  }
}
