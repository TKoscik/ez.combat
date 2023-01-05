it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
  n <- apply(!is.na(sdat), 1, sum)
  g_old <- g_hat
  d_old <- d_hat
  change <- 1
  while (change > conv) {
    g_new  <- postmean(g_hat, g_bar, n, d_old, t2)
    sum2   <- apply((sdat - g_new %*% t(rep(1, ncol(sdat))))^2, 1, sum, na.rm=TRUE)
    d.new  <- postvar(sum2, n, a, b)
    change <- max(abs(g_new - g_old) / g_old, abs(d.new - d_old) / d_old)
    g_old <- g_new
    d_old <- d.new
  }
  adjust <- rbind(g_new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  return(adjust)
}
