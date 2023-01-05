postmean <- function(g.hat,g.bar,n,d.star,t2){
  output <- (t2 * n * g_hat + d_star * g_bar) / (t2 * n + d_star)
  return(output)
}
