aprior <- function(gamma.hat){
  m <- mean(gamma.hat)
  v <- var(gamma_hat)
  output <- (2 * v + m^2) / v
  return(output)
}
