\name{ez.combat_internal_functions}
\alias{aprior}
\alias{bprior}
\alias{it.sol}
\alias{postmean}
\alias{postvar}
\title{Internal functions, not to be called by the user directly}
\description{Internal functions, not to be called by user directly}
\usage{
    aprior(gamma_hat)
    bprior(gamma_hat)
    it.sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=0.0001)
    postmean(g_hat, g_bar, n, d_star, t2)
    postvar(sum2, n, a, b)
}
\author{Timothy Koscik}
\keyword{internal}
