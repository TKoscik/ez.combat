\name{ez.combat}
\alias{ez.combat}
\title{ComBat Harmonization for Dataframe Objects}
\description{A easy to use function for applying ComBat Harmonization for batch effects on dataframe objects.}
\usage{
  ez.combat(df, batch.var, adjust.var = "all", exclude.var = NULL, model = NULL, opt = list(out.opt = c("overwrite", "append"), err.opt = c("abort", "continue"), eb = TRUE, verbose = FALSE)) 
}
\arguments{
  \item{df}{Dataframe containing variables to be harmonized as well as the batch variable to harmonize for.}
  \item{batch.var}{A character string or integer indicating the dataframe column containing the batch variable.}
  \item{adjust.var}{A character vector indicating the names of the variables in the dataframe to be harmonized. Default value is "all", indicating that all variables should be harmonized (excluding the batch variable and any independent variables given in a model}
  \item{exclude.var}{An optional character vector indicating which variables to exclude. Best used in conjunction with adjust.var = "all" if the number of variables to be excluded is fewer than the number to be included}
  \item{model}{a character vector specifying the model to be considered during harmonization. Models do not need to specify the dependent variable, and should follow R's standard modelling syntax, e.g., "~ group * age"}
  \item{opt}{A list containing additional options for ComBat Harmonization:
    opt$out.opt: Chnages the behaviour of the output dataframe in the output list. "overwrite" will replace the unadjusted variables with adjusted ones, "append" will add the adjusted variables to the end of the output dataframe with ".cb" appended to the variable name.
    opt$err.opt: Changes behaviour for handling certain input errors. "continue" harmonization will continue if varaibles have to be excluded, "abort" will discontinue harmonization
    opt$eb: Logical, Perform ComBat Harmonization with empirical bayes or without
    opt$verbose: Logical, print progress messages to console}
}
\value{A list containing:
    df: a dataframe with adjusted values
    gamma.hat and delta hat: Estimated location and shift (L/S) parameters before empirical Bayes.
    gamma.star and delta.star: Empirical Bayes estimated L/S parameters.
    gamma.bar, t2, a.prior and b.prior: esimated prior distributions parameters.
}
\author{
  Timothy R. Koscik <timkoscik-ezcombat@gmail.com>
}
    
