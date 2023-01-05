\name{ez.combat}
\alias{ez.combat}
\title{ComBat Harmonization for Dataframe Objects}
\description{A easy to use function for applying ComBat Harmonization for batch effects on dataframe objects.}
\usage{
  ez.combat(df,
            batch.var,
            adjust.var = "all",
            exclude.var = NULL,
            model = NULL,
            output = c("overwrite", "append"),
            error = c("abort", "continue"),
            use.eb = TRUE,
            verbose = FALSE)
}
\arguments{
  \item{df}{Dataframe containing variables to be harmonized as well as the batch variable to harmonize for.}
  \item{batch.var}{A character string or integer indicating the dataframe column containing the batch variable.}
  \item{adjust.var}{A character vector indicating the names of the variables in the dataframe to be harmonized. Default value is "all", indicating that all variables should be harmonized (excluding the batch variable and any independent variables given in a model}
  \item{exclude.var}{An optional character vector indicating which variables to exclude. Best used in conjunction with adjust.var = "all" if the number of variables to be excluded is fewer than the number to be included}
  \item{model}{a character vector specifying the model to be considered during harmonization. Models do not need to specify the dependent variable, and should follow R's standard modelling syntax, e.g., "~ group * age"}
  \item{output}{Select the behaviour of the output dataframe in the output list. "overwrite" will replace the unadjusted variables with adjusted ones, "append" will add the adjusted variables to the end of the output dataframe with ".cb" appended to the variable name.}
  \item{error}{Changes behaviour for handling certain input errors. "continue" harmonization will continue if varaibles have to be excluded, "abort" will discontinue harmonization}
  \item{use.eb}{Logical, Perform ComBat Harmonization with empirical bayes or without}
  \item{verbose}{Logical, print progress messages to console}
}
\value{A list containing:
    df: a dataframe with adjusted values
    gamma.hat and delta hat: Estimated location and shift (L/S) parameters before empirical Bayes.
    gamma.star and delta.star: Empirical Bayes estimated L/S parameters.
    gamma.bar, t2, a.prior and b.prior: esimated prior distributions parameters.
}
\author{
  Timothy R. Koscik <timothy-koscik@uiowa.edu>
}
\examples{
# Harmonize 'iris' data by removing species effects:
cb <- ez.combat(df = iris, batch.var = "Species")

## plots of raw iris data
plot(iris$Sepal.Length, iris$Petal.Length, col = iris$Species, pch = 16,
     xlab = "Sepal Length", ylab = "Petal Length",
     main = "Flower Characteristics by Species")
  legend(x = 4.5, y = 7, legend = levels(iris$Species), col = c(1:3), pch = 16)
boxplot(iris$Sepal.Length ~ iris$Species, notch = T, las = 1,
        xlab = "Species", ylab = "Sepal Length", main = "Sepal Length by Species",
        cex.lab = 1.5, cex.axis = 1.5,cex.main = 2)

## plots after dataset is harmonized acropss specie, i.e, species effects are removed
plot(cb$df$Sepal.Length, cb$df$Petal.Length, col = cb$df$Species,
     pch = 16,
     xlab = "Sepal Length", ylab = "Petal Length",
     main = "Flower Characteristics by Species - Harmonized")
  legend(x = 4.5, y = 7, legend = levels(iris$Species), col = c(1:3), pch = 16)
boxplot(cb$df$Sepal.Length ~ cb$df$Species, notch = T, las = 1,
        xlab = "Species", ylab = "Sepal Length",
        main = "Sepal Length by Species - Harmonized",
        cex.lab = 1.5, cex.axis = 1.5,cex.main = 2)
}
