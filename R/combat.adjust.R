# ------------------------------------------------------------------------------
# Usage and License:
# This code has been modified from the ComBatHarmonization code available from:
# https://github.com/Jfortin1/ComBatHarmonization, which in turn is a
# modification of the original ComBat function code from the sva package that
# can be found at: https://bioconductor.org/packages/release/bioc/html/sva.html
#
# The original, modifications, and present code is under the
# Artistic License 2.0.
# If using this code, make sure you agree and accept this license.
#
# Purpose:  To provide Combat Harmonization with the ease and utility of
#           dataframe structures in R.
# Author:   Timothy R. Koscik, timkoscik+combat@gmail.com
# Date:     July 7, 2019
#-------------------------------------------------------------------------------

combat.adjust <- function(df,
                          batch.var,
                          adjust.var = "all",
                          model = NULL,
                          opt = list(out.opt = c("append", "overwrite"),
                                     err.opt = c("abort", "continue"),
                                     eb = FALSE,
                                     verbose = FALSE)) {

  # Gather and check variable list ---------------------------------------------
  model <- switch(class(model),
                  `NULL` = NULL,
                  `character` = model <- model.matrix(as.formula(model)),
                  `formula` = model <- model.matrix(as.formula(model)),
                  otherwise = error("[EZ combat] Cannot parse model"))

  iv.ls <- colnames(model)[-1]

  if (adjust.var == "all") {
    dv.ls <- colnames(df)[!(colnames(df) %in% iv.ls)]
  } else {
    dv.ls <- adjust.var
  }

  # Check and remove constant variables (or abort) -----------------------------
  sd.chk <- logical(length(dv.ls))
  for (i in 1:length(dv.ls)) {
    sd.chk[i] <- sd(df[ ,dv.ls[i]]) == 0
  }
  if (any(sd.chk)) {
    if (opt$err.opt[1] == "continue") {
      dv.ls <- dv.ls[-which(sd.chk)]
      if (opt$verbose) {
        warning(sprintf("%s has been excluded as it is constant across samples.\n", dv.ls))
      }
    } else {
      stop(sprintf("ComBat aborted: %s is constant across samples.\n", dv.ls))
    }
  }

  # display method option ------------------------------------------------------
  if (opt$eb & opt$verbose) {
    cat("[EZ.combat] Performing ComBat with empirical Bayes\n")
  } else if (opt$verbose) {
    cat("[EZ.combat] Performing ComBat without empirical Bayes (L/S model)\n")
  }

  # make batch a factor and make a set of indicators for batch -----------------
  batch <- as.factor(df[ , batch.var])
  batchmod <- model.matrix(~-1+batch)
  if (opt$verbose) cat("[EZ.combat] Found",nlevels(batch),'batches\n')

  # A few other characteristics on the batches ---------------------------------
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch){ batches[[i]] <- which(batch == levels(batch)[i]) } # list of samples in each batch
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)

  # combine batch variable and covariates --------------------------------------
  design <- cbind(batchmod,model)
  # check for intercept in covariates, and drop if present ---------------------
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[,!check])

  # Number of covariates or covariate levels
  if (opt$verbose) {
    cat("[EZ combat] Adjusting for",
        ncol(design)-ncol(batchmod),
        "covariate(s) or covariate level(s)\n")
  }
}
