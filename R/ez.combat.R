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

ez.combat <- function(df,
                      batch.var,
                      adjust.var = "all",
                      exclude.var = NULL,
                      model = NULL,
                      opt = list(out.opt = c("overwrite", "append"),
                                 err.opt = c("abort", "continue"),
                                 eb = TRUE,
                                 verbose = FALSE)) {

  orig.f <- df # keep df in original state for output

  # Gather and check variable list ---------------------------------------------
  if (!is.null(model)) {
    iv.ls <- unlist(strsplit(as.character(model), split="[~+*: ]"))
    iv.ls <- iv.ls[-which(iv.ls == "")]
    for (i in 1:length(iv.ls)) {
      df[ ,iv.ls[i]] <- as.numeric(df[ ,iv.ls[i]])
    }
  } else {
    iv.ls <- NULL
  }

  model <- switch(class(model),
                  `NULL` = NULL,
                  `character` = model <- model.matrix(as.formula(model), df),
                  `formula` = model <- model.matrix(as.formula(model), df),
                  otherwise = error("[EZ combat] Cannot parse model"))

  if (adjust.var == "all") {
    if (!is.null(iv.ls)) {
      dv.ls <- colnames(df)[!(colnames(df) %in% iv.ls)]
    } else {
      dv.ls <- colnames(df)
    }
  } else {
    dv.ls <- adjust.var
  }
  dv.ls <- dv.ls[!dv.ls %in% batch.var]

  if (!is.null(exclude.var)) {
    dv.ls <- dv.ls[!dv.ls %in% exclude.var]
  }
  # Convert df to dat ----------------------------------------------------------
  dat <- t(as.matrix(df[ ,dv.ls]))

  # Check and remove constant variables (or abort) -----------------------------
  sd.chk <- logical(length(dv.ls))
  for (i in 1:length(dv.ls)) {
    sd.chk[i] <- sd(df[ ,dv.ls[i]]) == 0
  }
  if (any(sd.chk)) {
    if (opt$err.opt[1] == "continue") {
      dv.ls <- dv.ls[-which(sd.chk)]
      if (opt$verbose) {
        warning(sprintf("[EZ.combat] %s has been excluded as it is constant across samples.\n", dv.ls))
      }
    } else {
      stop(sprintf("[EZ.combat] aborted: %s is constant across samples.\n", dv.ls))
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

  # Check if the design is confounded ------------------------------------------
  if(qr(design)$rank<ncol(design)){
    if(ncol(design)==(n.batch+1)){
      stop("[EZ.combat] The covariate is confounded with batch. Remove the covariate and rerun ComBat.")
    }
    if(ncol(design)>(n.batch+1)){
      if((qr(design[,-c(1:n.batch)])$rank<ncol(design[,-c(1:n.batch)]))){
        stop('[EZ.combat] The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.')
      } else {
        stop("[EZ.combat] At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.")
      }
    }
  }

  ##Standardize Data across features
  if (opt$verbose) cat('[EZ.combat] Standardizing Data across features\n')
  B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(dat))

  #Standarization Model
  grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
  var.pooled <- ((dat-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)
  stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
  if(!is.null(design)){
    tmp <- design;tmp[,c(1:n.batch)] <- 0
    stand.mean <- stand.mean+t(tmp%*%B.hat)
  }
  s.data <- (dat-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

  ##Get regression batch effect parameters
  if (opt$eb){
    if (opt$verbose) cat("[EZ.combat] Fitting L/S model and finding priors\n")
  } else {
    if (opt$verbose) cat("[EZ.combat] Fitting L/S model\n")
  }

  batch.design <- design[,1:n.batch]
  gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))

  delta.hat <- NULL
  for (i in batches){
    delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var,na.rm=T))
  }

  # Empirical Bayes correction:
  gamma.star <- delta.star <- NULL
  gamma.bar <- t2 <- a.prior <- b.prior <- NULL
  if (opt$eb){
    ##Find Priors
    gamma.bar <- apply(gamma.hat, 1, mean)
    t2 <- apply(gamma.hat, 1, var)
    a.prior <- apply(delta.hat, 1, aprior)
    b.prior <- apply(delta.hat, 1, bprior)


    ##Find EB batch adjustments
    if (opt$verbose) cat("[EZ.combat] Finding parametric adjustments\n")
    for (i in 1:n.batch){
      temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
      gamma.star <- rbind(gamma.star,temp[1,])
      delta.star <- rbind(delta.star,temp[2,])
    }
  }

  ### Normalize the Data ###
  if (opt$verbose) cat("[EZ.combat] Adjusting the Data\n")
  bayesdata <- s.data
  j <- 1
  for (i in batches){
    if (opt$eb){
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
    } else {
      bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.hat))/(sqrt(delta.hat[j,])%*%t(rep(1,n.batches[j])))
    }
    j <- j+1
  }

  bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean

  if (opt$out.opt[1] == "append") {
    new.df <- data.frame(orig.f, t(bayesdata))
    colnames(orig.f) <- c(colnames(orig.f), paste0(dv.ls, ".cb"))
  } else if (opt$out.opt[1] == "overwrite") {
    new.df <- orig.f
    new.df[ ,dv.ls] <- t(bayesdata)
  }

  return(list(df=new.df,
              gamma.hat=gamma.hat, delta.hat=delta.hat,
              gamma.star=gamma.star, delta.star=delta.star,
              gamma.bar=gamma.bar, t2=t2, a.prior=a.prior, b.prior=b.prior, batch=batch, mod=model,
              stand.mean=stand.mean, stand.sd=sqrt(var.pooled)[,1])
  )
}
