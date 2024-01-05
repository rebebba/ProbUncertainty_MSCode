## last modified 22/07/2023

## Run with three arguments:
## first:  number of simulations
## second: experimental design (2 or 3)
## third:  proportion in the treatment group (probably 0.2 or 0.5)

## the four runs for the main article are:
## Rscript probCompMainScript.R 1000 2 0.2
## Rscript probCompMainScript.R 1000 2 0.5
## Rscript probCompMainScript.R 1000 3 0.2
## Rscript probCompMainScript.R 1000 3 0.5


install_github("rebebba/ProbUncertainity", build_vignettes = TRUE)

args <- commandArgs(trailingOnly = TRUE)
nSim <- as.numeric(args[1])
type <- as.numeric(args[2])
propTreatment <- as.numeric(args[3])

cat("\n\n")
cat("Running with ")
cat(nSim)
cat(" replicates of each scenario.\n")
cat("Experimental design is type: ")
cat(type)
cat(".\n")
cat("Proportion in treatment group: ")
cat(propTreatment)
cat(".\n\n")


pdfFigs <- TRUE

# experimental design (2 or 3)
type <- 2

# is the total sample size is split evenly between two groups?
propTreatment <- 0.5
#propTreatment <- 0.2

# a range of sample sizes
nVals <- c(10, 20, 30, 50, 70, 100)

# true values of the proportion of successes for the default group
p0Vals <- c(0.1, 0.5)

# difference in proportion of successes between the default group and the other group
deltapVals <- round(seq(-1, 1, by = 0.1), 1)

# any given run will be of this many replicate simulations of
# each combination of the above parameters
nSim <- 1000

# combine all the parameters
d <- expand.grid(
  n = nVals,
  p0 = p0Vals,
  deltap = deltapVals,
  rep = 1:nSim)

# remove all those that would correspond to expected
# proportions of successes outside of the permissible range
d <- subset(d, d$p0 + d$deltap >= 0 & d$p0 + d$deltap <= 1)

## THINGS TO KEEP TRACK OF
# p values from chi square test, with and without Yates' continuitiy correction, 
# also Fisher's exact test
d$chiSqP <- d$chiSqYP <- d$FisherExactP <- NA
# Berger and Boos correction for unconditional exact test
d$uncondExact_P <- NA
# estimate, SE and p for lm
d$lmEst <- d$lmSE <- d$lmP <- NA
# estimate, SE and p for glm
d$glmEst <- d$glmSE <- d$glmP <- NA
# glm estimate and standard error converted to the expected data scale
d$glmEstData <- d$glmSEData <- NA
# p value from LRT and profile CI bounds
d$glmLRTP <- d$profCI_L <- d$profCI_U <- NA
# Agresti confidence interval
d$scoreCI_L <- d$scoreCI_U <- NA
# for z test and associated CI on binomial variance
d$binomialSE <- NA

# OR and its CI from fisher.test
d$fisher.testCI_L <- d$fisher.testCI_U <- NA

library(epitools)
# three flavours of P-value provided by oddsratio in epitools
d$oddsratio_midp.exact_P <- d$oddsratio_fisher.exact_P <- d$oddsratio_chi.square.exact_P <- NA

library(exact2x2)
# estimate, CI and P_value from exact2x2
d$exact2x2_Est <- d$exact2x2_CI_L <- d$exact2x2_CI_U <- NA
d$exact2x2_P <- NA

cat("\n Size of simulation results dataframe:\n")
dim(d)
cat("\n")

# For type 3 experimental design:
# simulate random x [0,1] but make sure there are >=2 observations for each class of x
sim.x <- function(n, p) {
  x <- rbinom(n, 1, p)
  if (sum(x) == 0)
    x[1:2] <- 1
  if (sum(x) == 1 & x[1] == 0)
    x[1] <- 1
  if (sum(x) == 1 & x[2] == 0)
    x[2] <- 1
  if (sum(x) == n)
    x[1:2] <- 0
  if (sum(x) == (n - 1) & x[1] == 1)
    x[1] <- 0
  if (sum(x) == (n - 1) & x[2] == 1)
    x[2] <- 0
  return(x)
}





## do all simulations

for (i in 1:(dim(d)[1])) {
  # it takes a long time, so this output will be appreciated!
  if (i %% 1000 == 0)
    print(paste(
      "Done sim ",
      i,
      " of ",
      dim(d)[1],
      ", ",
      round(100 * i / dim(d)[1], 1),
      "% done."
    ))

  # simulate data where one group has a probability of p0
  # and the other has a probability of p0+deltap, with total
  # sample size of n
  x <- c(rep(0, d$n[i] * (1 - propTreatment)), rep(1, d$n[i] * propTreatment))
  if (type == 3)
    x <- sim.x(n = d$n[i], p = propTreatment)
  simDat <- data.frame(x = x)
  simDat$y <- rbinom(d$n[i], 1, d$p0[i] + simDat$x * d$deltap[i])

  # contingency table for the chi square and FETs, save the p-values
  ct <- table(factor(simDat$x, levels = 0:1), factor(simDat$y, levels = 0:1))

  # chi square and FETs, save the p values
  suppressWarnings(d$chiSqP[i] <- chisq.test(ct, correct = FALSE)$p.value)
  suppressWarnings(d$chiSqYP[i] <- chisq.test(ct, correct = TRUE)$p.value)
  d$FisherExactP[i] <- fisher.test(ct)$p.value
  d[i, c("fisher.testCI_L", "fisher.testCI_U")] <- fisher.test(ct)$conf.int
  
  # Berger and Boos correction for unconditional exact test using exact2x2 package
  n1 <- nrow(simDat[simDat$x == 0,])
  n2 <- nrow(simDat[simDat$x == 1,])
  x1 <- nrow(simDat[simDat$x == 0 & simDat$y == 0,])
  x2 <- nrow(simDat[simDat$x == 1 & simDat$y == 0,])
  d$uncondExact_P[i] <- exact2x2::uncondExact2x2(x1, n1, x2, n2, gamma = 10e-6)$p.value

  # linear model on the counts as a function of the treatment
  # save the estimated difference in probabilities, the standard error, and p value
  lm1 <- lm(y ~ x, data = simDat)
  d[i, c("lmEst", "lmSE", "lmP")] <-
    suppressWarnings(summary(lm1)$coefficients[2, c(1, 2, 4)])

  # logistic glm on the counts as a function of the treatment
  # save the estimated difference in probabilities (logit scale),
  # the standard error, and the p value
  glm1 <- glm(y ~ x, data = simDat, family = "binomial")
  d[i, c("glmEst", "glmSE", "glmP")] <- summary(glm1)$coefficients[2, c(1, 2, 4)]
  # translation from the logit / log-odds scale to the probability scale
  # then calculate difference in probabilities
  d$glmEstData[i] <- glm.prob.diff(glm1)
  # se, delta method
  d$glmSEData[i] <- glm.se.data(glm1)

  # null model for LRT
  glm0 <- glm(y ~ 1, data = simDat, family = "binomial")
  d$glmLRTP[i] <- as.numeric(1 - pchisq(2 * (logLik(glm1) - logLik(glm0)), 1))
  # profile likelihood CI
  d[i, c("profCI_L", "profCI_U")] <- profCI(x = simDat$x, y = simDat$y)

  # score CI
  d[i, c("scoreCI_L", "scoreCI_U")] <- scoreCI(x = simDat$x, y = simDat$y)

  # binomial sampling SE in delta p
  d$binomialSE[i] <-
    sqrt((ct[1, 2] / (ct[1, 1] + ct[1, 2]) * (1 - ct[1, 2] / (ct[1, 1] + ct[1, 2]))) /
           (ct[1, 1] + ct[1, 2]) +
           (ct[2, 2] / (ct[2, 1] + ct[2, 2]) * (1 - ct[2, 2] /
                                                  (ct[2, 1] + ct[2, 2]))) / (ct[2, 1] + ct[2, 2]))

  simDat$y <- factor(simDat$y, levels = c(0, 1))
  t <- table(simDat$x, simDat$y)

  # three flavours of P-value provided by oddsratio in epitools
  d[i, c(
    "oddsratio_midp.exact_P",
    "oddsratio_fisher.exact_P",
    "oddsratio_chi.square.exact_P"
  )] <-
    suppressWarnings(oddsratio(
      x = simDat$x,
      y = simDat$y,
      method = "fisher"
    )$p.value[2, ])

  e2x2 <- exact2x2(x = simDat$x, y = simDat$y)
  d$exact2x2_Est[i] <- e2x2$estimate
  d[i, c("exact2x2_CI_L", "exact2x2_CI_U")] <- e2x2$conf.int
  d$exact2x2_P[i] <- e2x2$p.value

}




## save batches of simulation results and their associated parameters

simParams <- list(
  type = type,
  propTreatment = propTreatment,
  nSim = nSim,
  nVals = nVals,
  p0Vals = p0Vals,
  deltapVals = deltapVals
)

save(
  d,
  simParams,
  file = paste(
    "~/GitHub/probabiliy-comparisons/saveFile_type",
    type,
    "_nsim",
    nSim,
    "_pt",
    propTreatment,
    ".RData",
    sep = ""
  )
)
