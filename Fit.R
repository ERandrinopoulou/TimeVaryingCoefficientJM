library(rjags)
library(R2WinBUGS)
library(JMbayes)

################################################
# Simulate data
source("Simulate.R")

################################################
# Prepare data
source("PrepareData.R")

################################################
# Save model as txt
source("ModelJAGS.R")
filename <- file.path("VCJM.txt")
write.model(model, filename)

################################################
# run model in jags
model.fit <- jags.model(file = "VCJM.txt", data = Data, n.chains = con$n.chains, 
                        n.adapt = con$n.adapt, quiet = con$quiet)


update(model.fit, con$n.burnin)
res <- coda.samples(model.fit, parms,  n.iter = con$n.iter - con$n.burnin, thin = con$n.thin)
codaFit <- as.mcmc.list(res)

################################################
# results
bss <- do.call(rbind,codaFit)
colnames(bss)
n.sims <- nrow(bss)
sims.list <- vector("list", length(parms))
names(sims.list) <- parms
for (p in seq_along(parms)) {
  ii <- grep(paste("^", parms[p], sep = ""), colnames(bss))
  sims.list[[p]] <- bss[, ii]
}

save(codaFit, file = "results.RData")

######### sims.list is a list with the results
# Obtain the mean/median/mode of each chain
# e.g.
apply(sims.list[[1]], 2, mean)

