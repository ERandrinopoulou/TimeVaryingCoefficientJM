model <- function ()
{
  for (i in 1:N) {
    # Longitudinal Part
    for (j in offset[i]:(offset[i+1] - 1)) {
      muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 1:ncZ], Z[j, 1:ncZ])
      y[j] ~ dnorm(muy[j], tau)
    }
    # Survival Part
    etaBaseline[i] <- inprod(gammas[1:ncW], W[i, 1:ncW])
    log.h0.T[i] <-  inprod(Bs.gammas[1:ncWBH], WBH[i, 1:ncWBH])
    f.T[i] <- inprod(betas[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 1:ncZ], Ztime[i, 1:ncZ])
    

    f.Lam[i] <- inprod(alphas[1:ncF], Lam[i, 1:ncF])
    
    
  
    log.hazard[i] <- log.h0.T[i] + etaBaseline[i] + f.Lam[i] * f.T[i]
    for (k in 1:K) {
      log.h0.s[i, k] <- inprod(Bs.gammas[1:ncWBH], WBHs[K * (i - 1) + k, 1:ncWBH])
      f.s[i, k] <- inprod(betas[1:ncX], Xs[K*(i - 1) + k, 1:ncX]) + inprod(b[i, 1:ncZ], Zs[K*(i - 1) + k, 1:ncZ])
    
      
      f.Lams[i, k] <- inprod(alphas[1:ncF], Lam.s[K * (i - 1) + k, 1:ncF])
 
      SurvLong[i, k] <- wk[k] * exp(log.h0.s[i, k] + f.Lams[i, k] * f.s[i, k])
    }
    log.survival[i] <- - exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ])
    phi[i] <- C - (event[i] * log.hazard[i]) - log.survival[i]
    zeros[i] ~ dpois(phi[i])
    # Random Effects Part
    b[i, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
  }
  # Priors
  # Longitudinal Part
  betas[1:ncX] ~ dmnorm(priorMean.betas[], priorTau.betas[, ])
  tau ~ dgamma(priorA.tau, priorB.tau)
  # Survival Part
  gammas[1:ncW] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, ])
  
  alphas[1:ncF] ~ dmnorm(priorMean.alphas[], Tau.alphas * priorTau.alphas[, ])
  Tau.alphas ~ dgamma(1, 0.005)

  Bs.gammas[1:ncWBH] ~ dmnorm(priorMean.Bs.gammas[], Tau.Bs.gammas * priorTau.Bs.gammas[, ])
  Tau.Bs.gammas ~ dgamma(1, 0.005)


  # Random Effects Part
  inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
