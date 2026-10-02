## ============================================================================
## Calculates the exact single-type marginal log-density for two hardcoded branches 
## to use as the ground-truth expected value in the unit test.
## ============================================================================

exact_marginal_logdensity <- function(spike, expNrHiddenEvents, spikeShape, isFakeParent, kmax = 2000) {
  k <- 0:kmax
  logpk <- dpois(k, lambda = expNrHiddenEvents, log = TRUE)
  nSpikes <- if (isFakeParent) k else k + 1
  
  iszero <- spike < 1e-9
  dens <- numeric(length(k))
  zero_idx  <- nSpikes == 0L
  nzero_idx <- !zero_idx
  
  if (iszero) {
    dens[zero_idx]  <- 1.0
    dens[nzero_idx] <- 0.0
  } else {
    dens[zero_idx]  <- 0.0
    dens[nzero_idx] <- dgamma(spike, shape = spikeShape * nSpikes[nzero_idx],
                              rate = spikeShape, log = FALSE)
  }
  
  total <- sum(exp(logpk) * dens)
  log(total)
}

# Define the exact parameters we will use in the JUnit test
spikeShape <- 2.0

# Branch 0 (Taxon 1)
spike0 <- 0.8
expEvents0 <- 0.5
logP0 <- exact_marginal_logdensity(spike0, expEvents0, spikeShape, isFakeParent = FALSE)

# Branch 1 (Taxon 2)
spike1 <- 2.5
expEvents1 <- 1.2
logP1 <- exact_marginal_logdensity(spike1, expEvents1, spikeShape, isFakeParent = FALSE)

# Add the root pseudo-prior
spike_root <- 0.5
logP_root <- dgamma(spike_root, shape=1.0, scale=1.0, log=TRUE)


total_logP <- logP0 + logP1 + logP_root


cat(sprintf("Branch 0 logP: %.8f\n", logP0))
cat(sprintf("Branch 1 logP: %.8f\n", logP1))
cat(sprintf("Root logP: %.8f\n", logP1))
cat(sprintf("Total Expected logP: %.8f\n", total_logP))


