## ============================================================================
## Calculates the exact multi-type marginal log-density for two hardcoded branches 
## to use as the ground-truth expected value in the unit test.
## ============================================================================

exact_multitype_marginal_logdensity <- function(spike, expNrHiddenEvents, pi_val, spikeShape, isFakeParent, kmax = 2000) {
  k <- 0:kmax
  logpk <- dpois(k, lambda = expNrHiddenEvents, log = TRUE)
  
  iszero <- spike < 1e-9
  
  # Initialize density vectors
  dens_obs0 <- numeric(length(k))
  dens_obs1 <- numeric(length(k))
  
  # Number of spikes for each condition
  nSpikes_obs0 <- k
  nSpikes_obs1 <- if (isFakeParent) k else k + 1
  
  # --- Condition 1: Observed Event == 0 (Weight = 1 - pi) ---
  zero_idx0 <- nSpikes_obs0 == 0L
  if (iszero) {
    dens_obs0[zero_idx0]  <- 1.0
    dens_obs0[!zero_idx0] <- 0.0
  } else {
    dens_obs0[zero_idx0]  <- 0.0
    dens_obs0[!zero_idx0] <- dgamma(spike, shape = spikeShape * nSpikes_obs0[!zero_idx0], rate = spikeShape, log = FALSE)
  }
  
  # --- Condition 2: Observed Event == 1 (Weight = pi) ---
  zero_idx1 <- nSpikes_obs1 == 0L
  if (iszero) {
    dens_obs1[zero_idx1]  <- 1.0
    dens_obs1[!zero_idx1] <- 0.0
  } else {
    dens_obs1[zero_idx1]  <- 0.0
    dens_obs1[!zero_idx1] <- dgamma(spike, shape = spikeShape * nSpikes_obs1[!zero_idx1], rate = spikeShape, log = FALSE)
  }
  
  # Combine probabilities using the mixture weights (pi and 1-pi)
  prob_k <- exp(logpk) * ((1 - pi_val) * dens_obs0 + pi_val * dens_obs1)
  
  log(sum(prob_k))
}

# --- Shared Parameters ---
spikeShape <- 2.0

# --- Node 0 (Taxon 1) ---
logP0_type0 <- exact_multitype_marginal_logdensity(spike = 0.8, expNrHiddenEvents = 0.5, pi_val = 0.8, spikeShape = spikeShape, isFakeParent = FALSE)
logP0_type1 <- exact_multitype_marginal_logdensity(spike = 0.3, expNrHiddenEvents = 0.2, pi_val = 0.6, spikeShape = spikeShape, isFakeParent = FALSE)

# --- Node 1 (Taxon 2) ---
logP1_type0 <- exact_multitype_marginal_logdensity(spike = 2.5, expNrHiddenEvents = 1.2, pi_val = 0.7, spikeShape = spikeShape, isFakeParent = FALSE)
logP1_type1 <- exact_multitype_marginal_logdensity(spike = 1.1, expNrHiddenEvents = 0.9, pi_val = 0.4, spikeShape = spikeShape, isFakeParent = FALSE)

# --- Node 2 (Root pseudo-priors) ---
# Multi-type assigns pseudo-priors to ALL types at the root
logP_root_type0 <- dgamma(0.5, shape = 2.0, scale = 0.5, log = TRUE)
logP_root_type1 <- dgamma(0.5, shape = 2.0, scale = 0.5, log = TRUE)


total_logP <- logP0_type0 + logP0_type1 + logP1_type0 + logP1_type1 + logP_root_type0 + logP_root_type1

cat(sprintf("Node 0 Type 0: %.8f\n", logP0_type0))
cat(sprintf("Node 0 Type 1: %.8f\n", logP0_type1))
cat(sprintf("Node 1 Type 0: %.8f\n", logP1_type0))
cat(sprintf("Node 1 Type 1: %.8f\n", logP1_type1))
cat(sprintf("Root Type 0: %.8f\n", logP_root_type0))
cat(sprintf("Root Type 1: %.8f\n", logP_root_type1))
cat(sprintf("---------------------------\n"))
cat(sprintf("Total Expected: %.8f\n", total_logP))