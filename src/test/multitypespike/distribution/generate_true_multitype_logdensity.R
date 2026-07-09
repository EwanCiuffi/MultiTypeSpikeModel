## ============================================================================
## Calculates the exact multi-type joint marginal log-density for a given branch
## across all its types to use as the ground-truth expected value in the unit test.
## ============================================================================

exact_multitype_node_logdensity <- function(spikes, expNrHiddenEvents, pi_vals, spikeShapes, isFakeParent, kmax = 2000) {
  n_types <- length(spikes)
  logP0 <- numeric(n_types)
  logP1 <- numeric(n_types)
  
  for (i in 1:n_types) {
    spike <- spikes[i]
    expHidden <- expNrHiddenEvents[i]
    shape <- spikeShapes[i]
    iszero <- spike < 1e-9
    
    k <- 0:kmax
    logpk <- dpois(k, lambda = expHidden, log = TRUE)
    pk <- exp(logpk)
    
    # --- P0: No observed event (nSpikes = k) ---
    prob0 <- 0.0
    if (iszero) {
      prob0 <- prob0 + pk[1] # k=0
    } else {
      # k > 0
      gamma_dens0 <- dgamma(spike, shape = shape * k[-1], rate = shape, log = FALSE)
      prob0 <- prob0 + sum(pk[-1] * gamma_dens0)
    }
    
    # --- P1: 1 observed event (nSpikes = k + 1) ---
    prob1 <- 0.0
    if (!iszero) {
      gamma_dens1 <- dgamma(spike, shape = shape * (k + 1), rate = shape, log = FALSE)
      prob1 <- sum(pk * gamma_dens1)
    }
    
    logP0[i] <- if (prob0 > 0) log(prob0) else -Inf
    logP1[i] <- if (prob1 > 0) log(prob1) else -Inf
  }
  
  if (isFakeParent) {
    return(sum(logP0))
  } else {
    # Joint combinatorial probability across types
    logTerms <- numeric(n_types)
    for (i in 1:n_types) {
      if (pi_vals[i] > 0) {
        term <- log(pi_vals[i]) + logP1[i]
        for (j in 1:n_types) {
          if (j != i) {
            term <- term + logP0[j]
          }
        }
        logTerms[i] <- term
      } else {
        logTerms[i] <- -Inf
      }
    }
    
    maxLogTerm <- max(logTerms)
    if (maxLogTerm == -Inf) {
      return(-Inf)
    } else {
      sumExp <- sum(exp(logTerms - maxLogTerm))
      return(maxLogTerm + log(sumExp))
    }
  }
}

# --- Shared Parameters ---
spikeShape <- 2.0
spikeShapes <- c(spikeShape, spikeShape) 

# --- Node 0 (Taxon 1) ---
spikes_node0  <- c(0.8, 0.3)
exp_node0     <- c(0.5, 0.2)
pi_node0      <- c(0.8, 0.6)
logP_node0 <- exact_multitype_node_logdensity(spikes_node0, exp_node0, pi_node0, spikeShapes, isFakeParent = FALSE)

# --- Node 1 (Taxon 2) ---
spikes_node1  <- c(2.5, 1.1)
exp_node1     <- c(1.2, 0.9)
pi_node1      <- c(0.7, 0.4)
logP_node1 <- exact_multitype_node_logdensity(spikes_node1, exp_node1, pi_node1, spikeShapes, isFakeParent = FALSE)

# --- Node 2 (Root pseudo-priors) ---
# Hardcoded to Gamma(1.0, 1.0) / Exponential(1)
logP_root_type0 <- dgamma(0.5, shape = 1.0, rate = 1.0, log = TRUE) 
logP_root_type1 <- dgamma(0.5, shape = 1.0, rate = 1.0, log = TRUE)
logP_root_total <- logP_root_type0 + logP_root_type1

total_logP <- logP_node0 + logP_node1 + logP_root_total

cat(sprintf("Total Expected: %.8f\n", total_logP))