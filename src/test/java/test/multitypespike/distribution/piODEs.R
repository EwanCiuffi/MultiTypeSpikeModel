
library(ape)
library(deSolve)




# ODE system
p0ge_ode <- function(t, y, parms) {
  nTypes <- parms$nTypes
  interval <- parms$interval
  b <- parms$b
  d <- parms$d
  s <- parms$s
  b_ij <- parms$b_ij
  M <- parms$M
  
  yDot <- numeric(2 * nTypes)
  
  for (i in 1:nTypes) {
    # p0 equations
    yDot[i] <- (b[interval, i] + d[interval, i] + s[interval, i] -
                  b[interval, i] * y[i]) * y[i] - d[interval, i]
    
    for (j in 1:nTypes) {
      if (i == j) next
      
      yDot[i] <- yDot[i] +
        b_ij[[interval]][i, j] * y[i] -
        b_ij[[interval]][i, j] * y[i] * y[j] +
        M[interval, i, j] * y[i] -
        M[interval, i, j] * y[j]
    }
    
    # ge equations
    ge_index <- nTypes + i
    yDot[ge_index] <- (b[interval, i] + d[interval, i] + s[interval, i] -
                         2 * b[interval, i] * y[i]) * y[ge_index]
    
    for (j in 1:nTypes) {
      if (i == j) next
      
      yDot[ge_index] <- yDot[ge_index] +
        b_ij[[interval]][i, j] * y[ge_index] -
        b_ij[[interval]][i, j] * (y[i] * y[nTypes + j] + y[j] * y[ge_index]) +
        M[interval, i, j] * y[ge_index] -
        M[interval, i, j] * y[nTypes + j]
    }
  }
  
  return(list(yDot))
}


# Parameter list
parms <- list(
  nTypes = n_types,
  interval = 1,
  b = birth_rate,
  d = death_rate,
  s = sampling_rate,
  b_ij = birth_rate_among,
  M = migration_rate
)

# Initial conditions

p0 <- rep(1, n_types)
ge <- c(0.5, 0.0)
y0 <- c(p0 = p0, ge = ge)

# Time grid
times <- seq(0, 1 , length.out = 100)

# Solve
p0ge <- ode(y = y0, times = times, func = p0ge_ode, parms = parms)
colnames(out)[2:5] <- c("p0_0", "p0_1", "ge_0", "ge_1")


# Parse Newick tree
newick <- "t1[&state=0]:1.0;"
tree <- read.tree(text = newick)

# Tip states
tip_states <- c(0, 1)  # t1 = 0, t2 = 1
names(tip_states) <- tree$tip.label

# Origin time
origin <- 0

# Model parameters
n_types <- 2
birth_rate <- matrix(rep(2.0, n_types), nrow = 1)         # [interval, type]
death_rate <- matrix(rep(1.0, n_types), nrow = 1)
sampling_rate <- matrix(rep(0.5, n_types), nrow = 1)
removal_prob <- matrix(rep(1.0, n_types), nrow = 1)

# Migration matrix: [interval, i, j]
migration_rate <- array(0, dim = c(1, n_types, n_types))
migration_rate[1,,] <- matrix(c(0.0, 0.2, 0.1, 0.0), nrow = n_types)

# Birth between types: [interval][i,j]
birth_rate_among <- list(2)
birth_rate_among[[1]] <- matrix(0, nrow = n_types, ncol = n_types)
birth_rate_among[[2]] <- matrix(0, nrow = n_types, ncol = n_types)


# ODE system
pi_ode <- function(t, y, parms) {
  nTypes <- parms$nTypes
  interval <- parms$interval
  b <- parms$b
  d <- parms$d
  s <- parms$s
  b_ij <- parms$b_ij
  M <- parms$M
  
  yDot <- numeric(nTypes)
  
  for (i in 1:nTypes) {
    yDot[i] = 0;
    
    for (j in 1:nTypes) {
      if (j == i) next;
      
      yDot[i] <-  ((b_ij[[interval]][i, j] * p0ge[j] + M[interval, i, j]) * (p0ge[nTypes + i] / max(p0ge[nTypes + j], 1e-12))) * y[j];
      - ((b_ij[[interval]][i, j] * p0ge[i] + M[interval, i, j]) * (p0ge[nTypes + j] / max(p0ge[nTypes + i], 1e-12))) * y[i];
    }
    
  }
  
  return(list(yDot))
}


# Parameter list
parms <- list(
  nTypes = n_types,
  interval = 1,
  b = birth_rate,
  d = death_rate,
  s = sampling_rate,
  b_ij = birth_rate_among,
  M = migration_rate
)

# Initial conditions

y0 <- c(0.5247090305225015, 0.4752909694774985)
# ge <- c(0.5, 0.0),
# y0 <- c(p0 = p0, ge = ge)

# Time grid
# times <- seq(2.5, 1 , length.out = 50)

# Solve
out <- ode(y = y0, times = times, func = pi_ode, parms = parms)
colnames(out)[2:5] <- c("p0_0", "p0_1", "ge_0", "ge_1")
