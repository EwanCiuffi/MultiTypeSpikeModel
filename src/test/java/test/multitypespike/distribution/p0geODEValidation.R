library(ape)
library(deSolve)

setwd("~/code/beast_and_friends/MultiTypeSpikeModel/src/test/multitypespike/distribution")

# Load BDMM-prime p0ge output
lines <- readLines("p0ge_output.csv")
# Remove header
lines <- lines[-1]
# Split each line by tab and convert to numeric
split_data <- strsplit(lines, "\t")
bdmm_matrix <- do.call(rbind, lapply(split_data, as.numeric))
# Convert to data frame
bdmm_p0ge <- as.data.frame(bdmm_matrix)
colnames(bdmm_p0ge) <- c("time", "p0_0", "p0_1", "ge_0", "ge_1")



# Parse Newick tree
newick <- "(t1[&state=0]:1.5, t2[&state=1]:0.5);"
tree <- read.tree(text = newick)

# Tip states
tip_states <- c(0, 1)  # t1 = 0, t2 = 1
names(tip_states) <- tree$tip.label

# Origin time
origin <- 2.5

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
times <- seq(2.5, 1 , length.out = 50)

# Solve
out <- ode(y = y0, times = times, func = p0ge_ode, parms = parms)
colnames(out)[2:5] <- c("p0_0", "p0_1", "ge_0", "ge_1")

# print(out)

# Plot results
# matplot(out[, "time"], out[, -1], type = "l", lty = 1, col = 1:4,
#         xlab = "Time", ylab = "Values", main = "p0 and ge trajectories")
# legend("topright", legend = c("p0_0", "p0_1", "ge_0", "ge_1"), col = 1:4, lty = 1)


# Compare ge ratio trajectories

# Compute ratio
ge_ratio_r <- out[, "ge_0"] / out[, "ge_1"]

ge_ratio_bdmm <- bdmm_p0ge$ge_0 / bdmm_p0ge$ge_1

## Plot ge_0/ge_1 ratio
plot(out[, "time"], ge_ratio_r, type = "l", col = "darkgreen", lwd = 2,
     xlab = "Time", ylab = "Values", main = "ge_0 / ge_1 ratio")

lines(out[, "time"], ge_ratio_bdmm, col = "orange", lwd = 2, lty = 2)

legend("topleft", legend = c("R output", "BDMM-prime output"),
       col = c("darkgreen", "orange"), lwd = 2, lty = c(1, 2))



# Compare p0 trajectories

p0_0_r <- out[, "p0_0"]
p0_1_r <- out[, "p0_1"]

p0_0_bdmm <- bdmm_p0ge$p0_0
p0_1_bdmm <- bdmm_p0ge$p0_1

# Plot p0_0
plot(out[, "time"], p0_0_r, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Values", main = "p0_0")

lines(out[, "time"], p0_0_bdmm, col = "red", lwd = 2, lty = 2)

legend("topleft", legend = c("R output", "BDMM-prime output"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))


# Plot p0_1
plot(out[, "time"], p0_1_r, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Values", main = "p0_1")

lines(out[, "time"], p0_1_bdmm, col = "red", lwd = 2, lty = 2)

legend("topleft", legend = c("R output", "BDMM-prime output"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))


#######################
## Birth among demes ##
#######################


# Load BDMM-prime p0ge output
lines <- readLines("p0ge_output_birthAmongDemes.csv")
# Remove header
lines <- lines[-1]
# Split each line by tab and convert to numeric
split_data <- strsplit(lines, "\t")
bdmm_matrix <- do.call(rbind, lapply(split_data, as.numeric))
# Convert to data frame
bdmm_p0ge <- as.data.frame(bdmm_matrix)
colnames(bdmm_p0ge) <- c("time", "p0_0", "p0_1", "ge_0", "ge_1")



# Parse Newick tree
newick <- "(t1[&state=0]:1.5, t2[&state=1]:0.5);"
tree <- read.tree(text = newick)

# Tip states
tip_states <- c(0, 1)  # t1 = 0, t2 = 1
names(tip_states) <- tree$tip.label

# Origin time
origin <- 2.5

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
birth_rate_among <- list(1)
birth_rate_among[[1]] <- matrix(c(0.0, 1.0, 1.5, 0.0), nrow = n_types)


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
times <- seq(2.5, 1 , length.out = 50)

# Solve
out <- ode(y = y0, times = times, func = p0ge_ode, parms = parms)
colnames(out)[2:5] <- c("p0_0", "p0_1", "ge_0", "ge_1")

# print(out)

# Plot results
# matplot(out[, "time"], out[, -1], type = "l", lty = 1, col = 1:4,
#         xlab = "Time", ylab = "Values", main = "p0 and ge trajectories")
# legend("topright", legend = c("p0_0", "p0_1", "ge_0", "ge_1"), col = 1:4, lty = 1)


# Compare ge ratio trajectories

# Compute ratio
ge_ratio_r <- out[, "ge_0"] / out[, "ge_1"]

ge_ratio_bdmm <- bdmm_p0ge$ge_0 / bdmm_p0ge$ge_1

## Plot ge_0/ge_1 ratio
plot(out[, "time"], ge_ratio_r, type = "l", col = "darkgreen", lwd = 2,
     xlab = "Time", ylab = "Values", main = "ge_0 / ge_1 ratio")

lines(out[, "time"], ge_ratio_bdmm, col = "orange", lwd = 2, lty = 2)

legend("topleft", legend = c("R output", "BDMM-prime output"),
       col = c("darkgreen", "orange"), lwd = 2, lty = c(1, 2))



# Compare p0 trajectories

p0_0_r <- out[, "p0_0"]
p0_1_r <- out[, "p0_1"]

p0_0_bdmm <- bdmm_p0ge$p0_0
p0_1_bdmm <- bdmm_p0ge$p0_1

# Plot p0_0
plot(out[, "time"], p0_0_r, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Values", main = "p0_0")

lines(out[, "time"], p0_0_bdmm, col = "red", lwd = 2, lty = 2)

legend("topleft", legend = c("R output", "BDMM-prime output"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))


# Plot p0_1
plot(out[, "time"], p0_1_r, type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Values", main = "p0_1")

lines(out[, "time"], p0_1_bdmm, col = "red", lwd = 2, lty = 2)

legend("topleft", legend = c("R output", "BDMM-prime output"),
       col = c("blue", "red"), lwd = 2, lty = c(1, 2))



###########
# 
# 
# # Tip states
# tip_states <- 0
# names(tip_states) <- "t1"
# 
# # Origin time
# origin <- 1
# 
# # Model parameters
# n_types <- 2
# birth_rate <- matrix(rep(2.0, n_types), nrow = 1)         # [interval, type]
# death_rate <- matrix(rep(1.0, n_types), nrow = 1)
# sampling_rate <- matrix(c(0.2,0), nrow = 1)
# removal_prob <- matrix(rep(1.0, n_types), nrow = 1)
# 
# # Migration matrix: [interval, i, j]
# migration_rate <- array(0, dim = c(1, n_types, n_types))
# migration_rate[1,,] <- matrix(c(0.0, 0.0, 0.0, 0.0), nrow = n_types)
# 
# # Birth between types: [interval][i,j]
# birth_rate_among <- list(1)
# birth_rate_among[[1]] <- matrix(c(0.0, 1.5, 0.0, 0.0), nrow = n_types)
# 
# 
# # Parameter list
# parms <- list(
#   nTypes = n_types,
#   interval = 1,
#   b = birth_rate,
#   d = death_rate,
#   s = sampling_rate,
#   b_ij = birth_rate_among,
#   M = migration_rate
# )
# 
# # Initial conditions
# 
# p0 <- c(1.0, 0.0)
# ge <- c(0.2, 0.0)
# y0 <- c(p0 = p0, ge = ge)
# 
# # Time grid
# times <- seq(1, 0 , length.out = 50)
# 
# # Solve
# out <- ode(y = y0, times = times, func = p0ge_ode, parms = parms)
# colnames(out)[2:5] <- c("p0_0", "p0_1", "ge_0", "ge_1")
# out
