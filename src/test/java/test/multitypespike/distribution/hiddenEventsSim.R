set.seed(42)

############################### SINGLE-TYPE ####################################

# Simulate birth-death trajectory
simulate_bd <- function(t_max = 1, birth_rate, death_rate, init_pop = 1) {
  time <- 0
  pop <- init_pop
  times <- c(time)
  pops <- c(pop)
  
  while (time < t_max && pop > 0) {
    rate <- birth_rate * pop + death_rate * pop
    # Sample time to next event
    dt <- rexp(1, rate)
    time <- time + dt
    if(time > t_max) {
      times <- c(times, t_max)
      pops <- c(pops, pop)
      break
    }
    
    # Choose event and update state
    if (runif(1) < birth_rate / (birth_rate + death_rate)) {
      pop <- pop + 1  # birth
    } else {
      pop <- pop - 1  # death
    }
    times <- c(times, time)
    pops <- c(pops, pop)
  }
  # Record state
  data.frame(time = times, population = pops)
}


lambda <- 3.0
mu <- 0.5
rho <- 0.2

# Simulate and plot
# bd_data <- simulate_bd(birth_rate = lambda, death_rate = mu)
# plot(bd_data$time, bd_data$population, type = "s", col = "blue",
#      xlab = "Time", ylab = "Population Size", main = "Birth-Death Trajectory")



hidden_events_df <- numeric(nMax)
nMax <- 100000
k <- 1

while(k <= nMax){
  bd_data <- simulate_bd(birth_rate = lambda, death_rate = mu)
  # Rho sampling: reject if rho samples > 1
  m <- rbinom(n=1, size = bd_data[nrow(bd_data),2], prob = rho)
  if(m != 1) {
    next
  } else {
    birth_indices <- which(diff(bd_data$population) > 0) + 1
    hidden_events <- 0
    for (i in birth_indices) {
      Nt <- bd_data$population[i]
      pAnc <- 2/Nt # Probability of birth event being ancestral to the sample
      if(runif(1) < pAnc) hidden_events <- hidden_events + 1
    }
    hidden_events_df[k] <- hidden_events
    k <- k + 1
  }
}


mean(hidden_events_df)


############################### MULTI-TYPE #####################################


# Simulate multi-type birth-death trajectory
simulate_bd_multitype <- function(
    birthA, deathA, birthB, deathB,
    migAB = 0, migBA = 0,
    birthAB = 0, birthBA = 0,
    startTypePriorProbA, t_max = 1) {
  
  time <- 0
  popA <- rbinom(1, 1, startTypePriorProbA)
  popB <- 1 - popA
  
  # Preallocate event storage
  events <- vector("list", 1000)
  i <- 1
  
  while (time < t_max && (popA || popB) > 0) {
    rates <- c(
      birthA * popA, deathA * popA,
      birthB * popB, deathB * popB,
      migAB * popA, migBA * popB,
      birthAB * popA, birthBA * popB
    )
    
    total_rate <- sum(rates)
    if (total_rate == 0) print("Warning total rate = 0")
    
    dt <- rexp(1, total_rate)
    next_time <- time + dt
    
    if (next_time > t_max) break  # Stop before exceeding t_max
    
    event <- sample.int(8, 1, prob = rates)
    
    # Apply event
    switch(event,
           { popA <- popA + 1 },  # birthA
           { popA <- popA - 1 },  # deathA
           { popB <- popB + 1 },  # birthB
           { popB <- popB - 1 },  # deathB
           { popA <- popA - 1; popB <- popB + 1 },  # migAB
           { popB <- popB - 1; popA <- popA + 1 },  # migBA
           { popB <- popB + 1 },  # birthAB
           { popA <- popA + 1 }   # birthBA
    )
    
    # Negative populations
    if (popA < 0 || popB < 0) print("Warning negative population")
    
    time <- next_time
    events[[i]] <- c(time, event, popA, popB)
    i <- i + 1
    
    # Expand list if needed
    if (i > length(events)) events <- c(events, vector("list", 1000))
  }
  
  # Add the final state at exactly t_max
  events[[i]] <- c(t_max, NA, popA, popB)
  n_events <- i
  
  mat <- do.call(rbind, events[1:n_events])
  colnames(mat) <- c("time", "event", "popA", "popB")
  as.data.frame(mat)
}

# -------------- No birth among demes -------------- #

lambdaA <- 3.0
lambdaB <- 3.0
muA <- 0.5
muB <- 0.5
mAB <- 0.9 
mBA <- 0.23
lambdaAB <- 0.0
lambdaBA <- 0.0
ProbA <- 0.5
rho <- 0.2
nMax <- 100000


# Simulate and plot
multitype_bd_data <- simulate_bd_multitype(
  birthA = lambdaA, deathA = muA,
  birthB = lambdaB, deathB = muB,
  migAB = mAB, migBA = mBA,
  birthAB = lambdaAB, birthBA = lambdaBA, startTypePriorProbA = ProbA
)

plot(multitype_bd_data$time, multitype_bd_data$popA, type = "s", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Population Size",
     ylim = range(c(multitype_bd_data$popA, multitype_bd_data$popB)),
     main = "Multi-Type Birth-Death Trajectory")
lines(multitype_bd_data$time, multitype_bd_data$popB, type = "s", col = "red", lwd = 2)
legend("topleft", legend = c("Type A", "Type B"),
       col = c("blue", "red"), lwd = 2, bty = "n")


hidden_events_typeA <- numeric(nMax)
hidden_events_typeB <- numeric(nMax)
k <- 1

while (k <= nMax) {
  sim <- simulate_bd_multitype(
    birthA = lambdaA, deathA = muA,
    birthB = lambdaB, deathB = muB,
    migAB = mAB, migBA = mBA,
    birthAB = lambdaAB, birthBA = lambdaBA,
    startTypePriorProbA = ProbA
  )
  
  # Final populations
  final_popA <- sim$popA[nrow(sim)]
  final_popB <- sim$popB[nrow(sim)]
  
  # Conditioning on a single tip in type A
  m <- rbinom(1, size = final_popA, prob = rho)
  if (m != 1) next  # Reject this trajectory if not 1 sampled tip
  current_type <- "A"
  
  hidden_eventsA <- 0
  hidden_eventsB <- 0
  
  # Traverse events backwards in time
  for (i in seq(nrow(sim), 1)) {
    event <- sim$event[i]
    Nt_A <- sim$popA[i]
    Nt_B <- sim$popB[i]
    
    if(is.na(event)) next
    
    if (event == 5 && current_type == "B") {  # migAB
      if (runif(1) < 1 / Nt_B)
        current_type <- "A"
    } 
    else if (event == 6 && current_type == "A") {  # migBA
      if (runif(1) < 1 / Nt_A)
        current_type <- "B"
    } 
    else if (event == 1 && current_type == "A") {  # birthA
      if (runif(1) < 2 / Nt_A)
        hidden_eventsA <- hidden_eventsA + 1
    } 
    else if (event == 3 && current_type == "B") {  # birthB
      if (runif(1) < 2 / Nt_B)
        hidden_eventsB <- hidden_eventsB + 1
    }
  }
  
  hidden_events_typeA[k] <- hidden_eventsA
  hidden_events_typeB[k] <- hidden_eventsB
  k <- k + 1
}

# Summaries
mean(hidden_events_typeA)
mean(hidden_events_typeB)



# ---------------------------- Birth among demes ----------------------------- #


lambdaA <- 3.0
lambdaB <- 3.0
muA <- 0.5
muB <- 0.5
mAB <- 0.2
mBA <- 0.4
lambdaAB <- 1.5
lambdaBA <- 1.0
ProbA <- 0.5
rho <- 0.2
nMax <- 100000


# Simulate and plot
multitype_bd_data <- simulate_bd_multitype(
  birthA = lambdaA, deathA = muA,
  birthB = lambdaB, deathB = muB,
  migAB = mAB, migBA = mBA,
  birthAB = lambdaAB, birthBA = lambdaBA, startTypePriorProbA = ProbA
)

plot(multitype_bd_data$time, multitype_bd_data$popA, type = "s", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Population Size",
     ylim = range(c(multitype_bd_data$popA, multitype_bd_data$popB)),
     main = "Multi-Type Birth-Death Trajectory")
lines(multitype_bd_data$time, multitype_bd_data$popB, type = "s", col = "red", lwd = 2)
legend("topleft", legend = c("Type A", "Type B"),
       col = c("blue", "red"), lwd = 2, bty = "n")



hidden_events_typeA <- numeric(nMax)
hidden_events_typeB <- numeric(nMax)
k <- 1

# Hidden events of type A refers to hidden speciation events that occur along a lineage of type A. 
# If a birth event occurs from type A to type B, then this is considered a hidden event of type A.
while (k <= nMax) {
  sim <- simulate_bd_multitype(
    birthA = lambdaA, deathA = muA,
    birthB = lambdaB, deathB = muB,
    migAB = mAB, migBA = mBA,
    birthAB = lambdaAB, birthBA = lambdaBA,
    startTypePriorProbA = ProbA
  )
  
  # Final population size
  final_popA <- sim$popA[nrow(sim)]

  # Conditioning on a single tip in type A
  m <- rbinom(1, size = final_popA, prob = rho)
  if (m != 1) next  # Reject this trajectory if not 1 sampled tip
  current_type <- "A"
  
  hidden_eventsA <- 0
  hidden_eventsB <- 0
  
  # Traverse events backwards in time
  for (i in seq(nrow(sim), 1)) {
    event <- sim$event[i]
    Nt_A <- sim$popA[i]
    Nt_B <- sim$popB[i]
    
    if(is.na(event)) next
    
    if (event == 5 && current_type == "B") {  # migAB
      if (runif(1) < 1 / Nt_B) current_type <- "A"
    } 
    else if (event == 6 && current_type == "A") {  # migBA
      if (runif(1) < 1 / Nt_A) current_type <- "B"
    } 
    else if (event == 1 && current_type == "A") {  # birthA
      if (runif(1) < 2 / Nt_A) hidden_eventsA <- hidden_eventsA + 1
    } 
    else if (event == 3 && current_type == "B") {  # birthB
      if (runif(1) < 2 / Nt_B) hidden_eventsB <- hidden_eventsB + 1
    }
    else if (event == 7 && current_type == "B") {  # birthAB
      if (runif(1) < 1 / Nt_B) {
        current_type <- "A"
        }
    }
    else if (event == 8 && current_type == "A") {  # birthBA
      if (runif(1) < 1 / Nt_A) {
        current_type <- "B"
        }
    }
  }
  
  hidden_events_typeA[k] <- hidden_eventsA
  hidden_events_typeB[k] <- hidden_eventsB
  k <- k + 1
}


# Summaries
mean(hidden_events_typeA)
mean(hidden_events_typeB)







