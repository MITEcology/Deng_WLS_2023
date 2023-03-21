# load necessary packages
library(tidyverse)
library(mvtnorm)
library(mgcv)
library(magrittr)
library(gtools) # permutation with repeats
library(dplyr)
library(geometry)
library(uniformly)
library(pracma)
library(deSolve)

# function that computes the normalized feasibility from an interaction matrix analytically
# inputs: alpha = interaction matrix
#         individual = TRUE: community-level Omega
#                      FALSE: individual-level omega
# output: out = the normalized feasibility
# !!! double check community or species level
Omega <- function(alpha, individual = FALSE) {
  alpha <- as.matrix(alpha)
  S <- nrow(alpha)
  # omega function
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    if(individual){
      out <- d[1]^(1/S) # species level
    }else{
      out <- d[1] # community level
    }
    return(out)
  }
  # rule out errors
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (all(f(alpha) == FALSE)) {
    return(0)
  }
  else {
    Sigma <- solve(t(alpha) %*% alpha)
    return(omega(S, Sigma))
  }
}

# function that normalizes a vector in the L2 norm
normalize <- function(a, norm = 2) {
  a / sqrt(sum(a^2))
}

# function that normalizes the spanning vectors of the feasibility domain in the L2 norm
# inputs: alpha = interaction matrix, diag(alpha) = -1
# output: the normalized spanning vectors
generate_span_vectors <- function(alpha) {
  apply(-alpha, 2, function(x) normalize(x))
}

# function that determine whether a vertex of one cone is inside another cone or not.
# inputs: Span = spanning vectors of the interaction matrix
#         vector = vertex to be checked
# outputs: 1 = inside, 0 = outside
inside_detection <- function(Span, vector) {
  lambda <- solve(Span, vector)
  # check whether vector is a positive linear combination of the column vectors of span
  if (sum(lambda >= -1e-10) == length(lambda)) {
    return(1)
  } else {
    return(0)
  }
}

# function that generates a random matrix (as in May, Nature (1972))
# inputs: num = number of species; stren = standard deviation of interaction strength; 
#         conne = connectance of the interaction matrix
# output: Inte = the generated random matrix
random_interaction_matrix <- function(num, stren = 1, conne = FALSE) {
  Inte <- rnorm(num * num, mean = 0, sd = stren)
  if(conne){
    zeroes <- sample(c(rep.int(1, floor(num * num * conne)), rep.int(0, (num * num - floor(num * num * conne)))))
    Inte[which(zeroes == 0)] <- 0
  }
  Inte <- matrix(Inte, ncol = num, nrow = num)
  diag(Inte) <- -1
  return(Inte)
}

# function that check the global stability
# inputs: A = S by S interaction matrix of the entire community under LV dynamics
# output: TRUE: globally stable
check_global_stability <- function(A){
  all(eigen(A+t(A))$values < 0)
}

# function that generate a random globally stable interaction matrix with `num` species
# inputs: num = number of species; stren = standard deviation of interaction strength
# output: Inte = a random globally stable interaction matrix
random_globally_stable_matrix <- function(num, stren = 1, conne = FALSE){
  Inte <- random_interaction_matrix(num, stren = stren, conne = conne)
  while(!check_global_stability(Inte)){
    Inte <- random_interaction_matrix(num, stren = stren, conne = conne)
  }
  return(Inte)
}

# Function that samples a d-dimensional vector uniformly 
# within or on the d-dimensional hypersphere
# d: dimension of the space to sample from
# positive: whether only vectors in the positive orthant are needed
# within: whether the points should also be sampled inside the hypersphere
hypersphere_sampling <- function(d, positive = FALSE, within = FALSE) {
  x <- rnorm(d, 0, 1)
  if (positive)
    x <- abs(x)
  norm_x <- sqrt(sum(x^2))
  x_scaled <- x / norm_x
  if (within) {
    u <- runif(1, 0, 1)^(1/d)
    x_scaled <- x_scaled * u
  }
  return(x_scaled)
}

# function that sampling initial condition
# input: d = dimension of the space to sample from
#        N0 = values of a non-random initial condition
#        random = TRUE - random initial condition uniformly sampled from [min,max]
#                FALSE - non-random initial condition with values N0
#        min / max = lower / upper bounds for the random initial condition
# output: generated sinitial condition
initial_condition <- function(d, N0 = 0.5, random = FALSE, min = 0, max = 1) {
  if(!random){
    N <- rep(N0, d)
  }else if(random){
    N <- runif(d, min, max)
    N <- N / sum(N)
  }
  return(N)
}

# Solves the system of ordinary differential equations
# given by the generalized Lotka-Volterra dynamics and
# returns the state variables over time 
# inputs: N0 = vector of initial population sizes (initial condition)
# r = vector of intrinsic growth rates
# K = vector of carrying capacities
# A = square interaction matrix
# times: sequence of times points to numerically integrate ODE
lotka_volterra <- function(N0, r, K, A, times = seq(0, 200, 0.01), 
                           formalism = "r", final_abundance = FALSE) {
  if (formalism == "r") {
    # list of parameters
    pars <- list(r = r, A = A)
    # function that returns rate of change
    model <- function(t, N0, pars) {
      dN_dt <- N0 * (pars$r + c(pars$A %*% N0))
      return(list(dN_dt))
    }
    # numerical integration ode45 Runge-Kutta method
    out <- ode(y = N0, times = times, func = model, parms = pars, method = "ode45")
  }
  if (formalism == "K") {
    pars <- list(r = r, K = K, A = A)
    model <- function(t, N0, pars) {
      dN_dt <- N0 * (pars$r / pars$K) * (pars$K + c(pars$A %*% N0))
      return(list(dN_dt))
    }
    out <- ode(y = N0, times = times, func = model, parms = pars, method = "ode45")
  }
  if (formalism == "K_typeII") {
    pars <- list(r = r, K = K, A = A)
    model <- function(t, N0, pars) {
      dN_dt <- pars$r * N0 * (1 + (1 / pars$K) * c(pars$A %*% diag(1 / (1 + N0)) %*% N0))
      return(list(dN_dt))
    }
    out <- ode(y = N0, times = times, func = model, parms = pars, method = "ode45")
  }
  if (formalism == "K_stochastic") {
    pars <- list(r = r, K = K, A = A)
    # defining deterministic part
    f <- function(u, p, t) {
      deterministic <- u * (p$r / p$K) * (p$K + c(p$A %*% u))
      return(deterministic)
    }
    # defining stochastic part
    g <- function(u, p, t) {
      s <- rep(1 / sqrt(length(p$K)), length(p$K))
      stochastic <- s * u * (p$r / p$K) * (p$K - c(p$A %*% u))
      return(stochastic)
    }
    # integration time steps
    time_step <- times[2] - times[1]
    # numerical integration
    sol <- sde.solve(f = f, g = g, u0 = N0, tspan = range(times), p = pars, saveat = time_step)
    out <- as.data.frame(cbind(sol$t, sol$u))
    names(out) <- paste("time", 1:length(pars$K))
  }
  out <- out[complete.cases(out),] # delete the NA cases (last row)
  if(final_abundance){
    out <- as.numeric(out[nrow(out), -1])
  }
  return(out)
}

# Function that solves the Lotka-Volterra dynamics and returns the 
# indexes of surviving species (i.e., abundance larger than threshold)
# N0: vector of initial population sizes
# r: vector of intrinsic growth rates
# K: vector of carrying capacities
# A: square interaction matrix
# times: sequence of times points to numerically integrate ODE
# extinct_tol: species with a final population size smaller than
# this value are considered extinct
# formalism: whether to use r or K Lotka-Volterra formalism
lv_pruning <- function(N0, r, K, A, times = seq(0, 200, 0.01), 
                       formalism = "r", extinct_tol = 0.000001) {
  # solve Lotka-Volterra dynamics
  eq <- lotka_volterra(N0, r, K, A, times, formalism, final_abundance = TRUE)
  # return pruned system
  which_surv_sp <- as.numeric(which(eq > extinct_tol))
  return(which_surv_sp)
}

# function that computes the **true** feasibility domain of a subset of species 
# in a large S(>2)-species community analytically
# only works for **globally stable** systems
# note that this is for just one region
# inputs: A = S by S interaction matrix of the entire community, diag(A) = -1 
#         species = index vector of the subset of species
# output: omega_comm_analytical = the normalized feasibility of the pair while others are transient
Omega_comm_analytical <- function(A, species){
  S <- nrow(A)
  if(is.null(species)){
    omega_comm_analytical <- Omega(-diag(S))
  }else{
    species <- sort(species)
    S <- nrow(A)
    B <- matrix(0, nrow = S, ncol = S)
    diag(B) <- 1
    
    other <- c(1:S)[-species]
    A_species <- cbind(A[,species], B[,other])
    omega_comm_analytical <- Omega(A_species)
  }
  omega_comm_analytical
}

# function that computes the **true overall** feasibility domain of a subset of species 
# in a large S(>2)-species community analytically
# only works for **globally stable** systems
# note that this is for a combination of regions
# inputs: A = S by S interaction matrix of the entire community, diag(A) = -1 
#         species = index vector of the subset of species
# output: omega_comm_analytical_all = the normalized feasibility of the pair across all possible comms
Omega_comm_analytical_all <- function(A, species){
  species <- sort(species)
  S <- nrow(A)
  S_sub <- length(species)
  S_other <- S - S_sub
  
  omega_comm_analytical_all <- Omega_comm_analytical(A, species)
  for(s in 1:S_other){
    other_comb <- combinations(S_other, s, c(1:S)[-species])
    n_other_comb <- nrow(other_comb)
    for(n in 1:n_other_comb){
      other <- other_comb[n,]
      comm <- sort(c(species, other))
      omega_comm_analytical_all <- omega_comm_analytical_all + Omega_comm_analytical(A, comm)
    }
  }
  omega_comm_analytical_all
}

# function that computes the probability of survival, augmentation, and exclusion
# on the *entire sphere* in a simultaneous assembly analytically
# inputs: A = S by S interaction matrix of the entire S-species community (invader + residents), diag(A) = -1 
#         invader = index vector of the invading species (only one)
# output: survival = probability of positive igr
#         persistence = probability of invader coexisting with residents
#         exclusion = probability of some (at least one) residents going extint
simultaneous_sph <- function(A, invader){
  survival <- Omega_comm_analytical_all(A, invader)
  persistence <- Omega(A)
  S <- nrow(A)
  res <- c(1:S)[-invader]
  exclusion <- survival - persistence
  list(survival = survival, persistence = persistence, exclusion = exclusion)
}

# function that computes the probability of colonization, augmentation, and replacement
# constrained in *the feasibility region of residents in isolation*
# in a sequential assembly using *brute force*
# inputs: A = S by S interaction matrix of the entire S-species community (invader + residents), diag(A) = -1 
#         invader = index vector of the invading species (only one)
# output: colonization = probability of positive igr
#         augmentation = probability of invader coexisting with residents
#         replacement = probability of some (at least one) residents going extint
sequential_iso <- function(A, invader, n_point = 5000){
  S <- nrow(A)
  res <- c(1:S)[-invader]
  # sample r
  r <- replicate(n_point, hypersphere_sampling(S), simplify = FALSE)
  spanA <- generate_span_vectors(A)
  B <- diag(S)
  # sample r from the geometric region of igr > 0 (assume global stability)
  check_igr <- mapply(inside_detection, vector = r,
                    MoreArgs = list(Span = cbind(B[, invader], spanA[, res])))
  res_only <- mapply(inside_detection, vector = r,
                         MoreArgs = list(Span = cbind(-B[, invader], spanA[, res])))
  r_inside <- r[sort(c(which(check_igr == 1), which(res_only == 1)))] # list of r inside the region
  n_point_inside <- length(r_inside) # number of r inside the region
  # calculate the colonization probability with positive igr 
  igr <- sum(check_igr)/n_point_inside
  N0 <- replicate(n_point_inside, initial_condition(S), simplify = FALSE) # sample initial conditions
  # check the surviving species within the feasibility region of residents in isolation
  surv_sp <- mapply(lv_pruning, r = r_inside, N0 = N0, MoreArgs = list(A = A), SIMPLIFY = FALSE)
  num_surv <- sapply(surv_sp, function(x) length(x))
  # augmentation probability
  augmentation <- sum(num_surv == S)/n_point_inside
  # replacement probability
  replacement <- (sum(num_surv < S) - sum(res_only))/n_point_inside
  # output
  list(colonization = igr, augmentation = augmentation, replacement = replacement)
}
