######################################################################################
# JAB model simulations                       ---------------------------------------
######################################################################################


# Library -----------------------------------------------------------------
library(deSolve)


# Different versions of the JAB model wrapped into functions ----------------------------------
## Here, we show three versions of the JAB model
## 1) the originally proposed iterateJAB()
## 2) an altered JAB model based on the reviewer's suggestions iterateJAB2()
## 3) for comparison an ODE version, calculateOdeJAB(), which can be numerically integrated with deSolve


## 1) Originally proposed JAB model
iterateJAB <- function(time, ## vector of times
                       state_0, ## vector of initial species states.
                       par){
  
  b <- par$b
  c_a <- par$c_a
  c_b <- par$c_b
  c_j <- par$c_j
  g <- par$g
  h <- par$h
  l <- par$l
  r <- par$r
  s <- par$s
  
  ba_a_avg <- 0.5
  ba_a_upper <- 0.7
  
  ## Set the count variables
  n <- length(r) # no. of species
  
  times_intern <- 1:max(time)
  n_times <- length(times_intern)
  
  # Prepare a state matrix
  whichstage <- rep(1:3, each = n)
  State <- matrix(rep(state_0, times = n_times), nrow = n_times, byrow = T)
  
  ## Here comes the model.
  for (t in 2:n_times) {
    
    ## States at t-1
    J <- State[t-1, whichstage == 1]
    A <- State[t-1, whichstage == 2]
    B <- State[t-1, whichstage == 3]
    
    ## The total basal area of big trees
    BA <- A * ba_a_avg + B
    BA_sum <- sum(BA)
    
    J_t <- (r*BA + l + J - g*J) / (1 + c_j*sum(J) + s*BA_sum) # count of juveniles J
    A_t <-  (g*J + A - h*A) / (1 + c_a*BA_sum) # count of small adults
    A_ba <- A * h * ba_a_upper # Basal area of small adults A. Conversion by multiplication with basal area of State exit (based on upper dhh boundary of the class)
    B_t <- (1+b)*(A_ba + B) / (1 + c_b*BA_sum)  # basal area of big adults B
    
    State[t, ] <- c(J_t, A_t, B_t)
  }
  
  whichtimes <- match(time, times_intern)
  State <- State[whichtimes,]
  State <- cbind(time = whichtimes, State)

  return(State)
}


## 2) Altered JAB model with limitation only acting on the state
iterateJAB2 <- function(time, ## vector of times
                        state_0, ## vector of initial species states.
                        par){
  
  b <- par$b
  c_a <- par$c_a
  c_b <- par$c_b
  c_j <- par$c_j
  g <- par$g
  h <- par$h
  l <- par$l
  r <- par$r
  s <- par$s
  
  ba_a_avg <- 0.5
  ba_a_upper <- 0.7
  
  ## Set the count variables
  n <- length(r) # no. of species
  
  times_intern <- 1:max(time)
  n_times <- length(times_intern)
  
  # Prepare a state matrix
  whichstage <- rep(1:3, each = n)
  State <- matrix(rep(state_0, times = n_times), nrow = n_times, byrow = T)
  
  ## Here comes the model.
  for (t in 2:n_times) {
    
    ## States at t-1
    J <- State[t-1, whichstage == 1]
    A <- State[t-1, whichstage == 2]
    B <- State[t-1, whichstage == 3]
    
    ## The total basal area of big trees
    BA <- A * ba_a_avg + B
    BA_sum <- sum(BA)
    
    lim_J <- (1 + c_j*sum(J) + s*BA_sum)
    J_t <- l + r*BA + (J - g*J)/lim_J
    
    lim_A <- (1 + c_a*BA_sum)
    A_t <-  g*J/lim_J + (A - h*A)/lim_A
    
    A_ba <- ba_a_upper * (A * h)/lim_A
    B_t <- A_ba + (1+b)*B/(1 + c_b*BA_sum)
    
    State[t, ] <- c(J_t, A_t, B_t)
  }
  
  whichtimes <- match(time, times_intern)
  State <- State[whichtimes,]
  State <- cbind(time = whichtimes, State)
  
  return(State)
}


## 3. ODE version of the JAB model for comparison
calculateOdeJAB <- function(time,
                             state, # vector of species states.
                             par){
  
  b <- par$b
  c_a <- par$c_a
  c_b <- par$c_b
  c_j <- par$c_j
  g <- par$g
  h <- par$h
  l <- par$l
  r <- par$r
  s <- par$s
  
  ba_a_avg <- 0.5 ## constants for basal area conversion
  ba_a_upper <- 0.7

  J <- state[1:2]
  A <- state[3:4]
  B <- state[5:6]
  BA <- (A * ba_a_avg) + B
  BA_sum <- sum(BA)

  dJ <- l + r * BA - (c_j*sum(J) + s*BA_sum + g)*J
  dA <- g * J - (c_a*BA_sum + h)*A
  dB <- A * h * ba_a_upper + B*b - (c_b*BA_sum)*B
  
  return(list(c(dJ, dA, dB)))
}


# Set up simulation -------------------------------------------------------------------
set.seed(1)
time <- seq(1, 80, 4)

## Equal parameters among species, loosely based on the estimates
par <- list(b = c(-4, -4),
            c_a = c(-5, -5),
            c_b = c(-7, -7),
            c_j = c(-7, -7),
            g = c(-5, -5),
            h = c(-3, -3),
            l = c(5, 5),
            r = c(4, 4),
            s = c(-5, -5))

## Equal states, loosely based on initial states in the data (change sd in rnorm for some randomness)
state_0 <- exp(rnorm(6, rep(c(8, 5, 2), each = 2), 0)) # state_0 <- rep(c(4100, 236, 15), each = 2)

## Wrapper for  plotting trajectories given a Matrix with observations as rows and columns 1: time, 2: j, 3: J, 4: a, 5: A, etc.
wrapMatplot <- function(Mat) {
  matplot(Mat[, 1], Mat[, -1],
          type = "p", pch = c("j", "J", "a", "A", "b", "B"), col = rep(c(1, 3), 3),
          xlab = "time", ylab = "abundance",
          log = "y")
}




# Simulations ------------------------------------------------------

## Equal states and equal parameters between species --------------

Sim_JAB <- iterateJAB(time, state_0, lapply(par, exp))
wrapMatplot(Sim_JAB)
## Here, letter case and color indicate species, while J, A, B indicate stage

## The ODE version has very similar dynamics
Sim_JAB_ode <- deSolve::ode(state_0, time, calculateOdeJAB, lapply(par, exp))
wrapMatplot(Sim_JAB_ode) 

## We propose JAB2 which addresses the reviewer's suggestions
## - better parameter interpretability through limitation acting only between steps
Sim_JAB2 <- iterateJAB2(time, state_0, lapply(par, exp))
wrapMatplot(Sim_JAB2)


## Different recruitment rates clearly alter equilibria -------------------------
## regardless of model version

par_r <- within(par, { r <- c(3.5, 4)}) ## change only r within par

Sim_JAB_r <- iterateJAB(time, state_0, lapply(par_r, exp))
wrapMatplot(Sim_JAB_r)

Sim_odeJAB_r <- deSolve::ode(state_0, time, calculateOdeJAB, lapply(par_r, exp))
wrapMatplot(Sim_odeJAB_r)

Sim_JAB2_r <- iterateJAB2(time, state_0, lapply(par_r, exp))
wrapMatplot(Sim_JAB2_r)


par_l <- within(par, { l <- c(6, 4)}) ## change l; Note that l is acting only linearly, therefore larger assumed differences for faster equilibria

Sim_JAB_l <- iterateJAB(seq(1, 1500, 100), state_0, lapply(par_l, exp))
wrapMatplot(Sim_JAB_l)

Sim_odeJAB_l <- deSolve::ode(state_0, seq(1, 1500, 100), calculateOdeJAB, lapply(par_l, exp))
wrapMatplot(Sim_odeJAB_l)

Sim_JAB2_l <- iterateJAB2(seq(1, 1500, 100), state_0, lapply(par_l, exp))
wrapMatplot(Sim_JAB2_l)



