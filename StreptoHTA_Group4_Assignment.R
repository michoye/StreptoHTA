#------------------------------------------------------------------------------#
####              Calculate cost-effectiveness outcomes                     ####
#------------------------------------------------------------------------------#
#' Calculate cost-effectiveness outcomes
#'
#' \code{calculate_ce_out} calculates costs and effects for a given vector of parameters using a simulation model.
#' @param l_params_all List with all parameters of decision model
#' @param n_wtp Willingness-to-pay threshold to compute net monetary benefits (
#' NMB)
#' @return A dataframe with discounted costs, effectiveness and NMB.
#' @export
calculate_ce_out <- function(l_params_all, n_wtp = 10000, verbose = FALSE){ # User defined
  with(as.list(l_params_all), {


## General setup 
set.seed(1)               # set the seed  
cycle_length   <- 1/12       # cycle length equal to one year (use 1/12 for monthly)
n_cycles       <- 60      # time horizon, number of cycles
n_i            <- 100000  # number of individuals

# the 4 health states of the model:
v_names_states  <- c("H",  # Healthy (H)
                     "S", # Sick (S)
                     "D")  # Dead (D)
v_names_cycles  <- paste("cycle", 0:n_cycles)    # cycle names
n_states        <- length(v_names_states) # number of health states               

### Discounting factors 
d_c <- 0.03 # monthly discount rate for costs 
d_e <- 0.03 # monthly discount rate for QALYs


### Strategies 
v_names_str   <- c("Standard of care",   # This is the RADT
                   "Strategy AB")        # Throat Culture
n_str         <- length(v_names_str)     # number of strategies

### Transition probabilities 
# (all non-probabilities are conditional on survival)

# if (df_X[, dist_Age$age] > 18) {
#  p_HS <- 0.1
# } else {
#  p_HS <- 0.4
# }

p_HS         <- 1 - (1 - 0.1)^(1/12)     # probability of becoming sick when healthy
p_SH_SoC         <- 1 - (1 - 0.88)^(1/12)      # probability of recovering to healthy when sick and diagnosed through SoC
p_SH_diagAB          <- 1 - (1 - 0.93)^(1/12)  # probability of recovering to healthy when sick and diagnosed through AB


# Monthly probabilities of death
# load age dependent probability
p_mort   <- read.csv("data/mortProb_age.csv")
# Adjust to monthly probabilities
p_mort$p_HD <- 1 - (1 - p_mort$p_HD)^(1/12)

# load age distribution
dist_Age <- read.csv("data/MyPopulation-AgeDistribution.csv") 

# Monthly probabilities of becoming sick
# We created a new matrix
df_p_sick <- data.frame(Age = 1:107)
df_p_sick$p_HS <- ifelse(df_p_sick$Age <= 18, 0.4, 0.1)

# probability to die in S by cycle (is increasing)
v_p_SD <- 1 - (1 - c(0.005, 0.007, 0.01, 0.02, 0.03, rep(0.04, n_cycles - 5)))^(1/12)


### State rewards 
#### Costs 
c_H     <- 500 / 12  # monthly cost of being Healthy
c_S     <- 750 / 12  # monthly cost of being Sick
c_D     <- 0     # monthly cost of being dead
c_diagSoC <- 30 # monthly (very inflated) cost (for this assignment) of receiving diagnosis SoC + treatment when in Sick
c_diagAB <- 70 # monthly (very inflated) cost (for this assignment) of receiving diagnosis AB + treatment when in Sick

#### Utilities 
u_H     <- 1     # monthly utility of being Healthy
u_S     <- 0.85  # monthly utility of being Sick
u_D     <- 0     # monthly utility of being dead
u_diagAB <- 0.95  # monthly utility when receiving diagnosis AB when in Sick


## ------------------------------------------------------------------------------------------------------------------------------------------------------
m_P_diag <- matrix(0, nrow = n_states, ncol = n_states, dimnames = list(v_names_states, v_names_states))
m_P_diag["H", "H" ] = ""
m_P_diag["H", "S" ] = ""
m_P_diag["H", "D" ] = ""

m_P_diag["S", "H"] = ""
m_P_diag["S", "S" ] = ""
m_P_diag["S", "D" ] = ""

m_P_diag["D", "D" ] = ""

layout.fig <- c(2, 1)

# Save plot as PNG
png("figs/transition_diagram.png", width = 800, height = 600)
# Generate the plot
plotmat(t(m_P_diag), t(layout.fig), self.cex = 0.5, curve = 0.04, arr.pos = 0.7,
        latex = T, arr.type = "curved", relsize = 0.85, box.prop = 0.8,
        cex = 0.8, box.cex = 0.7, lwd = 1)
# Close the graphical device
dev.off()


## ------------------------------------------------------------------------------------------------------------------------------------------------------
### Discount weight for costs and effects 
v_dwc <- (1 / (1 + d_c)^(1/12)) ^ (0:n_cycles)
v_dwe <- (1 / (1 + d_e)^(1/12)) ^ (0:n_cycles)

# Within-cycle correction (WCC) - method  options Simpson's 1/3 rule, "half-cycle" or "none" 
v_wcc    <- darthtools::gen_wcc(n_cycles = n_cycles, 
                                method = "Simpson1/3") # vector of wcc


## ------------------------------------------------------------------------------------------------------------------------------------------------------
# sample the treatment effect modifier at baseline
v_x     <- runif(n_i, min = 0.95, max = 1.05) 
# sample from the age distribution the initial age for every individual
v_age0  <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE) 


## ------------------------------------------------------------------------------------------------------------------------------------------------------
v_M_init  <- rep("H", times = n_i)   # Specify the initial health state of the individuals 
v_n_cycles_s_init <- rep(0, n_i)  # everyone begins in the healthy state (in this example)


## ------------------------------------------------------------------------------------------------------------------------------------------------------
# data frame with each individual's 
# ID number, treatment effect modifier, age and initial time in sick state and initial health state at the start of the simulation
df_X  <- data.frame(ID = 1:n_i, M_x = v_x, Age = v_age0, n_cycles = v_n_cycles_s_init, M_init = v_M_init)

head(df_X)  # print the first rows of the dataframe

# See also:
# Define the plot
p <- ggplot(df_X, aes(x = Age)) +
    geom_histogram()
# Save the plot as PNG
ggsave("figs/histogram_age.png", plot = p, width = 8, height = 6, dpi = 300)


## ------------------------------------------------------------------------------------------------------------------------------------------------------
Probs <- function(M_t, df_X, Diag = "SoC") { 
  # Arguments:
    # M_t:  health state occupied at cycle t (character variable)
    # df_X: data frame with individual characteristics data 
    # Diag:  diagnosis
  # Returns: 
    # transition probabilities for that cycle

# Diagnosis specific transition probabilities
  if (Diag == "SoC") {
    p_SH <- p_SH_SoC
  } else if (Diag == "AB") {
    p_SH <- p_SH_diagAB 
  }
   
  # create matrix of state transition probabilities
  m_p_t           <- matrix(data = 0, nrow = n_states, ncol = n_i)  
  # give the state names to the rows
  rownames(m_p_t) <-  v_names_states                               
  
  # lookup baseline probability and rate of dying based on individual characteristics
  p_HD_all <- inner_join(x = df_X, y = p_mort, by = c("Age"))
  p_HD     <- p_HD_all[M_t == "H", "p_HD"]
  
  # lookup baseline probability and rate of becoming sick based on individual characteristics (Age)  
  p_HS_all <- inner_join(x = df_X, y = df_p_sick, by = c("Age"))
  p_HS     <- p_HS_all[M_t == "H", "p_HS"]
 
  # update m_p_t with the appropriate probabilities 
  # (all non-death probabilities are conditional on survival)
  # transition probabilities when Healthy 
  m_p_t["H", M_t == "H"] <- (1 - p_HD) * (1 - p_HS)
  m_p_t["S", M_t == "H"] <- (1 - p_HD) * p_HS 
  m_p_t["D", M_t == "H"] <-      p_HD 
  
  
  # transition probabilities when Sick 
  m_p_t["H", M_t == "S"] <-  (1 - v_p_SD[df_X$n_cycles]) * p_SH
  m_p_t["S", M_t == "S"] <-  (1 - v_p_SD[df_X$n_cycles]) * (1 - p_SH)
  m_p_t["D", M_t == "S"] <-      v_p_SD[df_X$n_cycles]   
  

  # transition probabilities when Dead
  m_p_t["H", M_t == "D"] <- 0
  m_p_t["S", M_t == "D"] <- 0
  m_p_t["D", M_t == "D"] <- 1  
  
  return(t(m_p_t))
}       


## ------------------------------------------------------------------------------------------------------------------------------------------------------
Costs <- function (M_t, Diag = "SoC") {
  # Arguments:
    # M_t: health state occupied at cycle t (character variable)
    # Diag:  diagnosis
  # Returns: 
    # costs accrued in this cycle

    # Diagnosis specific costs
  if (Diag == "SoC") {
    c_diag <- c_diagSoC
  } else if (Diag == "AB") {
    c_diag <- c_diagAB
  } 
  
  c_t <- c()
  c_t[M_t == "H"] <- c_H  # costs accrued by being healthy this cycle
  c_t[M_t == "S"] <- c_S + c_diag         # costs accrued by being sick this cycle
  c_t[M_t == "D"] <- c_D          # costs at dead state
  
  return(c_t)  # return costs accrued this cycle
}



## ------------------------------------------------------------------------------------------------------------------------------------------------------
Effs <- function (M_t, cycle_length = 1, Diag = "SoC") {
  # Arguments:
    # M_t: health state occupied at cycle t (character variable)
    # cl:  cycle length (default is 1)
    # df_X: data frame with individual characteristics data (TA added)
    # Diag:  Diagnosis
  # Returns: 
    # QALYs accrued this cycle 

  u_t <- 0                          # by default the utility for everyone is zero
  u_t[M_t == "H"]    <- u_H         # update the utility if healthy
  
  if (Diag == "SoC") {  # update the utility if sick under standard of care
    u_t[M_t == "S"] <- u_S  
  } else if (Diag == "AB") {
    u_t[M_t == "S"] <- u_diagAB * df_X$x[M_t == "S"]  
  }
  
  u_t[M_t == "D"]    <- u_D         # update the utility if dead
  
  QALYs <- u_t * cycle_length  # calculate the QALYs during cycle t
  return(QALYs)      # return the QALYs accrued this cycle
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------
MicroSim <- function(n_i, df_X, Diag = "SoC", seed = 1, cycle_length = 1, verbose ) {
  # Arguments:  
    # n_i: number of individuals
    # df_X: data frame with individual data 
    # Diag: diagnosis
    # seed: seed for the random number generator, default is 1
    # cycle_length: cycle length 
  # Returns:
    # results: data frame with total cost and QALYs
  
  set.seed(seed) # set a seed to be able to reproduce the same results
  
  # Here, we create three matrices called m_M, m_C and m_E
  # number of rows is equal to the n_i, the number of columns is equal to n_cycles 
  # (the initial state and all the n_cycles cycles)
  # m_M is used to store the health state information over time for every individual
  # m_C is used to store the costs information over time for every individual
  # m_E is used to store the effects information over time for every individual
  
  m_M <- m_C <- m_E <-  matrix(nrow = n_i, ncol = n_cycles + 1, 
                               dimnames = list(paste("ind"  , 1:n_i, sep = " "), 
                                               paste("cycle", 0:n_cycles, sep = " ")))  
 
  m_M[, 1] <- as.character(df_X$M_init) # initial health state at cycle 0 for individual i
  m_C[, 1] <- Costs(m_M[, 1])           # costs per individual during cycle 0
  m_E[, 1] <- Effs( m_M[, 1], cycle_length = cycle_length)   # QALYs per individual during cycle 0
  
  # open a loop for time running cycles 1 to n_cycles 
  for (t in 1:n_cycles) {
    print(t)
    # calculate the transition probabilities for the cycle based on health state t
    m_P <- Probs(m_M[, t], df_X, Diag = Diag)
    # check if transition probabilities are between 0 and 1
    check_transition_probability(m_P, verbose = TRUE)
    # check if each of the rows of the transition probabilities matrix sum to one
    check_sum_of_transition_array(m_P, n_rows = n_i, n_cycles = n_cycles, verbose = TRUE)
    
    # sample the next health state and store that state in matrix m_M
    m_M[, t + 1]  <- samplev(m_P, 1)    
    # calculate costs per individual during cycle t + 1
    m_C[, t + 1]  <- Costs(m_M[, t + 1], Diag = Diag)  
    # calculate QALYs per individual during cycle t + 1
    m_E[, t + 1]  <- Effs (m_M[, t + 1], cycle_length = cycle_length)  
    
    # update time since illness onset for t + 1 
    # NOTE: this code has a "reset of history" for time being sick
    # once someone is not "Sick" anymore, we reset n_cycles (set back to zero)
    # when you don't want a "reset" replace the last zero with the current value
    df_X$n_cycles <- if_else(m_M[, t + 1] == "S", df_X$n_cycles + 1, 0) 
    # update the age of individuals that are alive
    df_X$Age[m_M[, t + 1] != "D"]  <- df_X$Age[m_M[, t + 1] != "D"] + 1
    
    # Display simulation progress
    if(t/(n_cycles/10) == round(t/(n_cycles/10), 0)) { # display progress every 10%
      cat('\r', paste(t/n_cycles * 100, "% done", sep = " "))
    }
    
  } # close the loop for the time points 
  
  # Discounted total expected QALYs and Costs per strategy and apply cycle correction
  tc      <- m_C %*% (v_dwc * v_wcc)  # total (discounted and cycle corrected) cost per individual
  te      <- m_E %*% (v_dwe * v_wcc)  # total (discounted and cycle corrected) QALYs per individual 
  
  tc_hat  <- mean(tc)       # average (discounted and cycle corrected) cost 
  te_hat  <- mean(te)       # average (discounted and cycle corrected) QALY  
  # store the results from the simulation in a list
  results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc , te = te, 
                  tc_hat = tc_hat, te_hat = te_hat)   
  
  return(results)  # return the results

}


outcomes_SoC  <- MicroSim(n_i = n_i, df_X = df_X, seed = 1, cycle_length = cycle_length, Diag = "SoC")  # Run for Standard of Care
outcomes_diagAB <- MicroSim(n_i = n_i, df_X = df_X, seed = 1, cycle_length = cycle_length, Diag = "AB") # Run simulation for strategy AB


# store the mean costs of each strategy in a new variable C (vector of costs)
v_C <- c(outcomes_SoC$tc_hat, outcomes_diagAB$tc_hat)
# store the mean QALYs of each strategy in a new variable E (vector of effects)
v_E <- c(outcomes_SoC$te_hat, outcomes_diagAB$te_hat)

# use dampack to calculate the ICER
df_cea <- calculate_icers(cost       = v_C,
                          effect     = v_E,
                          strategies = v_names_str)
return(df_cea)

  }
  )
}


#------------------------------------------------------------------------------#
####             Generate a PSA input parameter dataset                     ####
#------------------------------------------------------------------------------#
#' Generate parameter sets for the probabilistic sensitivity analysis (PSA)
#'
#' \code{generate_psa_params} generates a PSA dataset of the parameters of the 
#' cost-effectiveness analysis.
#' @param n_sim Number of parameter sets for the PSA dataset
#' @param seed Seed for the random number generation
#' @return A data.frame with a PSA dataset of the parameters of the 
#' cost-effectiveness analysis
#' @export
generate_psa_params <- function(n_sim = 1000, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results
  df_psa <- data.frame(
    # Transition probabilities
    # p_HD is defined in function calculate_ce_out
    # probability of becoming sick when healthy
    p_HS         <- 1 - (1 - (rlnorm(n_sim,meanlog =   log(0.1), sdlog = 0.1)))^(1/12), 
    # probability of recovering to healthy when sick and diagnosed through AB
    p_SH_SoC  = 1 - (1 - (pmin(rlnorm(n_sim, meanlog = log(0.88), sdlog = 0.05), 0.99)))^(1/12) ,      # standard of care
    p_SH_diagAB = 1 - (1 - (pmin(rlnorm(n_sim,meanlog =   log(0.93), sdlog = 0.04), 0.99)))^(1/12) ,
    
    ## State rewards
    # Costs
    c_H       = pmin(pmax(rgamma(1000, shape = 100, scale = 5), 400), 600) / 12,        # cost of one cycle in healthy state
    c_S       = pmin(pmax(rgamma(1000, shape = 100, scale = 5), 600), 900) / 12,       # cost of one cycle in sick state
    c_D       = 0,                                            # cost of one cycle in dead state
    c_diagSoC <- 30, # monthly (very inflated) cost (for this assignment) of receiving diagnosis SoC + treatment when in Sick
    c_diagAB <- 70, # monthly (very inflated) cost (for this assignment) of receiving diagnosis AB + treatment when in Sick
    
    # Utilities
    u_H       = rbeta(n_sim, shape1 =  1.5, shape2 = 0.0015), # utility when healthy 
    u_S       = pmin(pmax(rbeta(n_sim, shape1 = 172.55, shape2 = 30.45), 0.8), 0.9)
,   # utility when sick
    u_D       = 0,                                             # utility when dead
    u_diagAB <- pmin(pmax(rbeta(n_sim, shape1 = 172.55, shape2 = 30.45), 0.85), 0.95)  # monthly utility when receiving diagnosis AB when in Sick
  )
  return(df_psa)
}


