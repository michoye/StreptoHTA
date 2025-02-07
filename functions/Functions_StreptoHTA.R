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
calculate_ce_out <- function(l_params_all, n_wtp = 5000, verbose = FALSE){ # User defined
  with(as.list(l_params_all), {
    set.seed(1) 

# ## General setup 
    # Calculate model input parameters
    v_names_cycles  = paste("cycle", 0:n_cycles)    # cycle names
    n_states        = length(v_names_states)        # number of health states 
    n_str           = length(v_names_str)           # number of strategies
    
    ## 03.2 Calculate internal model parameters
    
    # Monthly probabilities of death
    # load age dependent probability
    p_mort   <- read.csv("data/mortProb_age.csv")
    # Adjust to monthly probabilities
    p_mort$p_HD <- 1 - (1 - p_mort$p_HD)^(1/12)
    
    df_p_HS     = data.frame(Sex = c("Child", "Adult"), p_HS = c(p_HS_child, p_HS_adult))
    
    # probability to die in S by cycle (is increasing)
    p_SD <- 1 - (1 - c(0.005, 0.007, 0.01, 0.02, 0.03, rep(0.04, n_cycles - 5)))^(1/12)
    
    
    ### Discount weight for costs and effects 
    v_dwc     <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
    v_dwe     <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
    
    # Within-cycle correction (WCC) - method  options Simpson's 1/3 rule, "half-cycle" or "none" 
    v_wcc    <- darthtools::gen_wcc(n_cycles = n_cycles, 
                                    method = "Simpson1/3") # vector of wcc
    
    ## 04.1 Static characteristics
    # randomly sample the sex of an individual (50% female)
    
    # load age distribution
    dist_Age <- read.csv("data/MyPopulation-AgeDistribution.csv") 
    # sample from the age distribution the initial age for every individual
    v_age  <- sample(x = dist_Age$age, prob = dist_Age$prop, size = n_i, replace = TRUE)
    
    
    ## 04.2 Dynamic characteristics 
    # Specify the initial health state of the individuals 
    # everyone begins in the healthy state (in this example)
    v_M_init  <- rep("H", times = n_i)   
    v_Ts_init <- rep(0, n_i)  # a vector with the time of being sick at the start of the model 
    # sample the treatment effect modifier at baseline
    v_x     <- runif(n_i, min = 0.95, max = 1.05)

    ## 04.3 Create a dataframe with the individual characteristics 
    # create a data frame with each individual's 
    # ID number, treatment effect modifier, age and initial time in sick state 
    df_X  <- data.frame(ID = 1:n_i, M_x = v_x, Age = v_age, n_cycles_s = v_Ts_init, M_init = v_M_init)


    # 05 Define Simulation Functions
    
    ## 05.1 Probability function
    #The `Probs` function updates the transition probabilities of every cycle is shown below 
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


## 05.2 Cost function
# The `Costs` function estimates the costs at every cycle.
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



## 05.3 Health outcome function
#The `Effs` function to update the utilities at every cycle.
Effs <- function (M_t, cycle_length = 1, Diag = "SoC") {
  # Arguments:
    # M_t: health state occupied at cycle t (character variable)
    # cl:  cycle length (default is 1)
    # df_X: data frame with individual characteristics data (TA added)
    # Diag:  Diagnosis
  # Returns:
    # QALYs accrued this cycle

  q_t <- 0                          # by default the utility for everyone is zero
  q_t[M_t == "H"]    <- u_H         # update the utility if healthy

  if (Diag == "SoC") {  # update the utility if sick under standard of care
    q_t[M_t == "S"] <- u_S
  } else if (Diag == "AB") {
    q_t[M_t == "S"] <- u_diagAB * df_X$x[M_t == "S"]
  }

  q_t[M_t == "D"]    <- u_D         # update the utility if dead

  QALYs <- q_t * cycle_length  # calculate the QALYs during cycle t
  return(QALYs)      # return the QALYs accrued this cycle
}


# Microsimulation
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
  m_E[, 1] <- Effs( m_M[, 1], cycle_length = 1)   # QALYs per individual during cycle 0
  
  # open a loop for time running cycles 1 to n_cycles 
  for (t in 1:n_cycles) {
    # calculate the transition probabilities for the cycle based on health state t
    m_P <- Probs(m_M[, t], df_X, Diag = Diag)
    # # check if transition probabilities are between 0 and 1
    # check_transition_probability(m_P, verbose = TRUE)
    # # check if each of the rows of the transition probabilities matrix sum to one
    # check_sum_of_transition_array(m_P, n_rows = n_i, n_cycles = n_cycles, verbose = TRUE)
    
    # sample the next health state and store that state in matrix m_M
    m_M[, t + 1]  <- samplev(m_P, 1)    
    # calculate costs per individual during cycle t + 1
    m_C[, t + 1]  <- Costs(m_M[, t + 1], Diag = Diag)  
    # calculate QALYs per individual during cycle t + 1
    m_E[, t + 1]  <- Effs (m_M[, t + 1], cycle_length = 1)  
    
    # update time since illness onset for t + 1 
    # NOTE: this code has a "reset of history" for time being sick
    # once someone is not "Sick" anymore, we reset n_cycles (set back to zero)
    # when you don't want a "reset" replace the last zero with the current value
    df_X$n_cycles <- if_else(m_M[, t + 1] == "S", df_X$n_cycles + 1, 0) 
    
    # update the age of individuals that are alive
    df_X$Age[m_M[, t + 1] != "D"]  <- df_X$Age[m_M[, t + 1] != "D"] + 1
    
    # # Display simulation progress
    # if (t %% 5 == 0 | t == n_cycles) {  # Print progress every 5 cycles and at the last cycle
    #   cat(sprintf("\rSimulation progress: %d/%d cycles (%.1f%% complete)", 
    #               t, n_cycles, (t / n_cycles) * 100))
    #   flush.console()  # Ensures immediate printing in RStudio
    # }
    
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
  
} # end of the `MicroSim` function 

# Run the simulation 
outcomes_SoC  <- MicroSim(n_i = n_i, df_X = df_X, seed = 1, cycle_length = cycle_length, Diag = "SoC")  # Run for Standard of Care
outcomes_diagAB <- MicroSim(n_i = n_i, df_X = df_X, seed = 1, cycle_length = cycle_length, Diag = "AB") # Run simulation for strategy AB

# Cost-Effectiveness Analysis
# store the mean costs of each strategy in a new variable C (vector of costs)
v_C <- c(outcomes_SoC$tc_hat, outcomes_diagAB$tc_hat)
# store the mean QALYs of each strategy in a new variable E (vector of effects)
v_E <- c(outcomes_SoC$te_hat, outcomes_diagAB$te_hat)

## Vector with discounted net monetary benefits (NMB)
v_nmb_d <- v_E * n_wtp - v_C

# Dataframe with discounted costs, effectiveness and NMB
df_ce <- data.frame(Strategy = v_names_str,
                    Cost     = v_C,
                    Effect   = v_E,
                    NMB      = v_nmb_d)
return(df_ce)
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
generate_psa_params <- function(n_sim = 100, seed = 071818){
  set.seed(seed) # set a seed to be able to reproduce the same results

  df_psa <- data.frame(
    # Transition probabilities
    # probability of becoming sick when healthy
    p_HS         = 1 - (1 - (rlnorm(n_sim,meanlog =   log(0.1), sdlog = 0.1)))^(1/12), 
    # probability of recovering to healthy when sick and diagnosed through AB
    p_SH_SoC  = 1 - (1 - (pmin(rlnorm(n_sim, meanlog = log(0.88), sdlog = 0.05), 0.99)))^(1/12) ,      # standard of care
    p_SH_diagAB = 1 - (1 - (pmin(rlnorm(n_sim,meanlog =   log(0.93), sdlog = 0.04), 0.99)))^(1/12) ,
    
    ## State rewards
    # Costs
    c_H       = pmin(pmax(rgamma(n_sim, shape = 100, scale = 5), 200), 800) / 12,        # cost of one cycle in healthy state
    c_S       = pmin(pmax(rgamma(n_sim, shape = 100, scale = 5), 400), 1100) / 12,       # cost of one cycle in sick state
    c_D       = 0,                                            # cost of one cycle in dead state
    c_diagSoC = pmin(pmax(rgamma(n_sim, shape = 100, scale = 5), 10), 90) / 12, # monthly (very inflated) cost (for this assignment) of receiving diagnosis SoC + treatment when in Sick
    c_diagAB = pmin(pmax(rgamma(n_sim, shape = 100, scale = 5), 20), 110) / 12, # monthly (very inflated) cost (for this assignment) of receiving diagnosis AB + treatment when in Sick

    # Utilities
    u_H       = rbeta(n_sim, shape1 =  1.5, shape2 = 0.0015), # utility when healthy 
    u_S       = pmin(pmax(rbeta(n_sim, shape1 = 172.55, shape2 = 30.45), 0.8), 0.9)
,   # utility when sick
    u_D       = 0,                                             # utility when dead
    u_diagAB = pmin(pmax(rbeta(n_sim, shape1 = 172.55, shape2 = 30.45), 0.85), 0.95)  # monthly utility when receiving diagnosis AB when in Sick
  )
  return(df_psa)
}


