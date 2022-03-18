library(purrr)
library(tidyverse)
## list of parameters to loop over for sensitivity
params <- c("cond_sp1_surv","cond_sp2_surv", "cond_sp2_col", "cond_sp1_col", "sp1_surv", "sp2_surv","sp1_col", "sp2_col")

## list of three occupancy states
#sensitivity_metric <- c("both", "sp1", "sp2")

## range of parameter values ## make fewer than equil bc takes a long time to run
param_vals <- seq(0,0.9, length.out = 25) 

#number of sites

nsites = 104

## parameter values for parameters that are not changing
## loop over a set of these values  # for each base value there are 4 initial conditions so don't want too many...
base_val = c(0.1,0.5,0.7)

##initial occupancy conditions

init_cond <- matrix(c(0.25,0.25,0.25,0.25, 
                      0.3,.3, 0.4, 0, 
                      0.5,0.1,0.3,0.1,
                      0.7, 0.2, 0.06, 0.04 ),
                    nrow = 4, ncol = 4, byrow = TRUE)


## try with sp1 and sp2 reversed and make sure get the same results

#init_cond <- matrix(c(0.25,0.25,0.25,0.25, 
#                     0.3,.4, 0.3, 0, 
#                    0.5,0.3,0.1,0.1,
#                   0.7, 0.06, 0.2, 0.04 ),
#                nrow = 4, ncol = 4, byrow = TRUE)



## create transition matrix
sites <- c("Empty", "Sp1_only", "Sp2_only", "Both")

trans <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(sites, sites) )

## number of times to repeat calculating occupancy for same param set (loop d)

repititions = seq(1,100)

## need array to hold reps of calculated occupancy
occup <- c("empty","both" ,"sp1" , "sp2" )

rep_array <- array(data = NA, dim = c(length(param_vals)-1, 4, length(repititions), length(param_vals), 4),
                   dimnames = list(NULL,occup , NULL, NULL, occup))


years = 500

## storing transition states

state_trans <- array(data = 0, dim = c(4, years, nsites, length(params), length(param_vals), nrow(init_cond), length(base_val)),
                     dimnames = list(sites, NULL, NULL, params,NULL, NULL, base_val))


## functions

sum_fun <- function(x) {
  x <- as.data.frame(x)
  names(x) <- "type"
  x$type <- factor(x$type, levels = c(
    "Empty", "Sp1_only", "Sp2_only", "Both"
  ))
  x <- x %>% count(type, .drop = FALSE) %>% as.data.frame() %>% dplyr::select(n) %>% c()
  x$n
}


## criteria for equilibrium

eps.all <- 0.02
eps.one <- 0.02

system.time(
  ## loop over base vals
  for(q in 1:length(base_val)){
    
    ## loop is over initial conditions 
    for(i in 1:nrow(init_cond)) { 
      
      ## set up initial condition vector
      past_occ <- c(rep("empty", round(init_cond[i,1] * nsites)), 
                    rep("sp1", round(init_cond[i,2] * nsites) ), 
                    rep("sp2",  round(init_cond[i,3] * nsites)),
                    rep("both",  round(init_cond[i,4] * nsites)))  
      
      
      if(length(past_occ) < nsites){
        past_occ <- c(past_occ, rep("empty", nsites-length(past_occ)))
      } else if (length(past_occ > nsites)){
        past_occ <- sample(past_occ, nsites)
      } else if (length(past_occ = nsites)){
        past_occ <- past_occ
      }
      
      
      
      
      
      ## the next loop is over the 8 different parameters  
      for(j in 1:length(params)){
        
        ## all the parameters are set to the base value that was set outside of the loops
        param_base <- c(cond_sp1_surv = base_val[q],
                        cond_sp2_surv = base_val[q],
                        cond_sp2_col = base_val[q],
                        cond_sp1_col = base_val[q],
                        sp1_surv = base_val[q],
                        sp2_surv = base_val[q],
                        sp1_col = base_val[q],
                        sp2_col = base_val[q])
        
        
        ## set up occupancy matrix for new occupancy 
        occ_sp <- matrix(data = NA, nrow = nsites, ncol = length(param_vals) )
        
        
        
        ## the next loop calculates occupancy in the next time step and 
        ## loops over all parameter values (set as param_vals outside of the loop)
        ## for a single parameter
        
        for(k in 1:length(param_vals)){
          
          ## set parameter of interest (parameter j) to the appropriate value
          
          param_base[params[j]] <- param_vals[k] 
          
          # fill transition matrix
          trans[1,1] <- 1- (param_base["sp1_col"] * (1- param_base["sp2_col"])) - 
            (param_base["sp2_col"] * (1- param_base["sp1_col"])) - 
            (param_base["sp2_col"] * param_base["sp1_col"])
          trans[2,1] <- param_base["sp1_col"] * (1- param_base["sp2_col"])
          trans[3,1] <- param_base["sp2_col"] * (1- param_base["sp1_col"])
          trans[4,1] <- param_base["sp2_col"] * param_base["sp1_col"]
          
          trans[1,2] <- 1 - (param_base["sp1_surv"] * (1-param_base["cond_sp2_col"])) - 
            ((1-param_base["sp1_surv"]) * param_base["cond_sp2_col"]) - 
            (param_base["sp1_surv"] * param_base["cond_sp2_col"])
          trans[2,2] <- param_base["sp1_surv"] * (1-param_base["cond_sp2_col"])
          trans[3,2] <- (1-param_base["sp1_surv"]) * param_base["cond_sp2_col"]
          trans[4,2] <- param_base["sp1_surv"] * param_base["cond_sp2_col"]
          
          trans[1,3] <- 1-((1-param_base["sp2_surv"]) * param_base["cond_sp1_col"])-
            (param_base["sp2_surv"] * (1-param_base["cond_sp1_col"]))-
            (param_base["sp2_surv"] * param_base["cond_sp1_col"])
          trans[2,3] <- (1-param_base["sp2_surv"]) * param_base["cond_sp1_col"]
          trans[3,3] <- param_base["sp2_surv"] * (1-param_base["cond_sp1_col"])
          trans[4,3] <- param_base["sp2_surv"] * param_base["cond_sp1_col"]
          
          trans[1,4] <- 1 - (param_base["cond_sp1_surv"] * (1- param_base["cond_sp2_surv"])) - 
            (param_base["cond_sp2_surv"] * (1- param_base["cond_sp1_surv"])) - 
            (param_base["cond_sp1_surv"] * param_base["cond_sp2_surv"])
          trans[2,4] <- param_base["cond_sp1_surv"] * (1- param_base["cond_sp2_surv"])
          trans[3,4] <- param_base["cond_sp2_surv"] * (1- param_base["cond_sp1_surv"])
          trans[4,4] <- param_base["cond_sp1_surv"] * param_base["cond_sp2_surv"]
          
          ## calculate stable state distribution   
          
          stable_state <- trans %*% (eigen(trans)$vectors[, 1] / sum(eigen(trans)$vectors[, 1]))
          
          for(d in 1:length(repititions)){ 
            ## calculate occupancy
            for(L in 1:length(past_occ)) {  ## calculate occupancy
              if (past_occ[L] == "both") {
                temp_sp <- rmultinom(1,1, trans[,4])
                occ_sp[L,k] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
                
              } else if (past_occ[L] == "sp1") {
                temp_sp <- rmultinom(1,1, trans[,2])
                occ_sp[L,k] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
                
              } else if (past_occ[L] == "sp2") {
                temp_sp <- rmultinom(1,1, trans[,3])
                occ_sp[L,k] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
              } else if (past_occ[L] == "empty") {
                temp_sp <- rmultinom(1,1, trans[,1])
                occ_sp[L,k] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
              }
              
            } ## close L
            
            tt  <- sum_fun(occ_sp[,k])/ nsites
            names(tt) <- c("Empty", "Sp1_only", "Sp2_only", "Both")
            
            year = 1
            trans_temp <- matrix(data = 0, nrow = 4, ncol = years)
            trans_temp[,1] <- tt
            temp_occ <- occ_sp[,1]
            
            #  while(sum(tt >= (stable_state[,1]-stable_state[,1]*0.5) & tt <= (stable_state[,1]+stable_state[,1]*0.5)) != 4){
            while(abs(sum(tt - stable_state[,1])) > eps.all | (sum(abs((tt - stable_state[,1])) >= eps.one) > 1) | years < 250) {
              
              
              for(L in 1:length(temp_occ)) {  ## calculate occupancy
                if (temp_occ[L] == "Both") {
                  temp_sp <- rmultinom(1,1, trans[,4])
                  temp_occ[L] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
                  
                } else if (temp_occ[L] == "Sp1_only") {
                  temp_sp <- rmultinom(1,1, trans[,2])
                  temp_occ[L] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
                  
                } else if (temp_occ[L] == "Sp2_only") {
                  temp_sp <- rmultinom(1,1, trans[,3])
                  temp_occ[L] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
                } else if (temp_occ[L] == "Empty") {
                  temp_sp <- rmultinom(1,1, trans[,1])
                  temp_occ[L] <- dimnames(temp_sp)[[1]][which(temp_sp == 1)]
                }
                
              } ## close L
              
              tt  <- sum_fun(temp_occ)/ nsites
              names(tt) <- c("Empty", "Sp1_only", "Sp2_only", "Both")
              
              trans_temp[ ,year+1] <- tt
              
              year <- year + 1
            }
            
            ## store time to stable state
            
            state_trans[ , , d, j, k, i,q]     <- trans_temp
            #print(d)
            
            
          } ## close d
          
        } ## close k 
        
        
      } ## close j
      
      print(j)
      print(i)
      print(q)
      
    } ## close i
  } ## close q
)

####################################################
saveRDS(state_trans, file = "Time_to_equilibrium.RDS")