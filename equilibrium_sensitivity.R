library(tidyverse)
library(stringr)

params <- c("cond_sp1_surv","cond_sp2_surv", "cond_sp2_col", "cond_sp1_col", "sp1_surv", "sp2_surv","sp1_col", "sp2_col")


## range of parameter values
param_vals <- seq(0,0.9, length.out = 50) 


## parameter values for parameters that are not changing
## could loop over a set of these values as well

base_val = c(0.1, 0.3,0.5,0.7,0.9)



## create transition matrix
sites <- c("Empty", "Sp1_only", "Sp2_only", "Both")

trans <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(sites, sites) )




## storing transition states


equil_state <- array(data = 0, dim = c(4, length(params), length(param_vals), length(base_val)),
                     dimnames = list(sites,  params, NULL, NULL)) 

   
for(i in 1:length(base_val)) {   
    ## the next loop is over the 8 different parameters  
for(j in 1:length(params)){
      
      ## all the parameters are set to the base value that was set outside of the loops
      param_base <- c(cond_sp1_surv = base_val[i],
                      cond_sp2_surv = base_val[i],
                      cond_sp2_col = base_val[i],
                      cond_sp1_col = base_val[i],
                      sp1_surv = base_val[i],
                      sp2_surv = base_val[i],
                      sp1_col = base_val[i],
                      sp2_col = base_val[i])
      

      
      
      ## the next loop calculates stable state
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
        
        equil_state[ ,j,k,i] <- as.numeric(stable_state[,1])
      } #close k
    } #close j
} # close i

## calculate sensitivity 
## calculate change in parameter value
delta_param <- numeric(length = length(param_vals)-1)

for(v in 1:length(param_vals)-1){
  delta_param[v] <- param_vals[v+1] - param_vals[v] ## lol these are all the same bc i used seq to make param vals
}

## well i suck at apply so here we co again

theme_set(theme_bw()) 
theme_update(axis.text.x = element_text(size = 8),
             axis.text.y = element_text(size = 8),
             axis.title.x = element_text(size = 10),
             axis.title.y = element_text(size = 10),
             legend.title = element_text(size = 10),
             legend.text = element_text(size = 8),
             legend.spacing = unit(0.25, "cm"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             panel.spacing = unit(0, "lines"),
             legend.key = element_rect(fill = "white"),
             panel.spacing.y = unit(-0.25, "lines"),
             panel.border = element_rect(colour = "black", 
                                         fill = NA, size = 1),
             strip.text.x = element_text(size = 10, colour = "black", 
                                         face = "bold"))

ss <- seq(1,50)
tt <- seq(1,5)
zz <- reshape2::melt(equil_state) %>% rename(Site = Var1, param = Var2) %>% 
  mutate(param_value = plyr::mapvalues(Var3, from = ss, to = param_vals)) %>%
  mutate(base_val =plyr::mapvalues(Var4, from = tt, to = base_val)) 

qq <- zz %>% 
  group_by(base_val, param, Site) %>%
  mutate(elas = (value - lag(value))/delta_param[1])

## separate conditional param from non-conditional so that it is easier to look at
site.lab <- c("Empty", "Species 1 Only", "Species 2 Only", "Co-occurence")
names(site.lab) <- c("Empty", "Sp1_only", "Sp2_only", "Both")
con.lab <- c("Conditional \n Survival (sp1)", "Conditional \n Survival (sp2)", 
             "Conditional \n Colonization (sp1)", "Conditional \n Colonization (sp2)" )
names(con.lab) <- c("cond_sp1_surv", "cond_sp2_surv", 
                    "cond_sp1_col", "cond_sp2_col" )

qq %>% filter(grepl("cond", param)) %>% {
ggplot(., aes(param_value, elas)) + 
  geom_line(aes(color = as.factor(base_val),
                linetype = as.factor(base_val)), lwd = 1.5)  +
  facet_grid(~param~Site, labeller = labeller(Site = site.lab, param = con.lab)) +
  xlab("Parameter Value") + theme(panel.spacing = unit(0.5, "lines"))+
  ylab("Elasticity") + ggtitle("Conditional Parameters") +
    scale_color_discrete(name = "Base Value") + 
    scale_linetype_discrete(name = "Base Value")
    

}

## ok but actually can just filter to one sp only because everything is symmetric

qq %>% filter(grepl("cond", param), Site != "Sp2_only") %>% {
  ggplot(., aes(param_value, elas)) + 
    geom_line(aes(color = as.factor(base_val),
                  linetype = as.factor(base_val)), lwd = 1.5)  +
    facet_grid(~Site~param, labeller = labeller(Site = site.lab, param = con.lab)) +
    xlab("Parameter Value") + theme(panel.spacing = unit(0.5, "lines"))+
    ylab("Elasticity") + ggtitle("Conditional Parameters") +
    scale_color_discrete(name = "Base Value") + 
    scale_linetype_discrete(name = "Base Value")
  
  
}

## non conditional parameters
par.lab <- c("Survival (sp1)","Colonization (sp1)" )
names(par.lab) <- c("sp1_surv", "sp1_col")

qq %>% filter(!grepl("cond", param), grepl("sp1", param)) %>% {
  ggplot(., aes(param_value, elas)) + 
    geom_line(aes(color = as.factor(base_val),
                  linetype = as.factor(base_val)), lwd = 1.5)  +
    facet_grid(~Site~param, labeller = labeller(Site = site.lab, 
                                                 param = par.lab)) +
    xlab("Parameter Value") + theme(panel.spacing = unit(0.5, "lines"))+
    ylab("Elasticity") + ggtitle("Non-conditional Parameters") +
    scale_color_discrete(name = "Base Value") + 
    scale_linetype_discrete(name = "Base Value")
  
  
}


##### graph occupancy instead of sensitivity

qq %>% filter(grepl("cond", param)) %>% {
  ggplot(., aes(param_value, value)) + 
    geom_line(aes(color = as.factor(base_val),
                  linetype = as.factor(base_val)), lwd = 1.5)  +
    facet_grid(~param~Site, labeller = labeller(Site = site.lab, param = con.lab)) +
    xlab("Parameter Value") + theme(panel.spacing = unit(0.5, "lines"))+
    ylab("Proportion Occupied") + ggtitle("Conditional Parameters") +
    scale_color_discrete(name = "Base Value") + 
    scale_linetype_discrete(name = "Base Value")
  
  
}


