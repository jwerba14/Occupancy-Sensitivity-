library(tidyverse)
library(ggplot2)
#set up the theme for all the graphs
theme_set(theme_bw())
theme_update(axis.text.x = element_text(size = 12),
             axis.text.y = element_text(size = 12),
             axis.title.x = element_text(size = 14),
             axis.title.y = element_text(size = 14),
             legend.title = element_text(size = 12),
             legend.text = element_text(size = 10),
             legend.spacing = unit(0.25, "cm"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             panel.spacing = unit(0, "lines"),
             legend.key = element_rect(fill = "white"),
             panel.spacing.y = unit(-0.25, "lines"),
             panel.border = element_rect(colour = "black",
                                         fill = NA, size = 1),
             strip.text.x = element_text(size = 12, colour = "black",
                                         face = "bold"))

#trans <- readRDS("Time_to_equilibrium.RDS")
## dim1 = site type; dim 2 = time steps allowed; dim 3 = repititions; dim 4 = parameter;
##  dim 5 = parameter values; dim 6 = initial conditions; dim 7 = base values 

## have to break it up because too large to melt on my computer


#trans1 <- trans[,,,,,,1]
#trans2 <- trans[,,,,,,2]
#trans3 <- trans[,,,,,,3]


#saveRDS(trans1, "base1.RDS")
#saveRDS(trans2, "base2.RDS")
#saveRDS(trans3, "base3.RDS")



#rm(trans)
#rm(trans1)
#rm(trans2)
#rm(trans3)
#gc()

## distance formula

dd <- function(init_empty, init_sp1,init_sp2, init_both, equil_empty, equil_sp1, equil_sp2, equil_both){
  
  sqrt((init_empty-equil_empty)^2 + (init_sp1-equil_sp1)^2 + (init_sp2-equil_sp2)^2 + (init_both-equil_both)^2)
}



init_cond <- matrix(c(0.25,0.25,0.25,0.25, 
                      0.3,.3, 0.4, 0, 
                      0.5,0.1,0.3,0.1,
                      0.7, 0.2, 0.06, 0.04 ),
                    nrow = 4, ncol = 4, byrow = TRUE)

### work with each base value one at a time ##filesize problems ###

param_vals <- seq(0,0.9, length.out = 25) 

##### base 1 (0.1)
trans1 <- readRDS("base1.RDS")
## output we care about
## average years to equilibrium given 4 starting conditions and changes in param values

## to calculate sensitivity I split it into four pieces (one for each initial condition because I struggled with apply statements that would do this more efficiently
## so the rest is all very repetitive)

avg_year <- array(data = NA, dim = c(100, 25, 8, 4 ),
                  dimnames = list(NULL, NULL, dimnames(trans1)[[4]], NULL))

## to calculate sensitivity I split it into four pieces (one for each initial condition because I struggled with apply statements that would do this more efficiently
## so the rest is all very repetitive)


for(k in 1:dim(trans1)[[6]]){ ## init conditions
for(j in 1:dim(trans1)[[4]]){ ## parameters
  for(i in 1:dim(trans1)[[5]]){ ## parameter values
    
    avg_year[,i,j,k] <-apply(trans1[, , , j,i,k], 3, FUN = function(x) min(which(colSums(x) == 0)))
    
  }
}
}


zz <- apply(avg_year, c(2,4), FUN = function(x) colMeans(x))

## add error? ## 100 repititions
qq <- apply(avg_year, c(3,2,4), FUN = function(x) sd(x))

dimnames(zz)[[2]] <- param_vals
tt <- reshape2::melt(zz) %>% 
  rename(parameter = Var1, param_value= Var2, init=Var3, time = value)
tt$param_value <- round(tt$param_value, digits = 4)


dimnames(qq)[[2]] <- param_vals
qq <- reshape2::melt(qq) %>% 
  rename(parameter = Var1, param_value= Var2, init=Var3, sd = value)
qq$param_value <- round(tt$param_value, digits = 4)

tt <- left_join(tt, qq, by = c("parameter", "param_value", "init"))


labs <- c("Sp1 cond \n survival", "Sp2 cond \n survival",
          "Sp1 cond \n colonization", "Sp2 cond \n colonization",
          "Sp1 \n survival", "Sp2 \n survival", "Sp1 \n colonization", "Sp2 \n colonization")


tt <- tt %>% mutate(labels = plyr::mapvalues(parameter, from = c(
  "cond_sp1_surv", "cond_sp2_surv",
  "cond_sp1_col" , "cond_sp2_col",
  "sp1_surv", "sp2_surv", 
  "sp1_col", "sp2_col"
), to = labs))


## graph raw time
ggplot(tt, aes(param_value, time)) + 
  geom_point(aes(color = as.factor(init))) + geom_errorbar(aes(ymin = time-sd, ymax = time+sd)) +
  facet_grid(labels~as.factor(init)) + ggtitle("Time to equilibrium \n Base Value 0.1") +
  xlab("Parameter Value") + ylab("Timesteps to equilibrium") +
  scale_color_discrete("Initial Value") 


## calculate sensitivity -- does not include error

zz2 <- tt %>% 
  group_by(init, parameter) %>%
  arrange(param_value) %>% 
  mutate(diff = lag(time)-time, param_diff = lag(param_value)- param_value) %>%
  mutate(sens = diff/param_diff)


## graph sensitivity
ggplot(zz2, aes(param_value, abs(sens))) + 
  geom_point(aes(color = as.factor(init))) + 
  facet_grid(parameter~as.factor(init), scales = "free") + 
  ggtitle("Sensitivity \n Base Value 0.1") +
  xlab("Parameter Value") + ylab("Absolute Value Sensitivity") +
  scale_color_discrete("Initial Conditions")

## graph distance (from initial starting) to equilibrium vs. mean time to equilibrium

## calculate distance
distout <- array(data = NA, dim = c(1,length(unique(tt$parameter)), length(unique(tt$param_value)),
                                    length(unique(tt$init))), 
                 dimnames = list(NULL,  unique(tt$parameter), NULL, NULL))



for(l in 1:nrow(init_cond)){
  for(j in 1:length(unique(tt$parameter))){
    for(i in 1:length(unique(tt$param_value))){
      
      distout[,j,i,l] <- dd(init_empty = init_cond[l,1], 
                            init_sp1 = init_cond[l,2], 
                            init_sp2 = init_cond[l,3] , 
                            init_both = init_cond[l,4],
                            equil_empty = mean(trans1[1,1, ,j,i,l]) ,
                            equil_sp1 = mean(trans1[2,1, ,j,i,l]) ,
                            equil_sp2 = mean(trans1[3,1, ,j,i,l]) ,
                            equil_both = mean(trans1[4,1, ,j,i,l]) )
      
    }
  }
}

dimnames(distout)[[3]] <- param_vals

dis <- reshape2::melt(distout) %>%
  rename(parameter = Var2, param_value = Var3, init = Var4, dist = value)
dis$param_value <- round(dis$param_value, digits = 4)

dis <- left_join(dis, tt)

## graph distance
ggplot(dis, aes(dist, time)) + 
  geom_point(aes(color = as.factor(param_value))) + 
  facet_grid(parameter~as.factor(init), scales = "free") +
  ggtitle("Distance vs time") +
  xlab("Distance") + ylab("Time")

(g1 <- ggplot(dis %>% filter(init == 1), aes(param_value, time)) +
      geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 1), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 1), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 1 \n Base value 0.1")


## init 2
(g1 <- ggplot(dis %>% filter(init == 2), aes(param_value, time)) +
  geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
  xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 2), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 2), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 2 \n Base value 0.1")

## init 3
(g1 <- ggplot(dis %>% filter(init == 3), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 3), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 3), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 3 \n Base value 0.1")

##init 4
(g1 <- ggplot(dis %>% filter(init == 4), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value")  )

(g2 <- ggplot(dis %>% filter(init == 4), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium") )

(g3 <- ggplot(dis %>% filter(init == 4), aes(dist, time)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    ylab("Time to equilibrium") + xlab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value")  )

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 4 \n Base value 0.1")

rm(trans1)
gc()

#################### Base value 0.5
trans2 <- readRDS("base2.RDS")
## output we care about
## average years to equilibrium given 4 starting conditions and changes in param values

## to calculate sensitivity I split it into four pieces (one for each initial condition because I struggled with apply statements that would do this more efficiently
## so the rest is all very repetitive)

avg_year <- array(data = NA, dim = c(100, 25, 8, 4 ),
                  dimnames = list(NULL, NULL, dimnames(trans2)[[4]], NULL))

## to calculate sensitivity I split it into four pieces (one for each initial condition because I struggled with apply statements that would do this more efficiently
## so the rest is all very repetitive)


for(k in 1:dim(trans2)[[6]]){  ## init conditions
  for(j in 1:dim(trans2)[[4]]){ ## parameters
    for(i in 1:dim(trans2)[[5]]){  ## parameter values
      
      avg_year[,i,j,k] <-apply(trans2[, , , j,i,k], 3, FUN = function(x) min(which(colSums(x) == 0)))
      
    }
  }
}


zz <- apply(avg_year, c(2,4), FUN = function(x) colMeans(x))

## add error? ## 100 repititions
qq <- apply(avg_year, c(3,2,4), FUN = function(x) sd(x))

dimnames(zz)[[2]] <- param_vals
tt <- reshape2::melt(zz) %>% 
  rename(parameter = Var1, param_value= Var2, init=Var3, time = value)
tt$param_value <- round(tt$param_value, digits = 4)


dimnames(qq)[[2]] <- param_vals
qq <- reshape2::melt(qq) %>% 
  rename(parameter = Var1, param_value= Var2, init=Var3, sd = value)
qq$param_value <- round(tt$param_value, digits = 4)

tt <- left_join(tt, qq, by = c("parameter", "param_value", "init"))


labs <- c("Sp1 cond \n survival", "Sp2 cond \n survival",
          "Sp1 cond \n colonization", "Sp2 cond \n colonization",
          "Sp1 \n survival", "Sp2 \n survival", "Sp1 \n colonization", "Sp2 \n colonization")


tt <- tt %>% mutate(labels = plyr::mapvalues(parameter, from = c(
  "cond_sp1_surv", "cond_sp2_surv",
  "cond_sp1_col" , "cond_sp2_col",
  "sp1_surv", "sp2_surv", 
  "sp1_col", "sp2_col"
), to = labs))

## graph raw time
ggplot(tt, aes(param_value, time)) + 
  geom_point(aes(color = as.factor(init))) + 
  facet_grid(labels~as.factor(init)) + ggtitle("Time to equilibrium \n Base Value 0.1") +
  xlab("Parameter Value") + ylab("Timesteps to equilibrium")


## calculate sensitivity -- does not include error

zz2 <- tt %>% 
  group_by(init, parameter) %>%
  arrange(param_value) %>% 
  mutate(diff = lag(time)-time, param_diff = lag(param_value)- param_value) %>%
  mutate(sens = diff/param_diff)


## graph sensitivity
ggplot(zz2, aes(param_value, abs(sens))) + 
  geom_point(aes(color = as.factor(init))) + 
  facet_grid(parameter~as.factor(init), scales = "free") + 
  ggtitle("Sensitivity \n Base Value 0.5") +
  xlab("Parameter Value") + ylab("Absolute Value Sensitivity")

## graph distance (from initial starting) to equilibrium vs. mean time to equilibrium

## calculate distance
distout <- array(data = NA, dim = c(1,length(unique(tt$parameter)), length(unique(tt$param_value)),
                                    length(unique(tt$init))), 
                 dimnames = list(NULL,  unique(tt$parameter), NULL, NULL))



for(l in 1:nrow(init_cond)){
  for(j in 1:length(unique(tt$parameter))){
    for(i in 1:length(unique(tt$param_value))){
      
      distout[,j,i,l] <- dd(init_empty = init_cond[l,1], 
                            init_sp1 = init_cond[l,2], 
                            init_sp2 = init_cond[l,3] , 
                            init_both = init_cond[l,4],
                            equil_empty = mean(trans2[1,1, ,j,i,l]) ,
                            equil_sp1 = mean(trans2[2,1, ,j,i,l]) ,
                            equil_sp2 = mean(trans2[3,1, ,j,i,l]) ,
                            equil_both = mean(trans2[4,1, ,j,i,l]) )
      
    }
  }
}

dimnames(distout)[[3]] <- param_vals

dis <- reshape2::melt(distout) %>%
  rename(parameter = Var2, param_value = Var3, init = Var4, dist = value)
dis$param_value <- round(dis$param_value, digits = 4)

dis <- left_join(dis, tt)

## graph distance
ggplot(dis, aes(dist, time)) + 
  geom_point(aes(color = as.factor(param_value))) + 
  facet_grid(parameter~as.factor(init), scales = "free") +
  ggtitle("Distance vs time") +
  xlab("Distance") + ylab("Time")

(g1 <- ggplot(dis %>% filter(init == 1), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 1), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 1), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 1 \n Base value 0.5")


## init 2
(g1 <- ggplot(dis %>% filter(init == 2), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 2), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 2), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 2 \n Base value 0.5")

## init 3
(g1 <- ggplot(dis %>% filter(init == 3), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 3), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 3), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 3 \n Base value 0.5")

##init 4
(g1 <- ggplot(dis %>% filter(init == 4), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value")  + scale_y_log10())

(g2 <- ggplot(dis %>% filter(init == 4), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium") + scale_y_log10())

(g3 <- ggplot(dis %>% filter(init == 4), aes(dist, time)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    ylab("Time to equilibrium") + xlab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value") + scale_x_log10() + scale_y_log10() )

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 4 \n Base value 0.5")

rm(trans2)
gc()

################### Base value 0.7
trans3 <- readRDS("base3.RDS")
## output we care about
## average years to equilibrium given 4 starting conditions and changes in param values

## to calculate sensitivity I split it into four pieces (one for each initial condition because I struggled with apply statements that would do this more efficiently
## so the rest is all very repetitive)

avg_year <- array(data = NA, dim = c(100, 25, 8, 4 ),
                  dimnames = list(NULL, NULL, dimnames(trans3)[[4]], NULL))

## to calculate sensitivity I split it into four pieces (one for each initial condition because I struggled with apply statements that would do this more efficiently
## so the rest is all very repetitive)


for(k in 1:dim(trans3)[[6]]){  ## init conditions
  for(j in 1:dim(trans3)[[4]]){ ## parameters
    for(i in 1:dim(trans3)[[5]]){  ## parameter values
      
      avg_year[,i,j,k] <-apply(trans3[, , , j,i,k], 3, FUN = function(x) min(which(colSums(x) == 0)))
      
    }
  }
}


zz <- apply(avg_year, c(2,4), FUN = function(x) colMeans(x))

## add error? ## 100 repititions
qq <- apply(avg_year, c(2,3,4), FUN = function(x) sd(x))

qq <- apply(avg_year, c(3,2,4), FUN = function(x) sd(x))

dimnames(zz)[[2]] <- param_vals
tt <- reshape2::melt(zz) %>% 
  rename(parameter = Var1, param_value= Var2, init=Var3, time = value)
tt$param_value <- round(tt$param_value, digits = 4)


dimnames(qq)[[2]] <- param_vals
qq <- reshape2::melt(qq) %>% 
  rename(parameter = Var1, param_value= Var2, init=Var3, sd = value)
qq$param_value <- round(tt$param_value, digits = 4)

tt <- left_join(tt, qq, by = c("parameter", "param_value", "init"))


labs <- c("Sp1 cond \n survival", "Sp2 cond \n survival",
          "Sp1 cond \n colonization", "Sp2 cond \n colonization",
          "Sp1 \n survival", "Sp2 \n survival", "Sp1 \n colonization", "Sp2 \n colonization")


tt <- tt %>% mutate(labels = plyr::mapvalues(parameter, from = c(
  "cond_sp1_surv", "cond_sp2_surv",
  "cond_sp1_col" , "cond_sp2_col",
  "sp1_surv", "sp2_surv", 
  "sp1_col", "sp2_col"
), to = labs))


## graph raw time
ggplot(tt, aes(param_value, time)) + 
  geom_point(aes(color = as.factor(init))) + 
  facet_grid(labels~as.factor(init)) + ggtitle("Time to equilibrium \n Base Value 0.7") +
  xlab("Parameter Value") + ylab("Timesteps to equilibrium")


## calculate sensitivity -- does not include error

zz2 <- tt %>% 
  group_by(init, parameter) %>%
  arrange(param_value) %>% 
  mutate(diff = lag(time)-time, param_diff = lag(param_value)- param_value) %>%
  mutate(sens = diff/param_diff)


## graph sensitivity
ggplot(zz2, aes(param_value, sens)) + 
  geom_point(aes(color = as.factor(init))) + 
  facet_grid(parameter~as.factor(init), scales = "free") + 
  ggtitle("Sensitivity \n Base Value 0.7") +
  xlab("Parameter Value") + ylab("Elasticity")

## graph distance (from initial starting) to equilibrium vs. mean time to equilibrium

## calculate distance
distout <- array(data = NA, dim = c(1,length(unique(tt$parameter)), length(unique(tt$param_value)),
                                    length(unique(tt$init))), 
                 dimnames = list(NULL,  unique(tt$parameter), NULL, NULL))



for(l in 1:nrow(init_cond)){
  for(j in 1:length(unique(tt$parameter))){
    for(i in 1:length(unique(tt$param_value))){
      
      distout[,j,i,l] <- dd(init_empty = init_cond[l,1], 
                            init_sp1 = init_cond[l,2], 
                            init_sp2 = init_cond[l,3] , 
                            init_both = init_cond[l,4],
                            equil_empty = mean(trans3[1,1, ,j,i,l]) ,
                            equil_sp1 = mean(trans3[2,1, ,j,i,l]) ,
                            equil_sp2 = mean(trans3[3,1, ,j,i,l]) ,
                            equil_both = mean(trans3[4,1, ,j,i,l]) )
      
    }
  }
}

dimnames(distout)[[3]] <- param_vals

dis <- reshape2::melt(distout) %>%
  rename(parameter = Var2, param_value = Var3, init = Var4, dist = value)
dis$param_value <- round(dis$param_value, digits = 4)

dis <- left_join(dis, tt)

## graph distance
ggplot(dis, aes(dist, time)) + 
  geom_point(aes(color = as.factor(param_value))) + 
  facet_grid(parameter~as.factor(init), scales = "free") +
  ggtitle("Distance vs time") +
  xlab("Distance") + ylab("Time")

(g1 <- ggplot(dis %>% filter(init == 1), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 1), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 1), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 1 \n Base value 0.7")


## init 2
(g1 <- ggplot(dis %>% filter(init == 2), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 2), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 2), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 2 \n Base value 0.7")

## init 3
(g1 <- ggplot(dis %>% filter(init == 3), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value"))

(g2 <- ggplot(dis %>% filter(init == 3), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium"))

(g3 <- ggplot(dis %>% filter(init == 3), aes(time, dist)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("Time to equilibrium") + ylab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value"))

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 3 \n Base value 0.7")

##init 4
(g1 <- ggplot(dis %>% filter(init == 4), aes(param_value, time)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value")  + scale_y_log10())

(g2 <- ggplot(dis %>% filter(init == 4), aes(param_value, dist)) +
    geom_point() + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    xlab("parameter value") + ylab("distance to equilibrium") + scale_y_log10())

(g3 <- ggplot(dis %>% filter(init == 4), aes(dist, time)) +
    geom_point((aes(color = param_value))) + facet_wrap(~parameter, ncol = 1, scales = "free") + 
    ylab("Time to equilibrium") + xlab("distance to equilibrium") +
    scale_color_continuous("Parameter \n Value") + scale_x_log10() + scale_y_log10() )

gridExtra::grid.arrange(g1, g2, g3, ncol = 3, top = "Initial Occ 4 \n Base value 0.7")
