source("equilibrium_sensitivity.R")

## distance formula

dd <- function(init_empty, init_sp1,init_sp2, init_both, equil_empty, equil_sp1, equil_sp2, equil_both){
  
  sqrt((init_empty-equil_empty)^2 + (init_sp1-equil_sp1)^2 + (init_sp2-equil_sp2)^2 + (init_both-equil_both)^2)
}


init_cond <- matrix(c(0.25,0.25,0.25,0.25, 
                      0.3,.3, 0.4, 0, 
                      0.5,0.1,0.3,0.1,
                      0.7, 0.2, 0.06, 0.04 ),
                    nrow = 4, ncol = 4, byrow = TRUE)


output <- array(data = NA, dim = c(1,length(params), length(param_vals), length(base_val), nrow(init_cond)), 
                dimnames = list(NULL,  params, NULL, NULL,NULL))

for(l in 1:nrow(init_cond)){
  for(k in 1:length(base_val)){                                  
    for(j in 1:length(params)){
      for(i in 1:length(param_vals)){
        
        output[,j,i,k,l] <- dd(init_empty = init_cond[l,1], 
                               init_sp1 = init_cond[l,2], 
                               init_sp2 = init_cond[l,3] , 
                               init_both = init_cond[l,4],
                               equil_empty =equil_state[ 1,j,i,k] ,
                               equil_sp1 =equil_state[ 2,j,i,k] ,
                               equil_sp2 =equil_state[ 3,j,i,k] ,
                               equil_both =equil_state[ 4,j,i,k] )
        
      }
    }
  }
}
dimnames(output)[[3]] <- param_vals
dimnames(output)[[4]] <- base_val





##### graph distance/ time to equilibrium
sens <- read.csv("sensitivity_base0.1.csv") %>% dplyr::select(-X) 
sens_0.3 <- read.csv("sensitivity_base0.3.csv") %>% dplyr::select(-X) 
sens_0.5 <- read.csv("sensitivity_base0.5.csv") %>% dplyr::select(-X) 
sens_0.7 <- read.csv("sensitivity_base0.7.csv") %>% dplyr::select(-X) 
sens_0.9 <- read.csv("sensitivity_base0.9.csv") %>% dplyr::select(-X) 
combined_sens <- rbind(sens, sens_0.3)
comb2 <- rbind(combined_sens, sens_0.5)
combined_sens <- rbind(comb2, sens_0.7)
comb2 <- rbind(combined_sens, sens_0.9)
comb2$param_value <- round(comb2$param_value, digits = 4)

## initial cond 1
sens_init1 <- comb2 %>% filter(init == "init1") %>% dplyr::select(-init)

base1_init1 <- read.csv("time_baseparam0.1.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.1)
base3_init1 <- read.csv("time_baseparam0.3.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.3)
base5_init1 <- read.csv("time_baseparam0.5.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.5)
base7_init1 <- read.csv("time_baseparam0.7.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.7)
base9_init1 <- read.csv("time_baseparam0.9.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.9)

time_all <- rbind(base1_init1, base3_init1)
time2 <- rbind(time_all, base5_init1) 
time_all <- rbind(time2, base7_init1)
time2 <- rbind(time_all, base9_init1) 
time2$param_value <- round(time2$param_value, digits = 4)

equil_out <- reshape2::melt(output[1,,-1, ,1]) %>% 
  rename(parameter = Var1, param_value= Var2, base_param = Var3, distance = value)

equil_out$param_value <- round(equil_out$param_value, digits = 4)

init1_all <- left_join(time2,sens_init1 )
init1_fin <- left_join(init1_all, equil_out) %>% filter(param_value != 0)

##
gg_1 <- ggplot(init1_fin, aes(param_value,distance)) + 
  geom_line(aes(color = as.factor(base_param))) + facet_wrap(~parameter) + ggtitle("Initial Condition 1")
print(gg_1)


gg_time <- ggplot(init1_fin, aes(time,distance)) + 
  geom_point(aes(color = param_value)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 1")
print(gg_time)

gg_time2 <- ggplot(init1_fin, aes(param_value, distance)) + 
  geom_point(aes(color = time)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 1")
print(gg_time2)




#initial cond 2
sens_init2 <- comb2 %>% filter(init == "init2") %>% dplyr::select(-init)


base1_init2 <- read.csv("time_baseparam0.12.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.1)
base3_init2 <- read.csv("time_baseparam0.32.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.3)
base5_init2 <- read.csv("time_baseparam0.52.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.5)
base7_init2 <- read.csv("time_baseparam0.72.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.7)
base9_init2 <- read.csv("time_baseparam0.92.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.9)

time_all2 <- rbind(base1_init2, base3_init2)
time22 <- rbind(time_all2, base5_init2) 
time_all2 <- rbind(time22, base7_init2)
time22 <- rbind(time_all2, base9_init2) 
time22$param_value <- round(time22$param_value, digits = 4)

equil_out2 <- reshape2::melt(output[1,,-1, ,2]) %>% 
  rename(parameter = Var1, param_value= Var2, base_param = Var3, distance = value)

equil_out2$param_value <- round(equil_out2$param_value, digits = 4)

init2_all <- left_join(time22,sens_init2 )
init2_fin <- left_join(init2_all, equil_out2) %>% filter(param_value != 0)

##
gg_2 <- ggplot(init2_fin, aes(param_value,distance)) + 
  geom_line(aes(color = as.factor(base_param))) + facet_wrap(~parameter) + ggtitle("Initial Condition 2")
print(gg_2)


gg_time2 <- ggplot(init2_fin, aes(time,distance)) + 
  geom_point(aes(color = param_value)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 2")
print(gg_time2)

gg_time22 <- ggplot(init2_fin, aes(param_value, distance)) + 
  geom_point(aes(color = time)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 2")
print(gg_time22)


# initial cond 3
sens_init3 <- comb2 %>% filter(init == "init3") %>% dplyr::select(-init)


base1_init3 <- read.csv("time_baseparam0.13.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.1)
base3_init3 <- read.csv("time_baseparam0.33.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.3)
base5_init3 <- read.csv("time_baseparam0.53.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.5)
base7_init3 <- read.csv("time_baseparam0.73.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.7)
base9_init3 <- read.csv("time_baseparam0.93.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.9)

time_all3 <- rbind(base1_init3, base3_init3)
time23 <- rbind(time_all3, base5_init3) 
time_all3 <- rbind(time23, base7_init3)
time23 <- rbind(time_all3, base9_init3) 
time23$param_value <- round(time23$param_value, digits = 4)

equil_out3 <- reshape2::melt(output[1,,-1, ,3]) %>% 
  rename(parameter = Var1, param_value= Var2, base_param = Var3, distance = value)

equil_out3$param_value <- round(equil_out3$param_value, digits = 4)

init3_all <- left_join(time23,sens_init3 )
init3_fin <- left_join(init3_all, equil_out3) %>% filter(param_value != 0)

##
gg_3 <- ggplot(init3_fin, aes(param_value,distance)) + 
  geom_line(aes(color = as.factor(base_param))) + facet_wrap(~parameter) + ggtitle("Initial Condition 3")
print(gg_3)


gg_time3 <- ggplot(init3_fin, aes(time,distance)) + 
  geom_point(aes(color = param_value)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 3")
print(gg_time3)

gg_time23 <- ggplot(init2_fin, aes(param_value, distance)) + 
  geom_point(aes(color = time)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 3")
print(gg_time23)



# initial cond 4
sens_init4 <- comb2 %>% filter(init == "init4") %>% dplyr::select(-init)


base1_init4 <- read.csv("time_baseparam0.14.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.1)
base3_init4 <- read.csv("time_baseparam0.34.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.3)
base5_init4 <- read.csv("time_baseparam0.54.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.5)
base7_init4 <- read.csv("time_baseparam0.74.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.7)
base9_init4 <- read.csv("time_baseparam0.94.csv") %>% dplyr::select(-X) %>% mutate(base_param = 0.9)

time_all4 <- rbind(base1_init4, base3_init4)
time24 <- rbind(time_all4, base5_init4) 
time_all4 <- rbind(time24, base7_init4)
time24 <- rbind(time_all4, base9_init4) 
time24$param_value <- round(time24$param_value, digits = 4)

equil_out4 <- reshape2::melt(output[1,,-1, ,4]) %>% 
  rename(parameter = Var1, param_value= Var2, base_param = Var3, distance = value)

equil_out4$param_value <- round(equil_out4$param_value, digits = 4)

init4_all <- left_join(time24,sens_init4 )
init4_fin <- left_join(init4_all, equil_out4) %>% filter(param_value != 0)

##
gg_4 <- ggplot(init4_fin, aes(param_value,distance)) + 
  geom_line(aes(color = as.factor(base_param))) + facet_wrap(~parameter) + ggtitle("Initial Condition 4")
print(gg_4)


gg_time4 <- ggplot(init4_fin, aes(time,distance)) + 
  geom_point(aes(color = param_value)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 4")
print(gg_time4)

gg_time24 <- ggplot(init4_fin, aes(param_value, distance)) + 
  geom_point(aes(color = time)) + facet_wrap(base_param~parameter) + ggtitle("Initial Condition 4")
print(gg_time24)


gg_1 +theme(legend.position = c(0.8,0.2)) + ylab("distance to equilibrium")
gg_2 +theme(legend.position = c(0.8,0.2)) + ylab("distance to equilibrium")
gg_3 +theme(legend.position = c(0.8,0.2)) + ylab("distance to equilibrium")
gg_4 +theme(legend.position = c(0.8,0.2)) + ylab("distance to equilibrium")

print(gg_time4)
