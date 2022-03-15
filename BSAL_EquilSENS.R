params <- data.frame(params = c (  "c_H", "c_S", "phi_Hb", 
                                 "phi_Hs", "phi_HL", "e_S", "e_L", "g_S", "d_L" ),
                     upr = c(0.95,0.95,0.98,0.2,0.15,0.1,0.1,0.95,0.1 ),
                     lwr = c(0.6,0.6,0.8,0.01,0.001,0.01,0.001,0.6,0.01),
                      mode = c( 0.8,0.8,0.95,0.08,0.05,0.05,0.01,0.8,0.06 ))

## try with larger difference in max/min
params <- data.frame(params = c (  "c_H", "c_S", "phi_Hb", 
                                   "phi_Hs", "phi_HL", "e_S", "e_L", "g_S", "d_L" ),
                     upr = c(0.99,0.99,0.98,0.5,0.5,0.5,0.5,0.99,0.5 ),
                     lwr = c(0.2,0.2,0.2,0.001,0.001,0.001,0.001,0.2,0.001),
                     mode = c( 0.8,0.8,0.95,0.08,0.05,0.05,0.01,0.8,0.06 ))


# There are 6 states:
# 1 = bh = No Bsal, no host
# 2 = bH = No Bsal, host is present
# 3 = sh = Small Bsal prevalence, no host
# 4 = sH = Small Bsal prevalence, host is present
# 5 = Lh = Large Bsal prevalence, no host
# 6 = LH = Large Bsal prevalence, host is present

sites <- c("bh", "bH", "sh", "sH", "Lh", "LH")


trans <- matrix(data = NA, nrow = length(sites), ncol = length(sites), dimnames = list(sites, sites) )


## storing transition states


equil_state <- array(data = 0, dim = c(length(sites), 15, nrow(params)),
                     dimnames = list(sites, NULL, params$params)) 

par_val <- array(data = 0, dim = c( 15, nrow(params)),
                     dimnames = list(NULL, params$params)) 

par_index = c("c_H", "c_S", "phi_Hb", 
            "phi_Hs", "phi_HL", "e_S", "e_L", "g_S", "d_L" )

# Transition matrix after time step 1
# Sites can be either bh or bH
# Pre-Bsal arrival
  
 for(i in 1:length(par_index)){
   
   #phi_bh <- params[which(params$params == "phi_bh"), "mode"]
   c_H <- params[which(params$params == "c_H"), "mode"]
   c_S <- params[which(params$params == "c_S"), "mode"]  
   phi_Hb <- params[which(params$params == "phi_Hb"), "mode"]
    
   phi_Hs <- params[which(params$params == "phi_Hs"), "mode"]
  phi_HL <- params[which(params$params == "phi_HL"), "mode"]
     e_S <- params[which(params$params == "e_S"), "mode"]
     e_L <- params[which(params$params == "e_L"), "mode"]
     g_S <- params[which(params$params == "g_S"), "mode"]
     d_L <- params[which(params$params == "d_L"), "mode"]
 
  lwr <- params$lwr[i]
  upr <- params$upr[i]
  
  val <- seq(lwr, upr, length.out = 15)
  
  for(j in 1:length(val)){
   
     assign(par_index[i], val[j] )
 ## fill in transition matrix
  
        trans[1,1] <- (1-c_H)*(1-c_S)           #bh --> bh (1-HostCol)*(1-BSalCol)
        trans[1,2] <- c_H*(1-c_S)               ## bh --> bH
        trans[1,3] <- (1-c_H)*c_S               ## bh --> sh
        trans[1,4] <- c_H * c_S                 ## bh --> sH
        trans[1,5] <- 0                         ## bh --> Lh
        trans[1,6] <- 0                         ## bh --> LH
        
        trans[2,1] <- (1-phi_Hb)*(1-c_S)   #bH -->bh (1- HostSurvNOBSal) *(1-BSalCol)
        trans[2,2] <- phi_Hb*(1-c_S)       #bH --> bH
        trans[2,3] <- (1-phi_Hb)*c_S       #bH -->sh
        trans[2,4] <- phi_Hb * c_S         #bH --> sH
        trans[2,5] <- 0                    #bH --> Lh                   
        trans[2,6] <- 0                    #bH --> LH
        
        
        trans[3,1] <- (1-c_H) * e_S                     #sh --> bh (1-HostCol)*BSALEXTSMALL
        trans[3,2] <- c_H * e_S                         #sh --> bH
        trans[3,3] <- (1-c_H) * (1-e_S) * (1-g_S)       #sh --> sh (1-HostCol)* (1-BSALEXT) * (1-BSALGROWTH)
        trans[3,4] <-  c_H * (1-e_S) * (1-g_S)          #sh --> sH
        trans[3,5] <-  (1-c_H) * g_S * (1-e_S)          #sh --> Lh 
        trans[3,6] <-  c_H * g_S * (1-e_S)              #sh --> LH
        
        trans[4,1] <- (1-phi_Hs) * e_S           #sH --> bh (1-HostSurvSmallBSAl) * BSALEXTSMALL
        trans[4,2] <- phi_Hs * e_S               #sH --> bH
        trans[4,3] <- (1-phi_Hs)*(1-e_S)*(1-g_S) #sH --> sh 
        trans[4,4] <- phi_Hs*(1-e_S)*(1-g_S)     #sH -->sH
        trans[4,5] <- (1-phi_Hs)*(1-e_S)*g_S     #sH --> Lh
        trans[4,6] <-  phi_Hs*(1-e_S)*g_S        #sH --> LH
        
        trans[5,1] <- (1-c_H)*e_L             #Lh --> bh(1-HostCol)* BSALEXTLARGE
        trans[5,2] <- c_H*e_L                 #Lh --> bH
        trans[5,3] <- (1-c_H)*d_L*(1-e_L)     #Lh --> sh(1-HostCol)*BSALDECLINE*(1-BSALEXTLARGE)
        trans[5,4] <- c_H*d_L*(1-e_L)         #Lh --> sH
        trans[5,5] <- (1-c_H)*(1-d_L)*(1-e_L) #Lh --> Lh
        trans[5,6] <- c_H*(1-d_L)*(1-e_L)     #Lh --> LH
        
        trans[6,1] <- (1-phi_HL)*e_L             #LH --> bh(1-HostSurvLargeBSal)*BSALEXTLARGE
        trans[6,2] <- phi_HL*e_L                 #LH --> bH
        trans[6,3] <- (1-phi_HL)*(1-e_L)*d_L     #LH --> sh 
        trans[6,4] <- phi_HL*(1-e_L)*d_L         #LH --> sH
        trans[6,5] <- (1-phi_HL)*(1-e_L)*(1-d_L) #LH --> Lh
        trans[6,6] <- phi_HL*(1-e_L)*(1-d_L)     #LH --> LH
    
        ## calculate stable state distribution   
        tt <- t(trans)
        
        stable_state <- tt %*% (eigen(tt)$vectors[, 1] / sum(eigen(tt)$vectors[, 1]))
        par_val[j,i] <- val[j] ## store parameter value
        equil_state[ ,j,i] <- as.numeric(stable_state[,1]) 
      
    
  } #close j (values)
 } #close i (parameters)


## calculate sensitivity



## change in equilibrium
zz <- reshape2::melt(equil_state) %>% rename(state = Var1, param = Var3, value= value, iter = Var2)

tt <- zz %>% 
  group_by(state, param) %>% 
  arrange(iter) %>%
  mutate(change =c(NA, diff(value, lag = 1)))


pp <- reshape2::melt(par_val) %>% rename(iter=Var1, param= Var2, par_val =value)

rr <- pp %>% 
  group_by(param) %>% 
  arrange(iter) %>%
  mutate(par_diff =c(NA, diff(par_val, lag = 1)))

qq <- left_join(tt, rr)
qq <- qq %>% mutate(sens = change/par_diff)


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


par.lab <- c("Host Colonization", "BSAL Colonization", "Host Survival",
             "Host Survival \n small BSAL", "Host Survival \n large BSAL load", 
             "BSAL Extinction \n when small", "BSAL Extinction \n when large",
             "Growth of BSAL infestation \n from small to large", 
             "Decrease of BSAL infestation \n from large to small ")


names(par.lab) <- c("c_H", "c_S", "phi_Hb", "phi_Hs", "phi_HL", "e_S", "e_L", "g_S", "d_L")


ggplot(qq, aes(par_val, sens)) + 
  geom_point(aes(color=state)) + 
  facet_wrap(~param, scales="free_x", labeller = labeller(param = par.lab)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(panel.spacing.x = unit(0.5, "lines")) + 
  ggtitle("BSAL sensitivity") +
  xlab("Parameter Values") + ylab("Elasticity") + 
  scale_color_discrete(labels = c("Empty", "Host Only", "BSAL small", "BSAL small +Host", 
                                  "BSAL large", "BSAL large +Host"))


## want a column that says proportion of sites with host regardless of disease state

nd <- qq 

dum <- grep("H", nd$state)
nd$dummy <- 0
nd$dummy[dum] <- 1

nd <- nd %>% filter(dummy == 1) %>% group_by(iter,param) %>% mutate(tot_host = sum(value))

ggplot(qq, aes(par_val, value)) + 
  geom_point(aes(color=state)) + 
  geom_point(data = nd, aes(par_val, tot_host), shape =17) +
  facet_wrap(~param, scales="free_x",  labeller = labeller(param = par.lab)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(panel.spacing.x = unit(0.5, "lines")) + 
  ggtitle("BSAL occupancy") +
  xlab("Parameter Values") + ylab("Occupancy") + 
  scale_color_discrete(labels = c("Empty", "Host Only", "BSAL small", "BSAL small +Host", 
                                  "BSAL large", "BSAL large +Host"))

## sensitivity of total host occupancy
dd <- nd %>%
  group_by(param,iter) %>%
  slice(1) %>%
  ungroup %>%
  group_by(param) %>%
  arrange(iter) %>%
  mutate(tot_change = c(NA, diff(tot_host, lag = 1))) %>%
  mutate(sens = tot_change/par_diff)




ggplot(dd, aes(par_val, sens)) + 
  geom_point(shape = 17) + 
  facet_wrap(~param, scales="free_x" , labeller = labeller(param = par.lab)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(panel.spacing.x = unit(0.5, "lines")) + 
  ggtitle("BSAL Host total occupancy") +
  xlab("Parameter Values") + ylab("Elasticity") 
