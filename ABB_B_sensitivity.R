library(tidyverse)
## abiotic-biotic (moss) sensitivity 
#setwd("C:/Users/jwerba/OneDrive - DOI/Desktop/SHENSAL_OUTPUT")
out <- readRDS("params.RDS")
out <- out %>% filter(model == "moss.Abiotic_Biotic", year == 15)



#par.med <- out %>% group_by(param) %>% 
# summarize(enframe(quantile(value, probs = c(0.5, 0.875, 0.125)), "quantile"))

par.med <- out %>% group_by(param) %>% 
  summarize(enframe(quantile(value, probs = c(0.5, 0.975, 0.025)), "quantile"))

par <- par.med %>% pivot_wider(names_from = quantile, values_from = value)
colnames(par) <- c("param", "med", "upr", "lwr")
#par.q75 <- reshape2::melt(out$q97.5[1:32]) %>% rename(param = L1, quant = value)
rm(out)
gc()
par_index <- c("alpha.colEmptyC",  "alpha.colEmptyS",  "alpha.colOccC" ,   "alpha.colOccS",                
               "alpha.survAloneC", "alpha.survAloneS", "alpha.survBothC",  "alpha.survBothS",
               "beta.colEmptyC1",  "beta.colEmptyC3" , "beta.colEmptyS1",  "beta.colEmptyS3",
               "beta.colOccC1",    "beta.colOccC3" ,   "beta.colOccS1" ,   "beta.colOccS3",
               "beta.survAloneC1", "beta.survAloneC3", "beta.survAloneS1", "beta.survAloneS3",
               "beta.survBothC1",  "beta.survBothC3" , "beta.survBothS1",  "beta.survBothS3" )


par <- par %>%
  ungroup() %>%
  filter(param %in% par_index) %>% droplevels()

par_index == unique(par$param)
#par <- left_join(par.med,par.q75) %>% mutate(upr= quantile(value, ), lwr = med-quant)


## create transition matrix
sites <- c("Empty", "Sp1_only", "Sp2_only", "Both")

trans <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(sites, sites) )


## storing transition states

## env parameters are centered, 
#so 0 is at mean tmp/moss, -1 is 1 std dv below and 1 is 1 std above
## combos assume other env variable is at mean
temp <- c(-1,1,0,0,0) 
moss <- c(0,0,0,-1,1)



equil_state <- array(data = 0, dim = c(length(sites), 15, nrow(par), length(temp)),
                     dimnames = list(sites, NULL, par$param, NULL)) 

par_val <- array(data = 0, dim = c( 15, nrow(par), length(temp)),
                 dimnames = list(NULL, par$param, NULL)) 


for(k in 1:length(temp)){ 
  tmp <- temp[k]
  mss <- moss[k]
  
  for(i in 1:length(par_index)){
    print(i)
    alpha.colEmptyS <- (par[which(par$param == "alpha.colEmptyS"), "med"]) %>% unlist()
    beta.colEmptyS1 <- (par[which(par$param == "beta.colEmptyS1"), "med"]) %>% unlist()
    beta.colEmptyS3 <- (par[which(par$param == "beta.colEmptyS3"), "med"]) %>% unlist()
    alpha.colEmptyC <- (par[which(par$param == "alpha.colEmptyC"), "med"])%>% unlist()
    beta.colEmptyC1 <- (par[which(par$param == "beta.colEmptyC1"), "med"])%>% unlist()
    beta.colEmptyC3 <- (par[which(par$param == "beta.colEmptyC3"), "med"])%>% unlist()
    alpha.survAloneS <- (par[which(par$param == "alpha.survAloneS"), "med"]) %>% unlist()
    beta.survAloneS1 <- (par[which(par$param == "beta.survAloneS1"), "med"])%>% unlist()
    beta.survAloneS3 <- (par[which(par$param == "beta.survAloneS3"), "med"])%>% unlist()
    alpha.survAloneC <- (par[which(par$param == "alpha.survAloneC"), "med"])%>% unlist()
    beta.survAloneC1 <- (par[which(par$param == "beta.survAloneC1"), "med"])%>% unlist()
    beta.survAloneC3 <- (par[which(par$param == "beta.survAloneC3"), "med"])%>% unlist()
    alpha.colOccS <- (par[which(par$param == "alpha.colOccS"), "med"])%>% unlist()
    beta.colOccS1 <- (par[which(par$param == "beta.colOccS1"), "med"])%>% unlist()
    beta.colOccS3 <- (par[which(par$param == "beta.colOccS3"), "med"])%>% unlist()
    alpha.colOccC <- (par[which(par$param == "alpha.colOccC"), "med"])%>% unlist()
    beta.colOccC1 <- (par[which(par$param == "beta.colOccC1"), "med"])%>% unlist()
    beta.colOccC3 <- (par[which(par$param == "beta.colOccC3"), "med"])%>% unlist()
    alpha.survBothS <- (par[which(par$param == "alpha.survBothS"), "med"])%>% unlist()
    beta.survBothS1 <- (par[which(par$param == "beta.survBothS1"), "med"])%>% unlist()
    beta.survBothS3 <- (par[which(par$param == "beta.survBothS3"), "med"])%>% unlist()
    alpha.survBothC <- (par[which(par$param == "alpha.survBothC"), "med"])%>% unlist()
    beta.survBothC1 <- (par[which(par$param == "beta.survBothC1"), "med"])%>% unlist()
    beta.survBothC3 <- (par[which(par$param == "beta.survBothC3"), "med"])%>% unlist()
    
    
    lwr <- par$lwr[i]
    upr <- par$upr[i]
    
    val <- seq(lwr, upr, length.out = 15)
    
    
    
    for(j in 1:length(val)){
      
      assign(par_index[i], val[j] )
      ## fill in transition matrix
      
      sp1_col <- plogis(alpha.colEmptyS + beta.colEmptyS1*tmp + beta.colEmptyS3*mss)
      print(sp1_col)
      sp2_col <- plogis(alpha.colEmptyC + beta.colEmptyC1*tmp + beta.colEmptyC3*mss)
      print(sp2_col)
      sp1_surv <- plogis(alpha.survAloneS + beta.survAloneS1*tmp + beta.survAloneS3*mss)
      print(sp1_surv)
      sp2_surv <- plogis(alpha.survAloneC + beta.survAloneC1*tmp + beta.survAloneC3*mss)
      print(sp2_surv)
      
      cond_sp1_col <- plogis(alpha.colOccS + beta.colOccS1*tmp + beta.colOccS3*mss)
      print(cond_sp1_col)
      cond_sp2_col <- plogis(alpha.colOccC + beta.colOccC1*tmp + beta.colOccC3*mss)
      print(cond_sp2_col)
      cond_sp1_surv <- plogis(alpha.survBothS + beta.survBothS1*tmp + beta.survBothS3*mss)
      print(cond_sp1_surv)
      cond_sp2_surv <-plogis(alpha.survBothC + beta.survBothC1*tmp + beta.survBothC3*mss) 
      print(cond_sp2_surv)
      
      trans[1,1] <- 1- (sp1_col * (1- sp2_col)) - (sp2_col * (1- sp1_col)) - (sp2_col * sp1_col)
      trans[2,1] <- sp1_col * (1- sp2_col)
      trans[3,1] <- sp2_col * (1- sp1_col)
      trans[4,1] <- sp2_col *  sp1_col
      
      
      trans[1,2] <- 1 - (sp1_surv * (1-cond_sp2_col)) - ((1-sp1_surv) * cond_sp2_col) - (sp1_surv * cond_sp2_col)
      trans[2,2] <- sp1_surv * (1-cond_sp2_col)
      trans[3,2] <- (1-sp1_surv) * cond_sp2_col
      trans[4,2] <- sp1_surv * cond_sp2_col
      
      trans[1,3] <- 1 - (sp2_surv * (1-cond_sp1_col)) - ((1-sp2_surv) * cond_sp1_col) - (sp2_surv * cond_sp1_col)
      trans[2,3] <- sp2_surv * (1-cond_sp1_col)
      trans[3,3] <- (1-sp2_surv) * cond_sp1_col
      trans[4,3] <- sp2_surv * cond_sp1_col
      
      trans[1,4] <- 1 - (cond_sp1_surv * (1- cond_sp2_surv)) - (cond_sp2_surv * (1- cond_sp1_surv)) - (cond_sp1_surv * cond_sp2_surv)
      trans[2,4] <- cond_sp1_surv * (1- cond_sp2_surv)
      trans[3,4] <- cond_sp2_surv * (1- cond_sp1_surv)
      trans[4,4] <- cond_sp1_surv * cond_sp2_surv
      
      ## calculate stable state distribution   
      
      stable_state <- trans %*% (eigen(trans)$vectors[, 1] / sum(eigen(trans)$vectors[, 1]))
      par_val[j,i,k] <- val[j] ## store parameter value
      equil_state[ ,j,i,k] <- as.numeric(stable_state[,1]) 
      print(j)
      
    } #close j (values)
  } #close i (parameters)
} #close k (temp)

## calculate sensitivity

## change in equilibrium
zz <- reshape2::melt(equil_state) %>% 
  rename(state = Var1, param = Var3, value= value, iter = Var2, env = Var4)

tt <- zz %>% 
  group_by(state, param, env) %>% 
  arrange(iter) %>%
  mutate(change =c(NA, diff(value, lag = 1)))


pp <- reshape2::melt(par_val) %>% rename(iter=Var1, param= Var2, par_val =value, env = Var3)

rr <- pp %>% 
  group_by(param, env) %>% 
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


qq2 <- qq %>% group_by(param, env) %>%
  mutate(elasc = abs(sens))

qq1 <- qq2 %>% group_by(param) %>% filter(env == 1, elasc > 0.01) %>% droplevels()
(me1 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic (moss)  \n mean moss cover -1 std temperature 
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 2, elasc > 0.01) %>% droplevels()
(me2 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic (moss)  \n mean moss cover +1 std temperature
  \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 3, elasc > 0.01) %>% droplevels()
(me3 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic (moss)  \n mean moss cover mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 4, elasc > 0.01) %>% droplevels()
(me4 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic (moss)  \n -1 std moss cover mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none") )



qq1 <- qq2 %>% group_by(param) %>% filter(env == 5, elasc > 0.01) %>% droplevels()
(me5 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic (moss)  \n +1 std moss cover mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none") )




#nd <- qq %>% filter(state == "Sp2_only" | state == "Both")
#nd <- nd%>% group_by(iter,param, env) %>% mutate(tot_host = sum(value))

nd1 <- qq %>% filter(state == "Sp1_only" | state == "Both")%>% 
  group_by(iter,param, env) %>%
  mutate(tot_host = sum(value))

nd1 <- nd1 %>% group_by(param, env) %>%
  mutate(occ_change = max(tot_host) - min(tot_host))

nd2 <- nd1 %>% group_by(param) %>% filter(env == 1, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 1, param %in% nd2$param)
(mo1 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy") 
)


q <- gridExtra::arrangeGrob(me1, mo1,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 2, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 2, param %in% nd2$param)
(mo2 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy")
)

q <- gridExtra::arrangeGrob(me2, mo2,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)

nd2 <- nd1 %>% group_by(param) %>% filter(env == 3, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 3, param %in% nd2$param)
(mo3 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))


q <- gridExtra::arrangeGrob(me3, mo3,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 4, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 4, param %in% nd2$param)
(mo4 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))


q <- gridExtra::arrangeGrob(me4, mo4,
                            layout_matrix = matrix(c(1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 5, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 5, param %in% nd2$param)
(mo5 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))

q <- gridExtra::arrangeGrob(me5, mo5,
                            layout_matrix = matrix(c(1, 1,1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)


### abiotic-biotic RMI
out <- readRDS("params.RDS")
out <- out %>% filter(model == "Abiotic_Biotic", year == 15)



#par.med <- out %>% group_by(param) %>% 
# summarize(enframe(quantile(value, probs = c(0.5, 0.875, 0.125)), "quantile"))

par.med <- out %>% group_by(param) %>% 
  summarize(enframe(quantile(value, probs = c(0.5, 0.975, 0.025)), "quantile"))

par <- par.med %>% pivot_wider(names_from = quantile, values_from = value)
colnames(par) <- c("param", "med", "upr", "lwr")
#par.q75 <- reshape2::melt(out$q97.5[1:32]) %>% rename(param = L1, quant = value)
rm(out)
gc()
par_index <- c("alpha.colEmptyC",  "alpha.colEmptyS",  "alpha.colOccC" ,   "alpha.colOccS",                
               "alpha.survAloneC", "alpha.survAloneS", "alpha.survBothC",  "alpha.survBothS",
               "beta.colEmptyC1",  "beta.colEmptyC2" , "beta.colEmptyS1",  "beta.colEmptyS2",
               "beta.colOccC1",    "beta.colOccC2" ,   "beta.colOccS1" ,   "beta.colOccS2",
               "beta.survAloneC1", "beta.survAloneC2", "beta.survAloneS1", "beta.survAloneS2",
               "beta.survBothC1",  "beta.survBothC2" , "beta.survBothS1",  "beta.survBothS2" )


par <- par %>%
  ungroup() %>%
  filter(param %in% par_index) %>% droplevels()

par_index == unique(par$param)
#par <- left_join(par.med,par.q75) %>% mutate(upr= quantile(value, ), lwr = med-quant)


## create transition matrix
sites <- c("Empty", "Sp1_only", "Sp2_only", "Both")

trans <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(sites, sites) )


## storing transition states

## env parameters are centered, 
#so 0 is at mean tmp/rmi env variable is at mean
temp <- c(-1,1,0,0,0) 
rmii <- c(0,0,0,-1,1)



equil_state <- array(data = 0, dim = c(length(sites), 15, nrow(par), length(temp)),
                     dimnames = list(sites, NULL, par$param, NULL)) 

par_val <- array(data = 0, dim = c( 15, nrow(par), length(temp)),
                 dimnames = list(NULL, par$param, NULL)) 


for(k in 1:length(temp)){ 
  tmp <- temp[k]
  rmi <- rmii[k]
  
  for(i in 1:length(par_index)){
    print(i)
    alpha.colEmptyS <- (par[which(par$param == "alpha.colEmptyS"), "med"]) %>% unlist()
    beta.colEmptyS1 <- (par[which(par$param == "beta.colEmptyS1"), "med"]) %>% unlist()
    beta.colEmptyS2 <- (par[which(par$param == "beta.colEmptyS2"), "med"]) %>% unlist()
    alpha.colEmptyC <- (par[which(par$param == "alpha.colEmptyC"), "med"])%>% unlist()
    beta.colEmptyC1 <- (par[which(par$param == "beta.colEmptyC1"), "med"])%>% unlist()
    beta.colEmptyC2 <- (par[which(par$param == "beta.colEmptyC2"), "med"])%>% unlist()
    alpha.survAloneS <- (par[which(par$param == "alpha.survAloneS"), "med"]) %>% unlist()
    beta.survAloneS1 <- (par[which(par$param == "beta.survAloneS1"), "med"])%>% unlist()
    beta.survAloneS2 <- (par[which(par$param == "beta.survAloneS2"), "med"])%>% unlist()
    alpha.survAloneC <- (par[which(par$param == "alpha.survAloneC"), "med"])%>% unlist()
    beta.survAloneC1 <- (par[which(par$param == "beta.survAloneC1"), "med"])%>% unlist()
    beta.survAloneC2 <- (par[which(par$param == "beta.survAloneC2"), "med"])%>% unlist()
    alpha.colOccS <- (par[which(par$param == "alpha.colOccS"), "med"])%>% unlist()
    beta.colOccS1 <- (par[which(par$param == "beta.colOccS1"), "med"])%>% unlist()
    beta.colOccS2 <- (par[which(par$param == "beta.colOccS2"), "med"])%>% unlist()
    alpha.colOccC <- (par[which(par$param == "alpha.colOccC"), "med"])%>% unlist()
    beta.colOccC1 <- (par[which(par$param == "beta.colOccC1"), "med"])%>% unlist()
    beta.colOccC2 <- (par[which(par$param == "beta.colOccC2"), "med"])%>% unlist()
    alpha.survBothS <- (par[which(par$param == "alpha.survBothS"), "med"])%>% unlist()
    beta.survBothS1 <- (par[which(par$param == "beta.survBothS1"), "med"])%>% unlist()
    beta.survBothS2 <- (par[which(par$param == "beta.survBothS2"), "med"])%>% unlist()
    alpha.survBothC <- (par[which(par$param == "alpha.survBothC"), "med"])%>% unlist()
    beta.survBothC1 <- (par[which(par$param == "beta.survBothC1"), "med"])%>% unlist()
    beta.survBothC2 <- (par[which(par$param == "beta.survBothC2"), "med"])%>% unlist()
    
    
    lwr <- par$lwr[i]
    upr <- par$upr[i]
    
    val <- seq(lwr, upr, length.out = 15)
    
    
    
    for(j in 1:length(val)){
      
      assign(par_index[i], val[j] )
      ## fill in transition matrix
      
      sp1_col <- plogis(alpha.colEmptyS + beta.colEmptyS1*tmp + beta.colEmptyS2*rmi)
      print(sp1_col)
      sp2_col <- plogis(alpha.colEmptyC + beta.colEmptyC1*tmp + beta.colEmptyC2*rmi)
      print(sp2_col)
      sp1_surv <- plogis(alpha.survAloneS + beta.survAloneS1*tmp + beta.survAloneS2*rmi)
      print(sp1_surv)
      sp2_surv <- plogis(alpha.survAloneC + beta.survAloneC1*tmp + beta.survAloneC2*rmi)
      print(sp2_surv)
      
      cond_sp1_col <- plogis(alpha.colOccS + beta.colOccS1*tmp + beta.colOccS2*rmi)
      print(cond_sp1_col)
      cond_sp2_col <- plogis(alpha.colOccC + beta.colOccC1*tmp + beta.colOccC2*rmi)
      print(cond_sp2_col)
      cond_sp1_surv <- plogis(alpha.survBothS + beta.survBothS1*tmp + beta.survBothS2*rmi)
      print(cond_sp1_surv)
      cond_sp2_surv <-plogis(alpha.survBothC + beta.survBothC1*tmp + beta.survBothC2*rmi) 
      print(cond_sp2_surv)
      
      trans[1,1] <- 1- (sp1_col * (1- sp2_col)) - (sp2_col * (1- sp1_col)) - (sp2_col * sp1_col)
      trans[2,1] <- sp1_col * (1- sp2_col)
      trans[3,1] <- sp2_col * (1- sp1_col)
      trans[4,1] <- sp2_col *  sp1_col
      
      
      trans[1,2] <- 1 - (sp1_surv * (1-cond_sp2_col)) - ((1-sp1_surv) * cond_sp2_col) - (sp1_surv * cond_sp2_col)
      trans[2,2] <- sp1_surv * (1-cond_sp2_col)
      trans[3,2] <- (1-sp1_surv) * cond_sp2_col
      trans[4,2] <- sp1_surv * cond_sp2_col
      
      trans[1,3] <- 1 - (sp2_surv * (1-cond_sp1_col)) - ((1-sp2_surv) * cond_sp1_col) - (sp2_surv * cond_sp1_col)
      trans[2,3] <- sp2_surv * (1-cond_sp1_col)
      trans[3,3] <- (1-sp2_surv) * cond_sp1_col
      trans[4,3] <- sp2_surv * cond_sp1_col
      
      trans[1,4] <- 1 - (cond_sp1_surv * (1- cond_sp2_surv)) - (cond_sp2_surv * (1- cond_sp1_surv)) - (cond_sp1_surv * cond_sp2_surv)
      trans[2,4] <- cond_sp1_surv * (1- cond_sp2_surv)
      trans[3,4] <- cond_sp2_surv * (1- cond_sp1_surv)
      trans[4,4] <- cond_sp1_surv * cond_sp2_surv
      
      ## calculate stable state distribution   
      
      stable_state <- trans %*% (eigen(trans)$vectors[, 1] / sum(eigen(trans)$vectors[, 1]))
      par_val[j,i,k] <- val[j] ## store parameter value
      equil_state[ ,j,i,k] <- as.numeric(stable_state[,1]) 
      print(j)
      
    } #close j (values)
  } #close i (parameters)
} #close k (temp)

## calculate sensitivity

## change in equilibrium
zz <- reshape2::melt(equil_state) %>% 
  rename(state = Var1, param = Var3, value= value, iter = Var2, env = Var4)

tt <- zz %>% 
  group_by(state, param, env) %>% 
  arrange(iter) %>%
  mutate(change =c(NA, diff(value, lag = 1)))


pp <- reshape2::melt(par_val) %>% rename(iter=Var1, param= Var2, par_val =value, env = Var3)

rr <- pp %>% 
  group_by(param, env) %>% 
  arrange(iter) %>%
  mutate(par_diff =c(NA, diff(par_val, lag = 1)))

qq <- left_join(tt, rr)
qq <- qq %>% mutate(sens = change/par_diff)


qq2 <- qq %>% group_by(param, env) %>%
  mutate(elasc = abs(sens))

qq1 <- qq2 %>% group_by(param) %>% filter(env == 1, elasc > 0.01) %>% droplevels()
(me1 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic  \n mean RMI -1 std temperature 
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 2, elasc > 0.01) %>% droplevels()
(me2 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic  \n mean RMI +1 std temperature
  \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 3, elasc > 0.01) %>% droplevels()
(me3 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic  \n mean RMI mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 4, elasc > 0.01) %>% droplevels()
(me4 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic \n -1 std RMI mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none") )



qq1 <- qq2 %>% group_by(param) %>% filter(env == 5, elasc > 0.01) %>% droplevels()
(me5 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic-Biotic  \n +1 std RMI mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none") )




#nd <- qq %>% filter(state == "Sp2_only" | state == "Both")
#nd <- nd%>% group_by(iter,param, env) %>% mutate(tot_host = sum(value))

nd1 <- qq %>% filter(state == "Sp1_only" | state == "Both")%>% 
  group_by(iter,param, env) %>%
  mutate(tot_host = sum(value))

nd1 <- nd1 %>% group_by(param, env) %>%
  mutate(occ_change = max(tot_host) - min(tot_host))

nd2 <- nd1 %>% group_by(param) %>% filter(env == 1, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 1, param %in% nd2$param)
(mo1 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy") 
)


q <- gridExtra::arrangeGrob(me1, mo1,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 2, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 2, param %in% nd2$param)
(mo2 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy")
)

q <- gridExtra::arrangeGrob(me2, mo2,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)

nd2 <- nd1 %>% group_by(param) %>% filter(env == 3, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 3, param %in% nd2$param)
(mo3 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))


q <- gridExtra::arrangeGrob(me3, mo3,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 4, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 4, param %in% nd2$param)
(mo4 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))


q <- gridExtra::arrangeGrob(me4, mo4,
                            layout_matrix = matrix(c(1, 1,1,1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 5, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 5, param %in% nd2$param)
(mo5 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))

q <- gridExtra::arrangeGrob(me5, mo5,
                            layout_matrix = matrix(c(1, 1,1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)





## abiotic moss
out <- readRDS("params.RDS")
out <- out %>% filter(model == "Abiotic_moss" , year == 15)



#par.med <- out %>% group_by(param) %>% 
# summarize(enframe(quantile(value, probs = c(0.5, 0.875, 0.125)), "quantile"))

par.med <- out %>% group_by(param) %>% 
  summarize(enframe(quantile(value, probs = c(0.5, 0.975, 0.025)), "quantile"))

par <- par.med %>% pivot_wider(names_from = quantile, values_from = value)
colnames(par) <- c("param", "med", "upr", "lwr")
#par.q75 <- reshape2::melt(out$q97.5[1:32]) %>% rename(param = L1, quant = value)
rm(out)
gc()
par_index <- c("alpha.colC",  "alpha.colS"  ,  "alpha.survC",
               "alpha.survS" ,"beta.colC1" , "beta.colC3",  "beta.colS1" , "beta.colS3" , 
               "beta.survC1" ,"beta.survC3",
               "beta.survS1", "beta.survS3" )


par <- par %>%
  ungroup() %>%
  filter(param %in% par_index) %>% droplevels()

par_index == unique(par$param)
#par <- left_join(par.med,par.q75) %>% mutate(upr= quantile(value, ), lwr = med-quant)


## create transition matrix
sites <- c("Empty", "Sp1_only", "Sp2_only", "Both")

trans <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(sites, sites) )


## storing transition states

## env parameters are centered, 
#so 0 is at mean tmp/moss, -1 is 1 std dv below and 1 is 1 std above
## combos assume other env variable is at mean
temp <- c(-1,1,0,0,0) 
moss <- c(0,0,0,-1,1)



equil_state <- array(data = 0, dim = c(length(sites), 15, nrow(par), length(temp)),
                     dimnames = list(sites, NULL, par$param, NULL)) 

par_val <- array(data = 0, dim = c( 15, nrow(par), length(temp)),
                 dimnames = list(NULL, par$param, NULL)) 


for(k in 1:length(temp)){ 
  tmp <- temp[k]
  mss <- moss[k]
  
  for(i in 1:length(par_index)){
    print(i)
    alpha.colC <- (par[which(par$param == "alpha.colC"), "med"]) %>% unlist()
    alpha.colS <- (par[which(par$param == "alpha.colS"), "med"]) %>% unlist()
    alpha.survC <-  (par[which(par$param == "alpha.survC"), "med"]) %>% unlist()
    alpha.survS <- (par[which(par$param == "alpha.survS"), "med"]) %>% unlist()
    beta.colC1 <- (par[which(par$param == "beta.colC1"), "med"]) %>% unlist()
    beta.colC3 <- (par[which(par$param == "beta.colC3"), "med"]) %>% unlist()
    beta.colS1 <- (par[which(par$param == "beta.colS1"), "med"]) %>% unlist()
    beta.colS3 <- (par[which(par$param == "beta.colS3"), "med"]) %>% unlist()
    beta.survC1 <- (par[which(par$param == "beta.survC1"), "med"]) %>% unlist()
    beta.survC3 <- (par[which(par$param == "beta.survC3"), "med"]) %>% unlist()
    beta.survS1 <- (par[which(par$param == "beta.survS1"), "med"]) %>% unlist()
    beta.survS3 <- (par[which(par$param == "beta.survS3"), "med"]) %>% unlist()
    
    
    lwr <- par$lwr[i]
    upr <- par$upr[i]
    
    val <- seq(lwr, upr, length.out = 15)
    
    
    
    for(j in 1:length(val)){
      
      assign(par_index[i], val[j] )
      ## fill in transition matrix
      
      sp1_col <- plogis(alpha.colS + beta.colS1*tmp + beta.colS3*mss)
      print(sp1_col)
      sp2_col <- plogis(alpha.colC + beta.colC1*tmp + beta.colC3*mss)
      print(sp2_col)
      sp1_surv <- plogis(alpha.survS + beta.survS1*tmp + beta.survS3*mss)
      print(sp1_surv)
      sp2_surv <- plogis(alpha.survC + beta.survC1*tmp + beta.survC3*mss)
      print(sp2_surv)
      
      
      trans[1,1] <- 1- (sp1_col * (1- sp2_col)) - (sp2_col * (1- sp1_col)) - (sp2_col * sp1_col)
      trans[2,1] <- sp1_col * (1- sp2_col)
      trans[3,1] <- sp2_col * (1- sp1_col)
      trans[4,1] <- sp2_col *  sp1_col
      
      
      trans[1,2] <- 1 - (sp1_surv * (1-sp2_col)) - ((1-sp1_surv) * sp2_col) - (sp1_surv * sp2_col)
      trans[2,2] <- sp1_surv * (1-sp2_col)
      trans[3,2] <- (1-sp1_surv) * sp2_col
      trans[4,2] <- sp1_surv * sp2_col
      
      trans[1,3] <- 1 - (sp2_surv * (1-sp1_col)) - ((1-sp2_surv) * sp1_col) - (sp2_surv * sp1_col)
      trans[2,3] <- sp2_surv * (1-sp1_col)
      trans[3,3] <- (1-sp2_surv) * sp1_col
      trans[4,3] <- sp2_surv * sp1_col
      
      trans[1,4] <- 1 - (sp1_surv * (1- sp2_surv)) - (sp2_surv * (1- sp1_surv)) - (sp1_surv * sp2_surv)
      trans[2,4] <- sp1_surv * (1- sp2_surv)
      trans[3,4] <- sp2_surv * (1- sp1_surv)
      trans[4,4] <- sp1_surv * sp2_surv
      
      ## calculate stable state distribution   
      
      stable_state <- trans %*% (eigen(trans)$vectors[, 1] / sum(eigen(trans)$vectors[, 1]))
      par_val[j,i,k] <- val[j] ## store parameter value
      equil_state[ ,j,i,k] <- as.numeric(stable_state[,1]) 
      print(j)
      
    } #close j (values)
  } #close i (parameters)
} #close k (temp)

## calculate sensitivity

## change in equilibrium
zz <- reshape2::melt(equil_state) %>% 
  rename(state = Var1, param = Var3, value= value, iter = Var2, env = Var4)

tt <- zz %>% 
  group_by(state, param, env) %>% 
  arrange(iter) %>%
  mutate(change =c(NA, diff(value, lag = 1)))


pp <- reshape2::melt(par_val) %>% rename(iter=Var1, param= Var2, par_val =value, env = Var3)

rr <- pp %>% 
  group_by(param, env) %>% 
  arrange(iter) %>%
  mutate(par_diff =c(NA, diff(par_val, lag = 1)))

qq <- left_join(tt, rr)
qq <- qq %>% mutate(sens = change/par_diff)



nd1 <- qq %>% filter(state == "Sp1_only" | state == "Both")%>% 
  group_by(iter,param, env) %>%
  mutate(tot_host = sum(value))

qq2 <- qq %>% group_by(param, env) %>%
  mutate(elasc = abs(sens))

qq1 <- qq2 %>% group_by(param) %>% filter(env == 1, elasc > 0.01) %>% droplevels()
(me1 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic (moss)  \n mean moss cover -1 std temperature 
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 2, elasc > 0.01) %>% droplevels()
(me2 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic (moss)  \n mean moss cover +1 std temperature
  \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 3, elasc > 0.01) %>% droplevels()
(me3 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic (moss)  \n mean moss cover mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none"))


qq1 <- qq2 %>% group_by(param) %>% filter(env == 4, elasc > 0.01) %>% droplevels()
(me4 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic (moss) \n -1 std moss cover mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none") )



qq1 <- qq2 %>% group_by(param) %>% filter(env == 5, elasc > 0.01) %>% droplevels()
(me5 <- ggplot( qq1, aes(par_val, sens)) + 
    geom_point(aes(color=state)) + 
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) +
    ggtitle("Abiotic (moss)  \n +1 std moss cover mean temperature
          \n \n A. Elasticity") +
    xlab("Parameter Value") +
    ylab("Elasticity") + guides(color = "none") )




#nd <- qq %>% filter(state == "Sp2_only" | state == "Both")
#nd <- nd%>% group_by(iter,param, env) %>% mutate(tot_host = sum(value))

nd1 <- qq %>% filter(state == "Sp1_only" | state == "Both")%>% 
  group_by(iter,param, env) %>%
  mutate(tot_host = sum(value))

nd1 <- nd1 %>% group_by(param, env) %>%
  mutate(occ_change = max(tot_host) - min(tot_host))

nd2 <- nd1 %>% group_by(param) %>% filter(env == 1, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 1, param %in% nd2$param)
(mo1 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy") 
)


q <- gridExtra::arrangeGrob(me1, mo1,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 2, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 2, param %in% nd2$param)
(mo2 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy")
)

q <- gridExtra::arrangeGrob(me2, mo2,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)

nd2 <- nd1 %>% group_by(param) %>% filter(env == 3, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 3, param %in% nd2$param)
(mo3 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))


q <- gridExtra::arrangeGrob(me3, mo3,
                            layout_matrix = matrix(c(1, 1, 1, 1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 4, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 4, param %in% nd2$param)
(mo4 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))


q <- gridExtra::arrangeGrob(me4, mo4,
                            layout_matrix = matrix(c(1, 1,2,2), byrow = T) )
gridExtra::grid.arrange(q)



nd2 <- nd1 %>% group_by(param) %>% filter(env == 5, occ_change > 0.05) %>% droplevels()
qq1 <- qq %>% group_by(param) %>% filter(env == 5, param %in% nd2$param)
(mo5 <- ggplot(qq1, aes(par_val, value)) + 
    geom_point(aes(color=state)) + 
    geom_point(data=nd2, aes(par_val, tot_host), shape = 17) +
    facet_wrap(~param, scales="free_x") +
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    theme(panel.spacing.x = unit(0.5, "lines")) + 
    ggtitle("B. P. shenandoah occupancy") +
    xlab("Parameter Values") + ylab("Occupancy"))

q <- gridExtra::arrangeGrob(me5, mo5,
                            layout_matrix = matrix(c(1, 1,1, 2,2), byrow = T) )
gridExtra::grid.arrange(q)

