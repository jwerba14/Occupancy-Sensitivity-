library(tidyverse)

##### invasive species scenarios ######
sites <- c("Empty", "Sp1_only", "Sp2_only", "Both")

## rare native species 
lab <- c("hS", "mS", "mlS", "lS")
prob.lab <- c("surv", "col")

native_param <- matrix(data = c(0.9,0.5,0.2,0.005,
                               0.005, 0.2,0.5,0.9), byrow = F,ncol = 2, 
                       dimnames = list(lab, prob.lab))   ## from high to low survival;

## native species conditional survival and colonization 5%, 40%, 9%
nat_par <- as.data.frame(native_param) %>%
  mutate(cond_surv.1 = surv-(surv*0.05), cond_col.1 = col- (col*0.05), 
         cond_surv.2 = surv-(surv*0.4), cond_col.2 = col-(col*0.4),
         cond_surv.3 = surv-(surv*0.9), cond_col.3 = col-(col*0.9)) %>%
  pivot_longer(-c("surv","col"), names_to= "cond_param", values_to = "conditional")

nat_par <- nat_par %>% separate(cond_param, sep =c("[.]"), into = c("cond_param", "scenario") )
  
## create transition matrix
sites <- c("Empty", "Sp1_only", "Sp2_only", "Both")

trans <- matrix(data = NA, nrow = 4, ncol = 4, dimnames = list(sites, sites) )

tr.names <- c("hS", "mhS", "mlS", "lS")

cond.names <- c("l", "m", "h")

## 
out_params <- data.frame(
  col = seq(0.9, 0.05, length.out = 50),
  surv = seq(0.05, 0.9, length.out = 50),
  cond_surv = seq(0.05, 0.9, length.out = 50),
  cond_col = seq(0.9, 0.05, length.out = 50)
)

## storing transition states

equil_state <- array(data = 0, dim = c(4, length(unique(nat_par$surv)), 3, nrow(out_params)),
                     dimnames = list(sites,tr.names , cond.names, NULL)) 

## fast invader relatively low survival

for(k in 1:nrow(out_params)){
  sp2_col = out_params$col[k]
  sp2_surv = out_params$surv[k]
  cond_sp2_col = out_params$cond_col[k]
  cond_sp2_surv = out_params$cond_surv[k]
for(i in 1:3) {   
  ## loop over 3 invasive scenarios
  temp <- nat_par %>% filter %>% filter(grepl(i, scenario)) %>%
    pivot_wider(., names_from = "cond_param", values_from = "conditional") 
  
  for(j in 1:length(unique(nat_par$surv))){  ## loop over native species param
    sp1_col = temp$col[j]
    sp1_surv = temp$surv[j]
    cond_sp1_col = temp$cond_col[j]
    cond_sp1_surv = temp$cond_surv[j]
   
    
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
      
      equil_state[ ,j,i,k] <- as.numeric(stable_state[,1])
  } #close j
 } # close i
} # close k



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



qq <- reshape2::melt(equil_state) %>% 
  rename(Site = Var1, k = Var4)  %>%
  mutate(surv = plyr::mapvalues(k, from = seq(1,50), to = out_params$surv),
          col = plyr::mapvalues(k, from = seq(1,50), to = out_params$col)) %>%
  group_by(Site, Var2, Var3)%>% ## group by site, type of native species, effect on nat. spec
  arrange(k) %>%
  mutate(diff = c(NA, diff(value, lag = 1)),
         surv_c = c(NA, diff(surv, lag = 1)),
         col_c = c(NA, diff(col, lag = 1))) %>%
   mutate(sens1 = diff/surv_c, sens2 = diff/col_c)


lab.var2 <- c("High Survival (0.9)", "Mid Survival (0.5)", 
              "Mid-Low Survival (0.2)", "Low Survival (0.005)")

names(lab.var2) <- c("hS", "mhS", "mlS", "lS")

lab.var3 <- c("5% decrease \n due to invasion", "40% decrease \n due to invasion",
              "90% decrease \n due to invasion")

names(lab.var3) <- c("l", "m", "h")

ggplot(qq, aes(surv, sens1)) + geom_point(aes(color = Site)) +
  facet_grid(Var2~Var3, labeller = labeller(Var2 = lab.var2, 
                                                     Var3 = lab.var3)) +
  xlab("Invader Survival Probability") + ylab("Elasticity") +
  scale_color_discrete(labels = 
                         c("Empty", "Native Species Only", "Invasive Species Only", "Co-occurence"))

lab.var2a <- c("Low Colonization (0.005)", "Mid-Low Colonization (0.2)", 
              "Mid Colonization (0.5)", "High Colonization (0.9)")

names(lab.var2a) <- c("hS", "mhS", "mlS", "lS")

ggplot(qq, aes(col, sens2)) + geom_point(aes(color = Site)) +
  facet_grid(Var2~Var3, labeller = labeller(Var2 = lab.var2a, 
                                            Var3 = lab.var3)) +
  xlab("Invader Colonization Probability") + ylab("Elasticity") +
  scale_color_discrete(labels = 
                         c("Empty", "Native Species Only", "Invasive Species Only", "Co-occurence"))


## graph occupancy

## native occupancy pre-invasion
nat_equil <- array(data = 0, dim = c(2, nrow(native_param)),
                   dimnames = list(c("empty", "occupied"),tr.names)) 

trans <- matrix(data = NA, nrow = 2, ncol = 2)

for(i in 1:nrow(native_param)){
  
  col <- native_param[i,2]
  surv <- native_param[i,1]
  
  trans[1,1] <- 1-col
  trans[1,2] <- 1-surv
  trans[2,1] <- col
  trans[2,2] <- surv
  
  stable_state <- trans %*% (eigen(trans)$vectors[, 1] / sum(eigen(trans)$vectors[, 1]))
  
  nat_equil[ ,i] <- as.numeric(stable_state[,1])
  
}

nat_equil1 <- reshape2::melt(nat_equil) %>% filter(Var1 == "occupied")

zz <- reshape2::melt(equil_state) %>% 
  rename(Site = Var1, k = Var4)  %>%
  mutate(surv = plyr::mapvalues(k, from = seq(1,50), to = out_params$surv),
         col = plyr::mapvalues(k, from = seq(1,50), to = out_params$col)) %>%
  group_by(k, Var2, Var3) %>% ## group by site, type of native species, effect on nat. spec
 pivot_wider(names_from = Site, values_from = value) %>%
  mutate(tot_occ = (Sp1_only + Both))

ggplot(zz, aes(surv, tot_occ)) + geom_point() + 
  geom_hline(data = nat_equil1, aes(yintercept = value), linetype = "dashed", color = "red") +
  facet_grid(Var2~Var3, labeller = labeller(Var2 = lab.var2a, 
                                            Var3 = lab.var3)) +
  xlab("Invader Survival Probability") + ylab("Total Native Occupancy") 


