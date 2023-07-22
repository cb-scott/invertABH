############################################
##########  LR DEPENDENT NEW RECRUITMENT ANALYSIS #######
############################################
library(tidyverse)

biomass_t <- function(maxt, B, l1, l2, vs, vl, bc){
  t = 1 #initial time
  while(t <= maxt){
    b1_t <- as.numeric(B[t,1])
    b2_t <- as.numeric(B[t,2])
    b1_tp1 = ((1-l1)*b1_t*vs) + (l1*b1_t*vl) + l1*b1_t*bc
    b2_tp1 = ((1-l2)*b2_t*vs) + (l2*b2_t*vl) + l2*b2_t*bc
    B <- rbind(B, c(b1_tp1, b2_tp1))
    t = t+1
  }
  colnames(B) = c("B1", "B2")
  return(B)
}

gr_calc <- function(vs, l, s, beta){
  gr_eq = vs*(1+l*(s + beta -1))
  return(gr_eq)
}

grow_steady = function(l, gamma, beta){
  g = 1/(1+(l*(gamma+beta-1)))
  return(g)
}

#lr_options = data.frame(ls = c(.1, .8), lr = c(.2, .9))
lr_options = seq(.1, .6, .1)
l_susceptible = c(.7) 
out_early <- data.frame()
for(i in 1:length(lr_options)){
  maxt = 250
  B <- matrix(0, ncol = 2, nrow = 1)
  #rows of P must sum to 1
  B[1,1] = 100 #initial value for type 1
  B[1,2] = 100 #initial value for type 2
  
  #l1 = 
  l1 = lr_options[i]
  l2 = l_susceptible[1] #must be higher than l1
  vs = 1.05
  vl = .735
  bc = .42
  s = vl/vs
  beta = bc/vs
  s
  beta
  
  l_early = data.frame(biomass_t(maxt, B, l1 = l1, l2 = l2 , vs, vl, bc))
  l_early$lr2 = l2
  l_early$lr1 = l1
  l_early$time = 1:nrow(l_early)
  
  l_recruit_fail = l_early %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)
  out_early <- rbind(out_early, l_recruit_fail)
  
}

out_late_nosave <- data.frame()
for(i in 1:length(lr_options)){
  maxt = 250
  B <- matrix(0, ncol = 2, nrow = 1)
  l1 = lr_options[i]
  l2 = l_susceptible[1] #must be higher than l1
  #rows of P must sum to 1
  B[1,1] = 0 #initial value for type 1
  B[1,2] = out_early[out_early$time == 250 & out_early$lr1 == l1 & out_early$lr2 == l2,]$B2 #initial value for type 2
  
  #l1 = 

  vs = 1.05
  vl = .735
  bc = .42*.25 #recruitfailure!
  s = vl/vs
  beta = bc/vs
  s
  beta
  
  l_early = data.frame(biomass_t(maxt, B, l1 = l1, l2 = l2 , vs, vl, bc))
  l_early$lr2 = l2
  l_early$lr1 = l1
  l_early$time = 250+(1:nrow(l_early))
  
  l_recruit_fail = l_early %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)
  out_late_nosave <- rbind(out_late_nosave, l_recruit_fail)
  
}

out_late_save <- data.frame()
for(i in 1:length(lr_options)){
  maxt = 250
  B <- matrix(0, ncol = 2, nrow = 1)
  l1 = lr_options[i]
  l2 = l_susceptible #must be higher than l1
  #rows of P must sum to 1
  B[1,1] = out_early[out_early$time == 250 & out_early$lr1 == l1 & out_early$lr2 == l2,]$B1
  B[1,2] = out_early[out_early$time == 250 & out_early$lr1 == l1 & out_early$lr2 == l2,]$B2 #initial value for type 2
  
  #l1 = 
  
  vs = 1.05
  vl = .735
  bc = .42*.25
  s = vl/vs
  beta = bc/vs
  s
  beta
  
  l_early = data.frame(biomass_t(maxt, B, l1 = l1, l2 = l2 , vs, vl, bc))
  l_early$lr2 = l2
  l_early$lr1 = l1
  l_early$time = 250+(1:nrow(l_early))
  
  l_recruit_fail = l_early %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)
  out_late_save <- rbind(out_late_save, l_recruit_fail)
  
}

#low = "#D6E4D7", high = "#263A28"
out_early$condition <- "noinvade"
out_late_nosave$condition <- "noinvade"
out_late_save$condition <- out_late_save$lr1
out_all <- rbind(out_early, out_late_nosave, out_late_save)

okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "blue")
#gradient scale for leaver color low = "#D6E4D7", high = "#263A28"
cc <- scales::seq_gradient_pal("#D6E4D7", "#263A28", "Lab")(seq(0,1,length.out=6))

pal <- c(cc, "gray")
gr_recruit_fail = out_all %>% ggplot(aes(x = time, y = log(total_biomass), col = interaction(condition), group = interaction(condition))) +
  geom_line(lwd = 2) + 
  theme_classic() + theme(legend.position = "none", text = element_text(size = 10)) +
  ylab("log(Coral Biomass)") + xlab("Time")  + scale_color_manual(values =pal) + 
  geom_vline(xintercept = 250, col = 'black', lty = 'dashed') + xlim(c(200, 500))
gr_recruit_fail

ggsave("figures/lr_depend_recruit.png", gr_recruit_fail, width = 2.5, height = 2.5, units = "in", dpi = 500)



