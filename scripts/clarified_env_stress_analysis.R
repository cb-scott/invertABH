############################################
##########  NEW ENV STRESS ANALYSIS #######
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

maxt = 500 #burn in for 250 timesteps
B <- matrix(0, ncol = 2, nrow = 1)
#rows of P must sum to 1
B[1,1] = 100 #initial value for type 1
B[1,2] = 100 #initial value for type 2

#l1 = 
l1 = .1
l2 = .8 #must be higher than l1
vs = 1.1
vl = .1
bc = .2
s = vl/vs
beta = bc/vs
s
beta

l_grow = data.frame(biomass_t(maxt, B, l1 = l1, l2 = l2 , vs, vl, bc))
l_grow$lr2 = l2
l_grow$lr1 = l1
l_grow$time = 1:nrow(l_grow)
l_grow = l_grow %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)

#take last 100 of grow
tail_lgrow105 <- l_grow[(nrow(l_grow)-50):nrow(l_grow),]

#spend 50 years at lower growth
maxt = 50 #burn in for 250 timesteps
B <- matrix(0, ncol = 2, nrow = 1)
#rows of P must sum to 1
B[1,1] = tail_lgrow105[nrow(tail_lgrow105),"B1"] #initial value for type 1
B[1,2] = tail_lgrow105[nrow(tail_lgrow105),"B2"]#initial value for type 2

#l1 = 
l1 = .1
l2 = .6 #must be higher than l1
vs = .95
vl = .8*vs
bc = .2
s = vl/vs
beta = bc/vs
s
beta

l_grow_95 = data.frame(biomass_t(maxt, B, l1 = l1, l2 = l2 , vs, vl, bc))
l_grow_95$lr2 = l2
l_grow_95$lr1 = l1
l_grow_95$time = 1:nrow(l_grow_95)
l_grow_95 = l_grow_95 %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)

#Reduce further

maxt = 50 #burn in for 250 timesteps
B <- matrix(0, ncol = 2, nrow = 1)
#rows of P must sum to 1
B[1,1] = l_grow_95[nrow(l_grow_95),"B1"] #initial value for type 1
B[1,2] = l_grow_95[nrow(l_grow_95),"B2"]#initial value for type 2

l1 = .1
l2 = .6 #must be higher than l1
vs = .9
vl = .8*vs
bc = .2
s = vl/vs
beta = bc/vs
s
beta

l_grow_9 = data.frame(biomass_t(maxt, B, l1 = l1, l2 = l2 , vs, vl, bc))
l_grow_9$lr2 = l2
l_grow_9$lr1 = l1
l_grow_9$time = 1:nrow(l_grow_9)
l_grow_9 = l_grow_9 %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)

##create master df
lgrow <- rbind(tail_lgrow105, l_grow_95, l_grow_9)
lgrow$time <- 1:nrow(lgrow)

plot(lgrow$time, lgrow$prop_unfaith)
##### SHOW SYSTEM AT g = 1.05 vs g = .95 ######
#### then show repeated disturbance ########

#repeated disturbance prevents exclusive symbiont associations

#want s.s. solution at each  ge
go = 1.1
b = .2
beta = b/go
#l1 = .3
#l2 = .6
out <- matrix(0, ncol = 5, nrow = 1)
for(ge in seq(0.1, 1, .01)){
  for(l1 in c(.01, .49, .98)){
    for(l2 in c(.02, .5, .99)){
  gamma = ge/go
  lwin = ifelse(gamma + beta > 1, l2, l1)
  gr_ge <- gr_calc(go, lwin, gamma, beta)
  out <- rbind(out, c(ge, lwin,gr_ge, l1, l2))
  }
  }
}
out <- out[-1,]
out.df <- data.frame(out)
colnames(out.df) <- c("ge", "win", "gr", "l1", "l2")
out.df <- out.df %>% filter(l1 < l2)
out.df <- out.df %>% mutate(winner = ifelse(win == l1, "Resistant", "Susceptible"))
fill.df <- out.df  %>% pivot_wider(names_from = win, values_from = gr, names_prefix = "win")
#need to recolor this.
worsebleach = ggplot() +
  geom_line(data = out.df, aes(1-ge, gr, col = as.factor(win),group = win), lwd = 1.2) + 
  geom_hline(yintercept = 1) + geom_vline(xintercept = 1-.905)+
  scale_y_log10() + scale_x_log10()+ #geom_ribbon(data = subset(out.df, .905 > ge), aes(ymin=subset(out.df, .905 > ge & win == .01), ymax = subset(out.df, .905 > ge & win == .99)))+
   theme_classic() + ylab("Total Biomass Growth Factor-log axis") + xlab("Mortality of Bleached Biomass: log10(1-ge)") +
  #scale_color_manual(values = c("#4F547D", "#7B936B")) +
  theme(text=element_text(size = 10), legend.position = c(.25, .18), 
          legend.background = element_rect(colour = "black", size = .2), 
          legend.text = element_text(size = 8)) +
  guides(color = guide_legend(title = "Dominant Strain"))
worsebleach
ggsave("figures/worsebleach.png", worsebleach, width = 3.1, height = 3, units = "in")


