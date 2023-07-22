#################################################
####### CREATE PROJECT FIGURES ##################
#################################################
library(tidyverse)


### regardless of leaver rate, s and beta are the only parameters that matter for 
#determining dominance

#stability only depends on whether s + b <= 1


#pair this plot with actual runs with different leaver rates
#GET OVER TIME SYSTEM BEHAVIOR
zoox_trajectory <- function(maxt, P, l1, l2, d, b){
  t = 1 #initial time
  while(t <= maxt){
    p1_t <- as.numeric(P[t,1])
    p2_t <- as.numeric(P[t,2])
    p1_tp1 <- round(p1_t*(1+(l1*(d + b - 1)))/(1+(((l1*p1_t) + (l2*p2_t))*(d + b - 1))), 5)
    p2_tp1 <- round(p2_t*(1+(l2*(d + b - 1)))/(1+(((l1*p1_t) + (l2*p2_t))*(d + b - 1))), 5)
    gr <- (1+(((l1*p1_t) + (l2*p2_t))*(d + b - 1)))
    P <- rbind(P, c(p1_tp1, p2_tp1, gr))
    if(sum(abs(P[t+1,1:2] - P[t,1:2])) < 1e-4){
      break
    }
    if(sum(abs(P[t+1,1:2] - P[t,1:2])) >= 1e-4){
      t = t+1
    }
    
  }
  found <- ifelse(t == (maxt+1), "notfound", "found")
  ss <- c(t, found, P[t,])
  names(ss) <- c("time_to_steady", "found", "p1", "p2", "gr")
  colnames(P) = c("p_stay", "p_leave", "gr_wrong")
  return(P)
}

start = .5
maxt = 500
P <- matrix(0, ncol = 3, nrow = 1)
#rows of P must sum to 1
P[1,1] = start #initial value for type 1
P[1,2] = 1-start #initial value for type 2
s = .5
#b = .5

l1 = .1
out_lowsb = data.frame()

for(b in c(0.6, 0.4)){
for(l2 in seq(.2, .9, .1)){
  temp = data.frame(zoox_trajectory(maxt, P, l1, l2, s, b))
  temp$l1 = l1
  temp$l2 = l2
  temp$s = s
  temp$b = b
  temp$time = 1:nrow(temp)
  out_lowsb = rbind(out_lowsb, temp)
}
}

sb_compare <- out_lowsb %>% ggplot(aes(x = time, y = p_leave, col = l2, group = interaction(l2, b))) + geom_line(lwd = 1.5) + theme_classic() + 
   scale_color_gradient(low = "#D6E4D7", high = "#263A28")+ ylab("Proportion Biomass Occupied \nby Unfaithful Symbiont") + ylim(c(0,1)) + 
  geom_hline(yintercept = 0.5, lty = 'dashed', col = 'grey') + xlab("Time") +   theme(text = element_text(size = 10)) #+scale_color_manual(values = c("#D6E4D7", "#6D9F71","#263A28")) + 

  sb_compare
#ggsave("figures/dominance_compare.png", sb_compare, width = 3.5, height = 3, units = "in")

#gradient scale for leaver color low = "#D6E4D7", high = "#263A28"
svec = seq(.01, 1, .01)
bvec = seq(.01, 1, .01)
viz_sb = data.frame(expand.grid(svec, bvec))
colnames(viz_sb) = c("s", "beta")
viz_sb = viz_sb %>% mutate(sum = s+beta) %>% filter(sum != 1)
#View(viz_sb)
viz_sb = viz_sb %>% mutate(g_t_1 = ifelse(sum >1, 1, 0))

concept_sb <- viz_sb %>% ggplot(aes(x = s, y = beta, fill = as.factor(g_t_1))) + geom_tile() + theme_classic() + 
  scale_fill_manual(values = c("#40476D", "#6D9F71")) + geom_abline(slope = -1, intercept = 1, col = 'white', lwd = 2) + 
  theme(legend.position="none") + theme(text = element_text(size = 10))
concept_sb

#ggsave("figures/conceptual_dominance.png", concept_sb, width = 3, height = 3, units = "in")
  

#######################################
### Think about growth rate ###########
######################################
#represent with biomass equation
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

maxt = 100
B <- matrix(0, ncol = 2, nrow = 1)
#rows of P must sum to 1
B[1,1] = 100 #initial value for type 1
B[1,2] = 100 #initial value for type 2

#l1 = 
l1 = .1
vs = .9
vl = .72
bc = .36
s = vl/vs
beta = bc/vs
s
beta

l_high = data.frame(biomass_t(maxt, B, l1, l2 = .6, vs, vl, bc))
l_high$lr = .6
l_high$time = 1:nrow(l_high)
l_low = data.frame(biomass_t(maxt, B, l1, l2 = .4, vs, vl, bc))
l_low$lr = .4
l_low$time = 1:nrow(l_low)
l_compare = rbind(l_high, l_low)

l_summary = l_compare %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)

#low = "#D6E4D7", high = "#263A28"
library(scales)
cc <- scales::seq_gradient_pal("#D6E4D7", "#263A28", "Lab")(seq(0,1,length.out=10))

gr_dis = l_summary %>% ggplot(aes(x = time, y = total_biomass, group = lr, col = as.factor(lr))) + geom_line(lwd = 2) + 
  scale_color_manual(values = rev(cc[c(9,2)])) + theme_classic() + theme(legend.position = "none") +
  ylab("System Biomass") + xlab("Time")+ ylim(c(0, max(l_summary))) + theme(text = element_text(size = 10))

gr_dis
#ggsave("figures/gr_disparity_leavers_win.png", gr_dis, width = 2, height = 2, units = "in")

gr_dis_prop = l_summary %>% ggplot(aes(x = time, y = prop_unfaith, group = lr, col = as.factor(lr))) + geom_line(lwd = 2) + 
  scale_color_manual(values = rev(cc[c(9,2)])) + theme_classic() + theme(legend.position = "none") +
  ylab("p_u") + xlab("Time")+ 
  ylim(c(0,1)) + geom_hline(yintercept = 0.5, col = "gray", lty = "dashed") + theme(text = element_text(size = 10))
gr_dis_prop
#ggsave("figures/gr_disparity_prop_leaver_wins.png", gr_dis_prop, width = 2, height = 2, units = "in")

#l_summary %>% ggplot(aes(x = time, y = log(total_biomass), group = lr)) + geom_line()

gr_calc <- function(vs, l, s, beta){
  gr_eq = vs*(1+l*(s + beta -1))
  return(gr_eq)
}
gr_calc(vs, l2, s, beta)


##### do it for stayer dominance
maxt = 100
B <- matrix(0, ncol = 2, nrow = 1)
#rows of P must sum to 1
B[1,1] = 100 #initial value for type 1
B[1,2] = 100 #initial value for type 2

#l1 = 
#l1 = .1
l2 = .9 #must be higher than l1
vs = 1.1
vl = .66
bc = .22
s = vl/vs
beta = bc/vs
s
beta

l_high = data.frame(biomass_t(maxt, B, l1 = .6, l2, vs, vl, bc))
l_high$lr = .6
l_high$time = 1:nrow(l_high)
l_low = data.frame(biomass_t(maxt, B, l1= .4, l2, vs, vl, bc))
l_low$lr = .4
l_low$time = 1:nrow(l_low)
l_compare = rbind(l_high, l_low)

l_summary = l_compare %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2)

#low = "#D6E4D7", high = "#263A28"

gr_dis_stay = l_summary %>% ggplot(aes(x = time, y = total_biomass, group = lr, col = as.factor(lr))) + geom_line(lwd = 2) + 
  scale_color_manual(values = rev(cc[c(9,2)])) + theme_classic() + theme(legend.position = "none") +
  ylab("System Biomass") + xlab("Time") + theme(text = element_text(size = 10))
gr_dis_stay

#ggsave("figures/gr_disparity_stayers_win.png", gr_dis_stay, width = 2, height = 2, units = "in")

gr_dis_prop = l_summary %>% ggplot(aes(x = time, y = prop_unfaith, group = lr, col = as.factor(lr))) + geom_line(lwd = 2) + 
  scale_color_manual(values = rev(cc[c(9,2)])) + theme_classic() + theme(legend.position = "none") +
  ylab("p_u") + xlab("Time")+ 
  ylim(c(0,1)) + geom_hline(yintercept = 0.5, col = "gray", lty = "dashed") + theme(text = element_text(size = 10))

gr_dis_prop
#ggsave("figures/gr_disparity_prop_stayers_win.png", gr_dis_prop, width = 2, height = 2, units = "in")

#zoox_trajectory(500, P,  l1 = .8, l2 = .9, d = s, b = beta)
# something is wrong with the biomass calculation?

####### Create a more conceptual figure showing vo, s, b interaction 

#don't need to run it all the way 
gr_calc <- function(vs, l, s, beta){
  gr_eq = vs*(1+l*(s + beta -1))
  return(gr_eq)
}

l_dom = .6
conmat = data.frame()
for(s in seq(0, 1, .05)){
  for(b in seq(0, 1, .05))
    for(vo in seq(0, 2, .05)){
      conmat = rbind(conmat, c(s, b, vo, gr_calc(vo, l_dom, s, b)))
    }
}
colnames(conmat) = c("gamma", "beta", "v_o", "gr")
conmat$gplb = conmat$gamma + conmat$beta
conmat = conmat %>% mutate(livedie = ifelse(gr < 1, "die", "live")) %>% mutate(dominance = ifelse(gplb >1, "sus", "res")) 
congg = conmat %>%filter(gplb != 1) %>% ggplot(aes(x=gplb, y = v_o)) + geom_tile(aes(fill=interaction(dominance, livedie)), width =.1) + theme_classic() + coord_equal() + 
  geom_hline(yintercept = 1, col = 'white') + geom_vline(xintercept = 1, col = 'white', lwd = 2) +
  scale_fill_manual(values = c("#A7ADC6","#BECDB7", "#3B3F61", "#759467")) + theme(legend.position="none", 
                                                                                     text = element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_blank())
congg
#ggsave("figures/gr_distparity_generalized.png", congg, width = 3, height = 3, units = "in")


#just plot the curve
grow_steady = function(l, gamma, beta){
  g = 1/(1+(l*(gamma+beta-1)))
  return(g)
}

g.df = data.frame()
for(g in seq(0, 1, .1)){
  for(b in seq(0, 1, .1)){
    for(ldom in seq(0, 1, .1)){
      g.df = rbind(g.df, c(g, b, g+b, ldom, grow_steady(ldom, g, b)))
    }
  }
}
colnames(g.df) = c("gamma", "beta", "gplb", "ldom", "vsteady")
g.df =g.df %>% filter(gplb != 0)
cc <- scales::seq_gradient_pal("#D6E4D7", "#263A28", "Lab")(seq(0,1,length.out=10))

grdeplr <- g.df %>% filter(ldom < 1) %>% filter (ldom >0) %>% ggplot(aes(x = gplb, y = vsteady, col = as.factor(ldom), group = ldom)) +
  geom_line(lwd = 2) +
  theme_classic() + geom_vline(xintercept = 1) + geom_hline(yintercept = 1) + scale_color_manual(values = cc) +
  theme(legend.position = "none", text=element_text(size = 12), axis.title.x=element_blank(), axis.title.y=element_blank())  + scale_y_continuous(trans = "log10")
grdeplr
#ggsave("figures/gr_depend_lr.png", grdeplr, width = 3, height = 3, units = "in")


