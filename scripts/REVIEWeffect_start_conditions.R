################################################
######## SHOW EFFECTS OF STARTING PROP #########
################################################
library(tidyverse)

#OVERALL, CORRESPONDING TO FIG2

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



##### NEW FIG - 3d explore of first figure under different starting proportions
l2=.5
l1 = .1
explore_start = data.frame()

for(b in seq(.1, 1, .1)){
  for(s in seq(.1, 1, .1)){
  for(p1 in seq(0.1, .9, .1)){
    P <- matrix(0, ncol = 3, nrow = 1)
    #rows of P must sum to 1
    P[1,1] = p1 #initial value for type 1
    P[1,2] = 1-p1 #initial value for type 2
    temp = data.frame(zoox_trajectory(maxt, P, l1, l2, s, b))
    temp$l1 = l1
    temp$l2 = l2
    temp$s = s
    temp$b = b
    temp$start1 <- p1
    temp$start2 <- 1-p1
    temp$time = 1:nrow(temp)
    explore_start = rbind(explore_start, temp)
  }
  }
}

#taking forever, it's a big boy
library(plotly)
sum_stab <- explore_start %>% group_by(start1, s, b) %>% mutate(stead=max(time), spb = s+b) %>% filter(time==stead, spb != 1) %>% mutate(whichdom=ifelse(p_stay>p_leave, "pstay", "pleave"))
lower_stab <- sum_stab %>% filter(spb < 1)
stab3d <- plot_ly(data = sum_stab, x=~s, y=~b, z=~start1, type = "scatter3d", mode = "markers", alpha = .8,
        color = ~whichdom, colors = c( "#6D9F71", "#40476D")) %>%
  layout(scene = list(xaxis=list(title="Growth Ratio", range=c(.1, .8)),yaxis=list(title="Birth Ratio"),zaxis=list(title="Starting Proportion"),
                      camera=list(eye=list(x = 1.75, y = -1.75, z = 1.75))))
stab3d

stab_time<- plot_ly(data = lower_stab, x=~s, y=~b, z=~start1, color=~log10(time), type = "scatter3d", alpha = .8) %>%
  layout(scene = list(xaxis=list(title="Growth Ratio", range=c(.1, .8)),yaxis=list(title="Birth Ratio"),zaxis=list(title="Starting Proportion"),
                      camera=list(eye=list(x = 2.75, y = -0.75, z = 1.75)))) 

library(orca)
#if (!require("processx")) install.packages("processx")
#fig <- plot_ly(z = ~volcano) %>% add_surface()
#takes a while, error?
orca(stab3d, "figures/3dstability.svg", width = 3, height = 3)
orca(stab_time, "figures/3dstability_ot.svg", width = 3, height = 3)




###################################################################
#### explore Fig 3, invasion under different starting proportions 
###################################################################
############################################
##########  LR DEPENDENT NEW RECRUITMENT ANALYSIS #######
############################################
library(tidyverse)

#I think I just need to loop and plot this through different starting conditions
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
start_bio = 200
out_early <- data.frame()
for(i in 1:length(lr_options)){
  for(j in seq(.1, .9, .1)){
      maxt = 250
      B <- matrix(0, ncol = 2, nrow = 1)
      #rows of P must sum to 1
      B[1,1] = j*200 #initial value for type 1
      B[1,2] = (1-j)*200 #initial value for type 2
      
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
      l_recruit_fail$startp1 <- j
      out_early <- rbind(out_early, l_recruit_fail)
  }
      
}

out_late_nosave <- data.frame()
for(i in 1:length(lr_options)){
  for(j in seq(.1, .9, .1)){
  maxt = 250
  B <- matrix(0, ncol = 2, nrow = 1)
  l1 = lr_options[i]
  l2 = l_susceptible[1] #must be higher than l1
  #rows of P must sum to 1
  B[1,1] = 0 #initial value for type 1
  B[1,2] = out_early[out_early$time == 250 & out_early$lr1 == l1 & out_early$lr2 == l2 & out_early$startp1 == j,]$B2 #initial value for type 2
  
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
  
  l_recruit_fail = l_early %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2, startp1=j)
  out_late_nosave <- rbind(out_late_nosave, l_recruit_fail)
  
  }
}

out_late_save <- data.frame()
for(i in 1:length(lr_options)){
  for(j in seq(.1, .9, .1)){
  maxt = 250
  B <- matrix(0, ncol = 2, nrow = 1)
  l1 = lr_options[i]
  l2 = l_susceptible #must be higher than l1
  #rows of P must sum to 1
  B[1,1] = out_early[out_early$time == 250 & out_early$lr1 == l1 & out_early$lr2 == l2 & out_early$startp1==j,]$B1
  B[1,2] = out_early[out_early$time == 250 & out_early$lr1 == l1 & out_early$lr2 == l2 & out_early$startp1==j,]$B2 #initial value for type 2
  
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
  
  l_recruit_fail = l_early %>% mutate(prop_unfaith = B2/(B1+B2), total_biomass = B1+B2, startp1=j)
  out_late_save <- rbind(out_late_save, l_recruit_fail)
  
  }
}

#low = "#D6E4D7", high = "#263A28"
out_early$condition <- "noinvade"
out_late_nosave$condition <- "noinvade"
out_late_save$condition <- out_late_save$lr1
out_all <- rbind(out_early, out_late_nosave, out_late_save)

#okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "blue")
#gradient scale for leaver color low = "#D6E4D7", high = "#263A28"
cc <- scales::seq_gradient_pal("#D6E4D7", "#263A28", "Lab")(seq(0,1,length.out=6))

pal <- c(cc, "gray")
gr_recruit_fail = out_all %>% ggplot(aes(x = time, y = log(total_biomass), col = interaction(condition), group = interaction(condition))) +
  geom_line(lwd = 2) + 
  theme_classic() + theme(legend.position = "bottom", text = element_text(size = 10)) +
  labs(color = "lr") +
  ylab("log(Coral Biomass)") + xlab("Time")  + scale_color_manual(values =pal) + 
  geom_vline(xintercept = 250, col = 'black', lty = 'dashed') + xlim(c(200, 500)) + facet_wrap(~startp1)
gr_recruit_fail
#ggsave("figures/lr_dependent_recruitment.png", width = 6.5, height = 7, units = "in", dpi = 500)

#doesn't really change anything
#code is working though:
out_early %>% filter(time==250) %>% ggplot(aes(x=startp1, y=log(total_biomass))) + geom_point(aes(col=lr1))



###################################################################
#### explore Fig 4, gr under different starting proportions 
###################################################################

#I don't know if there's a way to revise this for starting condition

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
