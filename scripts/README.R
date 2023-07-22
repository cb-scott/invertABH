#Clean Dataset Code From Source 
library(tidyverse)

######################################################################
######### PROCESS HUGHES DATA ########################################
######################################################################

bleach_DHW <- read_csv("data/Hughes2016-2017data_Bleach_DHW.csv", col_names = T)
recruit_hughes <- read_csv("data/RecruitmentData_Hughes.csv", col_names = T, skip = 1)
recruit_meta <- read_csv("data/ReefsSummary_Hughes.csv", col_names = T)


#hughes panels are 11cm/cm, so we want to calculate recruits per cm2 for everything 
#join the recruitment data and metadata by reef ID to combine lat/longs
recruit_join <- left_join(recruit_hughes, recruit_meta, by=c("ReefID" = "ReefID"))
#get data where reefIDs match
all_join <- left_join(recruit_join, bleach_DHW, by=c("ReefID" = "ReefID"))
#where is there no bleaching score for a reef?
missing_data <- all_join %>% filter(is.na(AerialScore))

bleach_meta <- bleach_DHW %>% dplyr::select(ReefID, ReefName, Longitude, Latitude) %>% distinct()


####### Reefs don't always have the same names and IDs - make sure that they do ########
#we're going to find the closest reef for all of these 
#get the unique missing reefs - this keeps us from needing to create a raster grid
NeedIDs <- unique(missing_data$ReefID) #reefs with recruitment data, without bleaching data
reefs <- c()
newid <- c()
for(reef in NeedIDs){
  latlon <- (missing_data %>% filter(ReefID == reef) %>% dplyr::select(Lat, Long))[1,]
  dist <- sqrt(((bleach_meta$Longitude - latlon$Long)^2) + ((bleach_meta$Latitude - latlon$Lat)^2))
  which_close <- which.min(dist)
  reefs <- c(reefs, reef)
  newid <- c(newid, as.character(bleach_meta[which_close,"ReefID"]))
}
missing_data_LUT <- data.frame(cbind(reefs, newid)) #assign data from nearest neighbor to each missing reef
colnames(missing_data_LUT) <- c("ReefIDOld", "ReefID")

#all reefs now have a valid ID from the closest reef with bleaching data 
#this is the ReefIDNew column - corresponds to the proper names from other data
recruit_newmeta <- left_join(recruit_meta, missing_data_LUT, by = c("ReefID" = "ReefIDOld")) %>% mutate(ReefIDNew = ifelse(is.na(ReefID.y), ReefID, ReefID.y)) %>%
  dplyr::select(-c(5:12), -ReefID.y)

#need yearly averages per cm sq for hughes
#sum over all tiles, then divide by cm/2 for all tiles 
sa_tile <- 11^2
recruit_sum_hughes <- recruit_hughes %>% group_by(ReefID, Year) %>% summarise(Total_sum = sum(Total), nobv = n()) %>% mutate(percm2 = Total_sum/(sa_tile*nobv))
recruit_newjoin <- left_join(recruit_sum_hughes, recruit_newmeta, by=c("ReefID" = "ReefID"))



############################################
######## COMBINE RECRUITMENT DATA #########
############################################
price_recruit <- read_csv("data/AllYear_Recruitment_Price_Renamed.csv", col_names = T)
#recruit_score is in /m2, need to convert to /cm2
price_recruit_mod <- price_recruit %>% mutate(recruit_cm2 = recruit_score*0.0001) #remarkably low

#check for duplicate observations in price and hughes (assume same year, same reef = same obs)
combine_recruitment <- full_join(price_recruit_mod, recruit_newjoin, by = c("ReefID" = "ReefIDNew", "Year" = "Year", "lat" = "Lat", "lon" = "Long", "recruit_cm2" = "percm2"))
combine_recruitment$ReefID.y <- NULL
#no overlapping data 

#how many unique reefs do we have data for? 
length(unique(combine_recruitment$ReefID))
#61 total data points

###########################################
######## COMBINE BLEACH DATA ###############
############################################

#first need to rename scores in hughes bleaching data
bleach_hughes <- bleach_only %>% mutate(BleachScore = case_when(
  #these are based on the percentile descriptions in both papers
  AerialScore == 0 ~ 0,
  AerialScore == 1 ~ 1,
  AerialScore == 2 ~ 2,
  AerialScore == 3 ~ 2.5,
  AerialScore == 4 ~3
))
bleach_hughes$AerialScore <- NULL

ay_bleach <- read_csv("data/AllYear_Bleaching_Renamed.csv", col_names = T) #comes from price?

ay_bleach <- ay_bleach %>% dplyr::select(ReefID, YEAR, bleach_score)
colnames(ay_bleach) <- c("ReefID", "Year", "BleachScore")

all_bleach <- rbind(bleach_hughes, ay_bleach)
#check for duplicates
dim(all_bleach)
dim(unique(all_bleach %>% dplyr::select(ReefID, Year))) #no overlapping data 

#assign each reefID a bleaching rank based on their over time average
bleach_rank <- all_bleach %>% group_by(ReefID) %>% summarise(b_o_t = mean(BleachScore), Nobvs_bleach = n(), var_o_t = var(BleachScore))
#
#We know that b_o_t is heavily influenced by Nobvs! Let's remove this, take the residual, then rank.
summary(lm(b_o_t~Nobvs_bleach, data = bleach_rank))
plot(bleach_rank$Nobvs_bleach, bleach_rank$b_o_t)
abline(b = -.07781, a = 2.02)
resid_bleach = residuals(lm(b_o_t~Nobvs_bleach, data = bleach_rank))
plot(bleach_rank$Nobvs_bleach, resid_bleach)
bleach_rank$residual_bleach = resid_bleach

bleach_rank <- bleach_rank[order(bleach_rank$residual_bleach),] 

getrank <- function(column){
  #column is an ordered vector
  column = as.vector(column)
  value = 0
  rank <- c(0)
  for(i in 1:(length(column)-1)){
    value = ifelse(column[i+1] > column[i], value + 1, value)
    rank <- c(rank, value)
  }
  return(rank)
}

#the rank of each reef for all time observations
bleach_rank$rank_bleach <- getrank(bleach_rank$residual_bleach)
#look at lon lat + bleach rank = there should be a pattern

bleach_rank_meta = left_join(bleach_rank, bleach_meta)
bleach_rank_meta %>% ggplot(aes(x=Longitude, y = Latitude, col = residual_bleach)) + geom_point() + scale_color_viridis_c()
#rem_sample_effect = residuals(lm(b_o_t ~ Nobvs_bleach, bleach_rank_meta))
#bleach_rank_meta$rank_bleach_rem_error <- rem_sample_effect
bleach_rank_meta %>% ggplot(aes(x=Longitude, y = Latitude, col = rank_bleach)) + geom_point() + scale_color_viridis_c()
bleach_rank_meta %>% ggplot(aes(x=Longitude, y = Latitude, col = Nobvs_bleach)) + geom_point() + scale_color_viridis_c()

#rank_bleach now represents the rank of the residuals of bleaching severity (to control for N observations)

########## MERGE THE RECRUITMENT DATA ####################

#combine_recruitment #reefs/years. calculate recruits_cm2

#Keep lon lat lookup
reef_lat_lon = bleach_rank_meta %>% select(ReefID, Longitude, Latitude) %>% group_by(ReefID) %>% summarise(lat = mean(Latitude), lon = mean(Longitude))

#combine it by years
recruitment_oy = combine_recruitment %>% group_by(ReefID) %>% summarise(mean_recruit_cm2 = mean(recruit_cm2), var_recruit = var(recruit_cm2), min_recruit = min(recruit_cm2), max_recruit = max(recruit_cm2), Nyear_obs = n())
#number of years observed not significant for this 
summary(lm(mean_recruit_cm2 ~ Nyear_obs, data = recruitment_oy))
#But variance in recruitment is? This is being driven by an outliers, don't worry about it.
#we could weigh these observations by their inverse variance? Later
summary(lm(mean_recruit_cm2 ~ var_recruit, data = recruitment_oy))

#how many totally overlapping observations are there?
#there ARE! bleaching obvservations for 61 reefs. 
recruit_bleach = inner_join(recruitment_oy, bleach_rank_meta)
recruit_bleach %>% ggplot(aes(x=mean_recruit_cm2, y = residual_bleach)) + geom_jitter() + geom_smooth(method = "lm")
summary(lm(residual_bleach ~mean_recruit_cm2, data = recruit_bleach))
lm_br = lm(residual_bleach ~  log(mean_recruit_cm2), data = recruit_bleach)
lm_br_sum = summary(lm_br)
lm_br_sum
plot(log(recruit_bleach$mean_recruit_cm2), recruit_bleach$residual_bleach )
abline(b = .2099, a = .481, col = 'red', lwd = 2)
plot(recruit_bleach$residual_bleach~log(recruit_bleach$mean_recruit_cm2))
#wait.... latitude is eating up my correlation? Leave it out.
cbind(lm_br_sum$coefficients[,1]-lm_br_sum$coefficients[,2]*qnorm(0.975), lm_br_sum$coefficients[,1]+lm_br_sum$coefficients[,2]*qnorm(0.975))

recruit_bleach %>% ggplot(aes(x = Longitude, y = Latitude, col = residual_bleach)) + geom_point(size = 3, alpha = .8) + scale_color_viridis_c() + 
  theme_classic() 

recruit_bleach %>% ggplot(aes(x = Longitude, y = Latitude, col = log(mean_recruit_cm2))) + geom_point(size = 3, alpha = .8) + scale_color_viridis_c() + 
  theme_classic()

lm_recruitment <- recruit_bleach %>% ggplot(aes(x = log(mean_recruit_cm2), y = residual_bleach)) + geom_point(size = 3, alpha = .8) + scale_color_viridis_c() + 
  theme_classic() + geom_smooth(method = "lm") + xlab("Log(Mean Recruitment/cm^2)") + ylab("Reef-Specific Bleaching Tendency") + theme(text = element_text (size = 10)) 
summary(lm(recruit_bleach$residual_bleach~recruit_bleach$mean_recruit_cm2))

lm_recruitment
ggsave("figures/lm_recruitment_redo.png", lm_recruitment, width = 2.75, height = 2.75, units = "in")

