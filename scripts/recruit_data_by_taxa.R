#recruit data 
suz08 <- read_csv("data/Suzuki_08_Recruit_Counts.csv")
suz12 <- read_csv("data/Suzuki_12_Recruit_Counts.csv")
transmit <- read_csv("data/Swain_2018_Zoox_In_Larvae.csv")
zoox_sub <- transmit %>% select(`Coral species`, `Symbiodinium in propagules`, `Mode of larval development`, `Sexual system`)
colnames(zoox_sub) <- c("Taxon", "ZooxIn", "ReproMode", "Sex")
zoox_sub <- zoox_sub %>% separate(Taxon, into = c("Genus", "Species")) 

gr <- read_csv("data/ctdb_growth_rate/data_20230425.csv")
gr <- gr %>% filter(trait_name == "Growth rate") %>% filter(standard_unit == "mm yr^-1")
gr_sum <- gr %>% mutate(value = as.numeric(value)) %>% group_by(specie_name) %>% summarise(mean_gr = mean(value)) %>% separate(specie_name, into = c("Genus", "Species"))


bleach <- read_csv("data/Swain_2016_bleaching_ranks_by_taxa.csv")
bleach <- bleach %>% separate(Taxon, into=c("Genus", "Species"))# %>% filter(Genus %in% HT)
bleach_sub <- bleach %>% select(Genus, Species, `Taxon-BRI (Bleaching Response Index) (%)                                                   
`)
colnames(bleach_sub) <- c("Genus", "Species", "TaxonBRI")



suz_across_habitat <- suz12 %>% group_by(Genus, Species) %>% mutate(RecruitCount = sum(RecruitCount)) %>% unique()

suz_bleach <- left_join(suz_across_habitat, bleach_sub)
summary(lm(suz_bleach$TaxonBRI~suz_bleach$RecruitCount))

library(ggrepel)
suz_bleach <- left_join(suz_bleach, gr_sum)
pacific <- suz_bleach %>% ggplot(aes(x=RecruitCount, y =TaxonBRI, col = mean_gr)) +  geom_smooth(method = "lm")+geom_point(size = 4, alpha = .8) + theme_classic() +
  ylab("Taxon-Specific Bleaching Index") + xlab("Total Recruits") + geom_label_repel(aes(label=paste(Species))) +
  theme(text=element_text(size = 10)) + scale_color_viridis_c()

pacific
ggsave("figures/acroporid_suziki_recruit.png", pacific, width = 3, height = 3, dpi = 400, units = "in" )
