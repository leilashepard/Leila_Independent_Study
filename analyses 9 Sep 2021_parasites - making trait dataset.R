library(tidyr)
library(dplyr)
library(ape)
library(nlme)
library(geiger)
library(ggplot2)
library(ggtree)
library(PerformanceAnalytics)

setwd("~/Git/Leila_Independent_Study")

#read in data from various sources
#tree from 10ktrees
tree = read.nexus("consensusTree_10kTrees_Primates_Version3_8.9.21.nex")

#absolute counts - wbc data from PAD :: start with totals 
leuks_all1 = read.csv("./PAD data/PAD_leukocytes_strepsirrhines_all.csv", header=T)
leuks_all = leuks_all1[-1]
colnames(leuks_all)[1] = "Common_Name"
leuks_all$Common_Name = tolower(leuks_all$Common_Name)

chart.Correlation(leuks_all[2:9])

#breeding seasonality data from DLC database
breed = read.csv("Zehr et al data edit names 8.5.21.csv", header=T)
str(breed)
breed$Species = gsub(pattern=" ", replacement="_", breed$Latin_Name)
breed$Common_Name = tolower(breed$Common_Name)
breed$Common_Name[breed$R_Pattern_Breeding=="NS"]


leuk_breed = leuks_all %>% left_join(breed)
leuk_breed$Species[leuk_breed$Common_Name=="golden-crowned sifaka"]="Propithecus_diadema"
leuk_breed$R_Pattern_Breeding[leuk_breed$Common_Name=="golden-crowned sifaka"]="S"


#breeding seasonality data from Heldstab et al. paper
held_data = read.csv("Heldstab et al. data.csv", header=T)

held_data$Species=gsub(" ","_",held_data$Species)
held_data$Species=gsub("\\*","",held_data$Species)

leuk_breed$Species[!leuk_breed$Species %in% held_data$Species]

str(held_data)

leuk_breed_2 = leuk_breed %>% left_join(held_data[c(1,14,16)])

#mating system data
matingsys = read.csv("Kling & Wright Mating System 2019.csv", header=T)
leuk_breed_sys = leuk_breed_2 %>% left_join(matingsys)


# other season data and mating system data ... 
leila_data = read.csv("Leila mating season system classification.csv", header=T)
leuk_breed_sys_leila = leuk_breed_sys %>% left_join(leila_data, by=c("Species" = "Scientific_Name_DLC"))
leuk_breed_sys_leila$Species

plot(Seasonality.natural.habitat~MatingSeasDur, data=leuk_breed_sys_leila)
comp_mod = lm(Seasonality.natural.habitat~MatingSeasDur, data=leuk_breed_sys_leila)
summary(comp_mod) #almost significant (p=0.06) positive correlation between Leila's classification 
#of mating season and Heldstab/van Schaik's, for the 9 overlapping species in the two datasets
nobs(comp_mod)

leuk_breed_sys_leila$Seasonality.natural.habitat

# add luepold data on testes mass
luepold = read.csv("Luepold_etal_Data_subset.csv", header=T)
luepold$Species = gsub(" ", "_", luepold$Species)

sort(leuk_breed_sys_leila$Species[leuk_breed_sys_leila$Species %in% luepold$Species])

leuk_breed_sys_test = leuk_breed_sys_leila %>% left_join(luepold, by="Species")




#leuk_breed_sys_leila$Species[leuk_breed_sys_leila$Species %in% luepold$Species]

#sexual dimorphism data from kappeler
#dimorph = read.csv("Kappeler 1993_sexual dimorphism only.csv", header=T)
#dimorph$Species = gsub(pattern=" ", replacement="_", dimorph$Species)
#length(dimorph$Relative_testes_size[!is.na(dimorph$Relative_testes_size)])
#dimorph$Relative_testes_size[dimorph$Species %in% leuk_breed_sys_leila$Species]

#going to use data from Propithecus verreauxi for Propithecus coquereli because they
#were considered the same species when the Kappeler paper was published
#dimorph$Species[which(dimorph$Species == "Propithecus_verreauxi")] = "Propithecus_coquereli"
#wbcdata2 = read.csv("sexual dimorphism for all species and all wbc types.csv", header=T)

#leuk_breed_sys_leila_dim = leuk_breed_sys_leila %>% left_join(dimorph, by="Species")
#leuk_breed_sys_leila_dim$Species[!leuk_breed_sys_leila_dim$Species %in% dimorph$Species]


#leuk_breed_sys_leila$Species[leuk_breed_sys_leila$Species %in% luepold$Species]
#sort(leuk_breed_sys_leila$Species[!leuk_breed_sys_leila$Species %in% luepold$Species])


#read in parasite data
parasites = read.csv("PSR data for lemurs.csv", header=T)
colnames(parasites)[2]="Species"
leuk_breed_sys_test2 = merge(leuk_breed_sys_test,parasites, by="Species",all=T)
str(leuk_breed_sys_test2)

traits_data = leuk_breed_sys_test2
View(traits_data)

str(traits_data)

#write.csv(traits_data, file="trait_data_strepsirrhines.csv")





