library(tidyr)
library(dplyr)
library(ape)
library(nlme)
library(geiger)

#downloaded the WBC data from PAD on 5 Aug 2021


#read in data from various sources
#tree from 10ktrees
tree = read.nexus("consensusTree_10kTrees_Primates_lumped.nex")


#sexual dimorphism data from kappeler
dimorph = read.csv("Kappeler 1993 dimorphism data.csv", header=T)
dimorph$Species = gsub(pattern=" ", replacement="_", dimorph$Species)

#going to use data from Propithecus verreauxi for Propithecus coquereli because they
#were considered the same species when the Kappeler paper was published
dimorph$Species[which(dimorph$Species == "Propithecus_verreauxi")] = "Propithecus_coquereli"


#breeding data from DLC database
breed = read.csv("Zehr et al data edit names 8.5.21.csv", header=T)
str(breed)
breed$Species = gsub(pattern=" ", replacement="_", breed$Latin_Name)

dimorph$Species[which(dimorph$Species %in% breed$Species)]
breed$Species[which(breed$Species %in% dimorph$Species)]

dimorph$Species[!which(dimorph$Species %in% breed$Species)]
breed$Species[!which(breed$Species %in% dimorph$Species)]


# WBC data from Primate Aging Database
wbcdata = read.csv("cumulative_data_finalll.csv", header=T)
#wbcdata2 = read.csv("sexual dimorphism for all species and all wbc types.csv", header=T)

str(wbcdata$Scientific_Name)

#WBC data is corrected to counts here -- rather than percentages
wbcdata_lump = wbcdata %>% separate(col= Scientific_Name, into=c("genus","species","subspecies"), "_") %>%
  unite(col = "genus_species", c("genus", "species")) %>% group_by(genus_species) %>%
  summarize(MatingSeasDur = mean(MatingSeasDur), Mating_System_Leila=unique(Mating_System), AllWBC = mean(AllWBC),
            BodyMass = mean(BodyWeight), Lymphocytes=mean(Lymphocytes)/(AllWBC*100), Monocytes = mean(Monocytes)/(AllWBC*100), 
            Basophils = mean(Basophils)/(AllWBC*100), Eosinophils = mean(Eosinophils)/(AllWBC*100), Neutrophils = mean(Neutrophils)/(AllWBC*100), 
            wbc_n = sum(wbc_n_ind), bm_n = sum(body_weight_n_ind), lymph_n=sum(lymphocytes_n_ind),
            mono_n = sum(monocytes_n_ind), baso_n = sum(basophils_n_ind), eosin_n = sum(eosinophils_n_ind), 
            neut_n = sum(neutrophils_seg_n_ind))
wbcdata_lump
colnames(wbcdata_lump)[1] = "Species"

trait_data1 = merge(breed, wbcdata_lump, by="Species")
traits_data = merge(trait_data1, dimorph, by="Species", all.x=T)

#Assigning "G" as the mating system of Varecia rubra based on report from primary literature
#Vasey, N. 2007. "The breeding system of wild red ruffed lemurs (Varecia rubra): a preliminary report"

traits_data$Mating_system[traits_data$Species=="Varecia_rubra"] = "G"

plot(AllWBC~BodyMass, data=traits_data)

plot(AllWBC~factor(R_Pattern_Breeding), data=traits_data)
plot(AllWBC~factor(Mating_system), data=traits_data)
points(AllWBC~factor(Mating_system), data=traits_data)

plot(AllWBC~CRA_Peak_Breeding_Season_Duration, data=traits_data)
summary(lm(AllWBC~CRA_Peak_Breeding_Season_Duration + BodyMass, data=traits_data))

plot(AllWBC~Relative_testes_size, data=traits_data)
summary(lm(AllWBC~Relative_testes_size, data=traits_data))


# have to get absolute counts of the subtypes

plot(Relative_testes_size~factor(Mating_system), data=traits_data)

points(Relative_testes_size~factor(Mating_system), data=traits_data)




# phylogenetic least squares regressions 

tree$tip.label[tree$tip.label=="Eulemur_macaco_macaco"] = "Eulemur_macaco"
tree$tip.label[tree$tip.label=="Varecia_variegata_variegata"] = "Varecia_variegata"
nameck = name.check(tree, traits_data, data.names=traits_data$Species)
nameck
tree_sub=drop.tip(tree, nameck$tree_not_data)
tree=tree_sub

#WBC~body mass
pgls1 = gls(scale(AllWBC)~scale(BodyMass), data=traits_data, correlation = corPagel(0.7, tree, form=~Species, fixed=FALSE))
summary(pgls1)

traits_data$Mating_system_binary=as.character(traits_data$Mating_system)
traits_data$Mating_system_binary[which(traits_data$Mating_system_binary=="S")] = "G"
traits_data$Mating_system_binary[which(traits_data$Mating_system_binary=="G")] = 2 # commonly more than one mate
traits_data$Mating_system_binary[which(traits_data$Mating_system_binary=="P")] = 1 # only one mate usually
traits_data$Mating_system_binary = factor(traits_data$Mating_system_binary)

# WBC ~ mating system
pgls2 = gls(scale(AllWBC)~factor(Mating_system), data=traits_data, correlation = corPagel(0.7, tree, form=~Species, fixed=FALSE))
summary(pgls2)

# WBC ~ relative testes mass - incomplete data
traits_data_sub1 = subset(traits_data, !is.na(Relative_testes_size))
pgls3 = gls(scale(AllWBC)~scale(Relative_testes_size), data=traits_data_sub1, correlation = corPagel(0.7, tree, form=~Species, fixed=FALSE))
summary(pgls3)

# WBC ~ mating season yes/no
pgls4 = gls(scale(AllWBC)~factor(R_Pattern_Breeding), data=traits_data, correlation = corPagel(0.7, tree, form=~Species, fixed=FALSE), method="REML")
summary(pgls4)

# WBC ~ mating season duration in seasonal species only
traits_data_sub2 = subset(traits_data, !is.na(CRA_Peak_Breeding_Season_Duration))
pgls4 = gls(scale(AllWBC)~CRA_Peak_Breeding_Season_Duration, data=traits_data_sub2, correlation = corPagel(0.7, tree, form=~Species, fixed=FALSE), method="REML")
summary(pgls4)

# intra-specific patterns ... 


# Next steps: 
# process the absolute wbc measurement data - units are x10^3 per mm3
#(Except red blood cells are x10^6 per mm3)








