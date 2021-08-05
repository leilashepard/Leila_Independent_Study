library(tidyr)
library(dplyr)
library(ape)
library(nlme)

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
wbcdata_lump = wbcdata %>% separate(col= Scientific_Name, into=c("genus","species","subspecies"), "_") %>%
  unite(col = "genus_species", c("genus", "species")) %>% group_by(genus_species) %>%
  summarize(MatingSeasDur = mean(MatingSeasDur), Mating_System_Leila=unique(Mating_System), AllWBC = mean(AllWBC),
            BodyMass = mean(BodyWeight), Lymphocytes=mean(Lymphocytes), Monocytes = mean(Monocytes), 
            Basophils = mean(Basophils), Eosinophils = mean(Eosinophils), Neutrophils = mean(Neutrophils), 
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

plot(AllWBC~CRA_Peak_Breeding_Season_Duration, data=traits_data)
summary(lm(AllWBC~CRA_Peak_Breeding_Season_Duration + BodyMass, data=traits_data))

plot(AllWBC~Relative_testes_size, data=traits_data)
summary(lm(AllWBC~Relative_testes_size, data=traits_data))


# have to get actual counts (% x Total wbcs) of the subtypes

plot(Relative_testes_size~factor(Mating_system), data=traits_data)

points(Relative_testes_size~factor(Mating_system), data=traits_data)




# phylogenetics



