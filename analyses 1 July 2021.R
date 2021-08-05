library(tidyr)
library(dplyr)
library(ape)
library(nlme)

#read in data from various sources
#tree from 10ktrees
tree = read.nexus("consensusTree_10kTrees_Primates_lumped.nex")

#trait data - mating season dur from the literature, wbc and body mass data from primate aging database
traits = read.csv("cumulative_data_finalll.csv", header=T)

#sexual dimorphism data from kappeler
dimorph = read.csv("Kappeler 1993 dimorphism data.csv", header=T)

str(dimorph)
dimorph$Species = gsub(pattern=" ", replacement="_", dimorph$Species)

dimorph$Species[! dimorph$Species %in% traits$Scientific_Name]
traits$Scientific_Name[! traits$Scientific_Name %in% dimorph$Species]

str(traits$Scientific_Name)
traits_lump = traits %>% separate(col= Scientific_Name, into=c("genus","species","subspecies"), "_") %>%
  unite(col = "genus_species", c("genus", "species")) %>% group_by(genus_species) %>%
  summarize(MatingSeasDur = mean(MatingSeasDur), Mating_System_Leila=unique(Mating_System), AllWBC = mean(AllWBC),
            BodyMass = mean(BodyWeight), Lymphocytes=mean(Lymphocytes), Monocytes = mean(Monocytes), 
            Basophils = mean(Basophils), Eosinophils = mean(Eosinophils), Neutrophils = mean(Neutrophils), 
            wbc_n = sum(wbc_n_ind), bm_n = sum(body_weight_n_ind), lymph_n=sum(lymphocytes_n_ind),
            mono_n = sum(monocytes_n_ind), baso_n = sum(basophils_n_ind), eosin_n = sum(eosinophils_n_ind), 
            neut_n = sum(neutrophils_seg_n_ind))
traits_lump

#if we lump the subspecies into species, most species are represented in the Kappeler data
traits_lump$genus_species[traits_lump$genus_species %in% dimorph$Species]
traits_lump$genus_species[!traits_lump$genus_species %in% dimorph$Species]

#traits (subspecies lumped) with kappeler dimorphism data merged 
colnames(traits_lump)[1]="Species"
traits_lp2 = merge(traits_lump, dimorph, by="Species", all.x=T, all.y=F)

#kappeler's assignment of mating system
boxplot(Relative_testes_size~Mating_system, data=traits_lp2)

#leila's assignment of mating system
boxplot(Relative_testes_size~Mating_System_Leila, data=traits_lp2)

linmod1 = lm(Relative_testes_size~Mating_system, data=traits_lp2)
summary(linmod1)

linmod2 = lm(Relative_testes_size~Mating_System_Leila, data=traits_lp2)
summary(linmod2)

tree$tip.label[15] = "Varecia_variegata"
tree$tip.label[5] = "Eulemur_macaco"


traits_lp3 = subset(traits_lp2, !is.na(Relative_testes_size))
tree.cplt = drop.tip(tree,setdiff(tree$tip.label,traits_lp3$Species))
tree.cplt$tip.label
traits_lp3$Species

olinmod1= lm(Relative_testes_size~Mating_system, data=traits_lp3)


pglinmod1 = gls(Lymphocytes~Relative_testes_size, 
                data=traits_lp3, correlation = corPagel(value=0.5, phy=tree.cplt, form=~Species, fixed=FALSE))
summary(pglinmod1)




