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

#mating system data
matingsys = read.csv("Kling & Wright Mating System 2019.csv", header=T)
leuk_breed_sys = leuk_breed %>% left_join(matingsys)


# other season data and mating system data ... 
leila_data = read.csv("Leila mating season system classification.csv", header=T)
leuk_breed_sys_leila = leuk_breed_sys %>% left_join(leila_data, by=c("Species" = "Scientific_Name_DLC"))


#sexual dimorphism data from kappeler
dimorph = read.csv("Kappeler 1993_sexual dimorphism only.csv", header=T)
dimorph$Species = gsub(pattern=" ", replacement="_", dimorph$Species)

#going to use data from Propithecus verreauxi for Propithecus coquereli because they
#were considered the same species when the Kappeler paper was published
dimorph$Species[which(dimorph$Species == "Propithecus_verreauxi")] = "Propithecus_coquereli"
#wbcdata2 = read.csv("sexual dimorphism for all species and all wbc types.csv", header=T)

leuk_breed_sys_leila_dim = leuk_breed_sys_leila %>% left_join(dimorph, by="Species")
leuk_breed_sys_leila_dim$Species[!leuk_breed_sys_leila_dim$Species %in% dimorph$Species]

#read in parasite data
parasites = read.csv("PSR data for lemurs.csv", header=T)
colnames(parasites)[2]="Species"
leuk_breed_sys_leila2 = merge(leuk_breed_sys_leila,parasites, by="Species",all=T)
str(leuk_breed_sys_leila2)

traits_data = leuk_breed_sys_leila2[c(35,c(1:34,36:41,43:49))]
View(traits_data)

str(traits_data)

plot(mean_value_White.Blood.Cells~mean_value_Body.Weight, data=traits_data)
plot(mean_value_Red.Blood.Cells~mean_value_Body.Weight, data=traits_data)
plot(mean_value_Lymphocytes.Abs~mean_value_Body.Weight, data=traits_data)
plot(mean_value_Neutrophil.Seg.Abs~mean_value_Body.Weight, data=traits_data)
plot(mean_value_Basophils.Abs~mean_value_Body.Weight, data=traits_data)
plot(mean_value_Eosinophils.Abs~mean_value_Body.Weight, data=traits_data)
plot(mean_value_Monocytes.Abs~mean_value_Body.Weight, data=traits_data)


plot(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding), data=traits_data) # more seasonal than non-seasonal
plot(mean_value_White.Blood.Cells~factor(Males_per_female_CRA), data=traits_data) # many more "multiple" than "single"

plot(mean_value_Red.Blood.Cells~factor(R_Pattern_Breeding), data=traits_data) # more seasonal than non-seasonal
plot(mean_value_Red.Blood.Cells~factor(Males_per_female_CRA), data=traits_data) # many more "multiple" than "single"

plot(mean_value_Lymphocytes.Abs~factor(R_Pattern_Breeding), data=traits_data) # more seasonal than non-seasonal
plot(mean_value_Lymphocytes.Abs~factor(Males_per_female_CRA), data=traits_data) # many more "multiple" than "single"

plot(mean_value_Basophils.Abs~factor(R_Pattern_Breeding), data=traits_data) # more seasonal than non-seasonal
plot(mean_value_Basophils.Abs~factor(Males_per_female_CRA), data=traits_data) # many more "multiple" than "single"

plot(mean_value_Eosinophils.Abs~factor(R_Pattern_Breeding), data=traits_data) # more seasonal than non-seasonal
plot(mean_value_Eosinophils.Abs~factor(Males_per_female_CRA), data=traits_data) # many more "multiple" than "single"

plot(mean_value_Monocytes.Abs~factor(R_Pattern_Breeding), data=traits_data) # more seasonal than non-seasonal
plot(mean_value_Monocytes.Abs~factor(Males_per_female_CRA), data=traits_data) # many more "multiple" than "single"




# phylogenetic least squares regressions 
tree$tip.label[tree$tip.label=="Eulemur_fulvus_albifrons"] = "Eulemur_albifrons"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_collaris"] = "Eulemur_collaris"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_fulvus"] = "Eulemur_fulvus"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_rufus"] = "Eulemur_rufus"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_sanfordi"] = "Eulemur_sanfordi"
tree$tip.label[tree$tip.label=="Eulemur_macaco_flavifrons"] = "Eulemur_flavifrons"
tree$tip.label[tree$tip.label=="Eulemur_macaco_macaco"] = "Eulemur_macaco"
tree$tip.label[tree$tip.label=="Varecia_variegata_variegata"] = "Varecia_variegata"


nameck = name.check(tree, traits_data, data.names=traits_data$Species)
nameck
#tree_sub=drop.tip(tree, nameck$tree_not_data)
#tree=tree_sub


p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
p1 = facet_plot(p, panel="Part", data=traits_data, geom=geom_point,
                aes(x=0, color=factor(Males_per_female_CRA))) 
p2 = facet_plot(p1, panel = "Seas", data=traits_data, geom=geom_point,
                      aes(x=0, color=factor(R_Pattern_Breeding)))
p2 + theme_tree2(panel.spacing = unit(c(5,0), "lines")) + theme(legend.position="none") +
  scale_color_manual(values=c("gray","black","gray","black"))


#WBC ~ males per female
pgls1 = gls(scale(mean_value_White.Blood.Cells)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1) # p = 0.18
plot(scale(mean_value_White.Blood.Cells)~factor(Males_per_female_CRA), data=traits_data)

#Figure 1
ggplot(data=traits_data, aes(x=factor(Males_per_female_CRA, labels=c("Multiple\n(n=20 spp)","Single\n(n=4 spp)")), y=mean_value_White.Blood.Cells, fill=factor(Males_per_female_CRA))) +
  geom_boxplot() + theme_bw() + ylab("Mean White Blood Cells (10^3/mm^3)") + xlab("Sexual Partners") + 
  theme(legend.position="none") +  scale_fill_brewer(palette="Paired") 

pgls1 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) # * p=0.02

plot(scale(mean_value_Neutrophil.Seg.Abs)~factor(Males_per_female_CRA), data=traits_data)
p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
p1 = facet_plot(p, panel="Neutrophils", data=traits_data, geom=geom_point,
                aes(x=mean_value_Neutrophil.Seg.Abs, color=factor(Males_per_female_CRA))) 
p1 + theme_tree2(panel.spacing = unit(5, "lines")) + theme(legend.position='none')

pgls1 = gls(scale(mean_value_Lymphocytes.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls1)


pgls1 = gls(scale(mean_value_Basophils.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Eosinophils.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls1)# p = 0.08
p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
p1 = facet_plot(p, panel="Eosinophils", data=traits_data, geom=geom_point,
                aes(x=mean_value_Eosinophils.Abs, color=factor(Males_per_female_CRA))) 
p1 + theme_tree2(panel.spacing = unit(5, "lines")) + theme(legend.position='none')


pgls1 = gls(scale(mean_value_Monocytes.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls1)


pgls1 = gls(scale(mean_value_Red.Blood.Cells)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE),
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls1)# *  p = 0.0003


# WBC ~ relative testes size -- sample size gets a lot smaller
trait_data_dim = leuk_breed_sys_leila_dim
trait_data_dimsub = trait_data_dim[-which(is.na(trait_data_dim$Relative_testes_size)),]


nameck = name.check(tree, trait_data_dimsub, data.names=trait_data_dimsub$Species)
nameck
tree_sub=drop.tip(tree, nameck$tree_not_data)


pgls1 = gls(scale(mean_value_White.Blood.Cells)~scale(Relative_testes_size), data=trait_data_dimsub, correlation = corPagel(0.5, tree_sub, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1)

pgls1 = gls(scale(mean_value_Neutrophil.Seg.Abs)~scale(Relative_testes_size), data=trait_data_dimsub, correlation = corPagel(0.5, tree_sub, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) # p = 0.1

pgls1 = gls(scale(mean_value_Lymphocytes.Abs)~scale(Relative_testes_size), data=trait_data_dimsub, correlation = corPagel(0.5, tree_sub, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Basophils.Abs)~scale(Relative_testes_size), data=trait_data_dimsub, correlation = corPagel(0.5, tree_sub, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Monocytes.Abs)~scale(Relative_testes_size), data=trait_data_dimsub, correlation = corPagel(0.5, tree_sub, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Eosinophils.Abs)~scale(Relative_testes_size), data=trait_data_dimsub, correlation = corPagel(0.5, tree_sub, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Red.Blood.Cells)~scale(Relative_testes_size), data=trait_data_dimsub, correlation = corPagel(0.5, tree_sub, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls1) # p = 0.01

plot(scale(mean_value_White.Blood.Cells)~scale(Relative_testes_size), data=trait_data_dimsub)
plot(scale(mean_value_Lymphocytes.Abs)~scale(Relative_testes_size), data=trait_data_dimsub)
plot(scale(mean_value_Monocytes.Abs)~scale(Relative_testes_size), data=trait_data_dimsub)
plot(scale(mean_value_Basophils.Abs)~scale(Relative_testes_size), data=trait_data_dimsub)
plot(scale(mean_value_Eosinophils.Abs)~scale(Relative_testes_size), data=trait_data_dimsub)
plot(scale(mean_value_Red.Blood.Cells)~scale(Relative_testes_size), data=trait_data_dimsub)
plot(scale(mean_value_Neutrophil.Seg.Abs)~scale(Relative_testes_size), data=trait_data_dimsub)



# WBC ~ seasonal vs nonseasonal
#traits_data_sub2 = subset(traits_data, !is.na(R_Pattern_Breeding))
#nameck = name.check(tree, traits_data_sub2, data.names=traits_data_sub2$Species)
#nameck
#tree_sub2=drop.tip(tree, nameck$tree_not_data)

pgls2 = gls(scale(mean_value_White.Blood.Cells)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0002

plot(scale(mean_value_White.Blood.Cells)~factor(R_Pattern_Breeding), data=traits_data)
p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
p1 = facet_plot(p, panel="White blood cells", data=traits_data, geom=geom_point,
                aes(x=mean_value_White.Blood.Cells, color=factor(R_Pattern_Breeding))) +
  #theme(legend.position="none") + 
  scale_color_manual(values=c("#ff7f00", "#33a02c"), name="Seasonal in\n captivity?",labels=c("No","Yes"))
p1 + theme_tree2(panel.spacing = unit(7, "lines"))



pgls2 = gls(scale(mean_value_Lymphocytes.Abs)~factor(R_Pattern_Breeding), data=traits_data_sub2, 
            correlation = corPagel(0.7, tree_sub2, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0004


pgls2 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(R_Pattern_Breeding), data=traits_data_sub2, 
            correlation = corPagel(0.7, tree_sub2, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls2)

pgls2 = gls(scale(mean_value_Eosinophils.Abs)~factor(R_Pattern_Breeding), data=traits_data_sub2, 
            correlation = corPagel(0.7, tree_sub2, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls2)

pgls2 = gls(scale(mean_value_Monocytes.Abs)~factor(R_Pattern_Breeding), data=traits_data_sub2, 
            correlation = corPagel(0.7, tree_sub2, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls2)

pgls2 = gls(scale(mean_value_Basophils.Abs)~factor(R_Pattern_Breeding), data=traits_data_sub2, 
            correlation = corPagel(0.7, tree_sub2, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs))
summary(pgls2) #lambda 0


pgls2 = gls(scale(mean_value_Red.Blood.Cells)~factor(R_Pattern_Breeding), data=traits_data_sub2, 
            correlation = corPagel(0.7, tree_sub2, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls2)


# WBC ~ duration of mating season - in the WILD 
traits_data_sub3 = subset(traits_data, !is.na(MatingSeasDur))
nameck = name.check(tree, traits_data_sub3, data.names=traits_data_sub3$Species)
nameck
tree_sub3=drop.tip(tree, nameck$tree_not_data)

traits_data_sub3$MatingSeasDur = as.numeric(traits_data_sub3$MatingSeasDur)
hist(traits_data_sub3$MatingSeasDur)



#
pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur, data=traits_data_sub3, 
            correlation = corPagel(0.7, tree_sub3, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # *p=0.03
plot(mean_value_White.Blood.Cells~MatingSeasDur, data=traits_data_sub3)

traits_data_sub3$preds = predict(pgls3)
ggplot(data=traits_data_sub3, aes(x=MatingSeasDur, y=mean_value_White.Blood.Cells)) + 
  geom_smooth(method="lm", formula = y~x, aes(y=preds), color="#b2df8a") +
  geom_point(color="#33a02c", size=2) +
  xlab("Mating season duration in wild (months)") +
  ylab("Mean white blood cell count (10^3/mm^3)") + theme_bw() 



p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
p1 = facet_plot(p, panel="White blood cells", data=traits_data, geom=geom_point,
                aes(x=scale(mean_value_White.Blood.Cells)), pch=1) + geom_point(data=traits_data, aes(x=scale(MatingSeasDur)), pch=2)
  facet_plot(p, panel="Mating season - wild", data=traits_data, geom=geom_point,
             aes(x=scale(MatingSeasDur)), pch=2) +
  theme(legend.position="none") 
p1

p1 = facet_plot(p, panel='WBC-blue, MatSeas-green (z-transf)', data=traits_data, geom=geom_point, 
                     aes(x=scale(mean_value_White.Blood.Cells)),pch=1, color='blue3') 
p2 = facet_plot(p1, panel='WBC-blue, MatSeas-green (z-transf)', data=traits_data, geom=geom_point, 
                aes(x=scale(MatingSeasDur)),pch=2, color='green') 
p2 + theme_tree2(panel.spacing = unit(5.5, "lines"))    

myplot2 = facet_plot(myplot1, panel="Trait", data=data, geom=geom_point, 
                     aes(x=Body_Mass*20 ), col="orange")

myplot2 + theme_tree2(panel.spacing = unit(5, "lines"))


pgls3 = gls(scale(mean_value_Neutrophil.Seg.Abs)~scale(MatingSeasDur), data=traits_data_sub3, 
            correlation = corPagel(0.7, tree_sub3, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls3) #p=0.03, lambda 0


pgls3 = gls(scale(mean_value_Lymphocytes.Abs)~scale(MatingSeasDur), data=traits_data_sub3, 
            correlation = corPagel(0.7, tree_sub3, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls3)


pgls3 = gls(scale(mean_value_Eosinophils.Abs)~scale(MatingSeasDur), data=traits_data_sub3, 
            correlation = corPagel(0.7, tree_sub3, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls3) 

pgls3 = gls(scale(mean_value_Basophils.Abs)~scale(MatingSeasDur), data=traits_data_sub3, 
            correlation = corPagel(0.7, tree_sub3, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs))
summary(pgls3) 

pgls3 = gls(scale(mean_value_Monocytes.Abs)~scale(MatingSeasDur), data=traits_data_sub3, 
            correlation = corPagel(0.7, tree_sub3, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls3) # * p = 0.03


pgls3 = gls(scale(mean_value_Red.Blood.Cells)~scale(MatingSeasDur), data=traits_data_sub3, 
            correlation = corPagel(0.7, tree_sub3, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Red.Blood.Cells)) 
summary(pgls3) 



#body mass

#parasites ~ mating season duration
plot(prop_close~MatingSeasDur, data=traits_data)
plot(prop_direct~MatingSeasDur, data=traits_data)

plot(PSR_close~MatingSeasDur, data=traits_data)
plot(PSR_direct~MatingSeasDur, data=traits_data)


#parasites ~ wbcs
plot(prop_close~mean_value_White.Blood.Cells, data=traits_data)
plot(prop_direct~mean_value_White.Blood.Cells, data=traits_data)

plot(PSR_close~mean_value_White.Blood.Cells, data=traits_data)
plot(PSR_direct~mean_value_White.Blood.Cells, data=traits_data)


### proportion directly transmitted parasites
traits_data_sub4 = traits_data[! (is.na(traits_data$prop_direct) | is.na(traits_data$mean_value_White.Blood.Cells)),]
nameck = name.check(tree, traits_data_sub4, data.names=traits_data_sub4$Species)
nameck
tree_sub4=drop.tip(tree, nameck$tree_not_data)

pgls4 = gls(scale(mean_value_White.Blood.Cells)~scale(prop_direct), data=traits_data_sub4, 
            correlation = corPagel(0.1, tree_sub4, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls4) 

pgls4 = gls(scale(mean_value_Basophils.Abs)~scale(prop_direct), data=traits_data_sub4, 
            correlation = corPagel(0.1, tree_sub4, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs)) 
summary(pgls4) 

pgls4 = gls(scale(mean_value_Eosinophils.Abs)~scale(prop_direct), data=traits_data_sub4, 
            correlation = corPagel(0.1, tree_sub4, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs)) 
summary(pgls4) #p=0.062

pgls4 = gls(scale(mean_value_Lymphocytes.Abs)~scale(prop_direct), data=traits_data_sub4, 
            correlation = corPagel(0.1, tree_sub4, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs)) 
summary(pgls4) #p=0.03

pgls4 = gls(scale(mean_value_Monocytes.Abs)~scale(prop_direct), data=traits_data_sub4, 
            correlation = corPagel(0.1, tree_sub4, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls4) 

pgls4 = gls(scale(mean_value_Neutrophil.Seg.Abs)~scale(prop_direct), data=traits_data_sub4, 
            correlation = corPagel(0.1, tree_sub4, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs)) 
summary(pgls4) 



#parasites ~ wbcs
plot(prop_close~mean_value_White.Blood.Cells, data=traits_data)
plot(prop_direct~mean_value_White.Blood.Cells, data=traits_data)

plot(PSR_close~mean_value_White.Blood.Cells, data=traits_data)
plot(PSR_direct~mean_value_White.Blood.Cells, data=traits_data)



### proportion close transmitted
traits_data_sub5 = traits_data[! (is.na(traits_data$prop_close) | is.na(traits_data$mean_value_White.Blood.Cells)),]
nameck = name.check(tree, traits_data_sub5, data.names=traits_data_sub5$Species)
nameck
tree_sub5=drop.tip(tree, nameck$tree_not_data)

pgls5 = gls(scale(mean_value_White.Blood.Cells)~scale(prop_close), data=traits_data_sub5, 
            correlation = corPagel(0.7, tree_sub5, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls5) 

pgls5 = gls(scale(mean_value_Basophils.Abs)~scale(prop_close), data=traits_data_sub5, 
            correlation = corPagel(0.7, tree_sub5, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs)) 
summary(pgls5) 

pgls5 = gls(scale(mean_value_Eosinophils.Abs)~scale(prop_close), data=traits_data_sub5, 
            correlation = corPagel(0.7, tree_sub5, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs)) 
summary(pgls5) #p=0.03

pgls5 = gls(scale(mean_value_Lymphocytes.Abs)~scale(prop_close), data=traits_data_sub5, 
            correlation = corPagel(0.7, tree_sub5, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs)) 
summary(pgls5)

pgls5 = gls(scale(mean_value_Monocytes.Abs)~scale(prop_close), data=traits_data_sub5, 
            correlation = corPagel(0.7, tree_sub5, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls5) #0.059

pgls5 = gls(scale(mean_value_Neutrophil.Seg.Abs)~scale(prop_close), data=traits_data_sub5, 
            correlation = corPagel(0.7, tree_sub5, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs)) 
summary(pgls5) 







traits_data_sub6 = traits_data[! (is.na(traits_data$prop_close) | is.na(traits_data$mean_value_White.Blood.Cells) | is.na(traits_data$MatingSeasDur)),]
nameck = name.check(tree, traits_data_sub6, data.names=traits_data_sub6$Species)
nameck
tree_sub6=drop.tip(tree, nameck$tree_not_data)

pgls6 = gls(scale(mean_value_White.Blood.Cells)~scale(prop_close) + scale(MatingSeasDur), data=traits_data_sub6, 
            correlation = corPagel(0.7, tree_sub6, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls6) 






