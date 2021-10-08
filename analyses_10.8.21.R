library(tidyr)
library(dplyr)
library(ape)
library(nlme)
library(geiger)
library(ggplot2)
library(ggtree)
library(PerformanceAnalytics)

# read in trait data, constructed in script: analyses 9 Sep 2021_parasites - making trait dataset.R
traits_data = read.csv("trait_data_strepsirrhines.csv", header=T)

# read in tree from 10ktrees
tree = read.nexus("consensusTree_10kTrees_Primates_Version3_8.9.21.nex")

# names in tree are subspecies format, rename to match dataset
tree$tip.label[tree$tip.label=="Eulemur_fulvus_albifrons"] = "Eulemur_albifrons"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_collaris"] = "Eulemur_collaris"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_fulvus"] = "Eulemur_fulvus"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_rufus"] = "Eulemur_rufus"
tree$tip.label[tree$tip.label=="Eulemur_fulvus_sanfordi"] = "Eulemur_sanfordi"
tree$tip.label[tree$tip.label=="Eulemur_macaco_flavifrons"] = "Eulemur_flavifrons"
tree$tip.label[tree$tip.label=="Eulemur_macaco_macaco"] = "Eulemur_macaco"
tree$tip.label[tree$tip.label=="Varecia_variegata_variegata"] = "Varecia_variegata"

#check that names match
nameck = name.check(tree, traits_data, data.names=traits_data$Species)
nameck

#if names don't match...
#tree_sub=drop.tip(tree, nameck$tree_not_data)
#tree=tree_sub

#WBC ~ body mass
pgls_bm = gls(scale(mean_value_White.Blood.Cells)~scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls_bm) # p = 0.72; lambda = 0.55

#neutrophils ~ BM
pgls_bm = gls(scale(mean_value_Neutrophil.Seg.Abs)~scale(mean_value_Body.Weight), data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls_bm) # p = 0.009; lambda = -0.196... but lambda is negative... so I re-fit the model fixing lambda at 0.

pgls_bm = gls(scale(mean_value_Neutrophil.Seg.Abs)~scale(mean_value_Body.Weight), data=traits_data, 
              correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
              method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls_bm) #p= 0.008

#lymphocytes ~ BM
pgls_bm = gls(scale(mean_value_Lymphocytes.Abs)~scale(mean_value_Body.Weight), data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls_bm) #p= 0.08

#basophils ~ BM
pgls_bm = gls(scale(mean_value_Basophils.Abs)~scale(mean_value_Body.Weight), data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls_bm) #p= 0.70, but lambda is negative, so re-fit?

pgls_bm = gls(scale(mean_value_Basophils.Abs)~scale(mean_value_Body.Weight), data=traits_data, 
              correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
              method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls_bm) #p= 0.35

#eosinophils ~ BM
pgls_bm = gls(scale(mean_value_Eosinophils.Abs)~scale(mean_value_Body.Weight), data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls_bm) #p= 0.38

#monocytes ~ BM
pgls_bm = gls(scale(mean_value_Monocytes.Abs)~scale(mean_value_Body.Weight), data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls_bm) #p= 0.16

#RBC ~ BM
pgls_bm = gls(scale(mean_value_Red.Blood.Cells)~factor(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE),
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls_bm)# *  p = 0.16


######################
#WBC ~ males per female
pgls1 = gls(scale(mean_value_White.Blood.Cells)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1) # p = 0.20; lambda = 0.55

#Figure 1A
ggplot(data=traits_data, aes(x=factor(Males_per_female_CRA, labels=c("Multiple\n(n=20 spp)","Single\n(n=4 spp)")), y=mean_value_White.Blood.Cells, fill=factor(Males_per_female_CRA))) +
  geom_boxplot() + theme_bw() + ylab("Mean White Blood Cells (10^3/mm^3)") + xlab("Sexual Partners per Female") + 
  theme(legend.position="none") +  scale_fill_brewer(palette="Paired") 

#Neutrophils ~ males per female
pgls1 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) # * p=0.002, but lambda is negative... so I re-fit the model fixing lambda at 0.

pgls1 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) # * p=0.02 when lambda is fixed at 0; should report this result I think


#Figure 1B
ggplot(data=traits_data, aes(x=factor(Males_per_female_CRA, labels=c("Multiple\n(n=20 spp)","Single\n(n=4 spp)")), 
                             y=mean_value_Neutrophil.Seg.Abs/1000, fill=factor(Males_per_female_CRA))) +
  geom_boxplot() + theme_bw() + ylab("Mean Neutrophils (10^3/per mm^3)") + xlab("Sexual Partners per Female") + 
  theme(legend.position="none") +  scale_fill_brewer(palette="Paired") 

#p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
#p1 = facet_plot(p, panel="Neutrophils", data=traits_data, geom=geom_point,
#                aes(x=scale(mean_value_Neutrophil.Seg.Abs), color=factor(Males_per_female_CRA))) 
#p1 + theme_tree2(panel.spacing = unit(5, "lines")) + theme(legend.position='none')

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
summary(pgls1)


pgls1 = gls(scale(mean_value_Monocytes.Abs)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls1)


pgls1 = gls(scale(mean_value_Red.Blood.Cells)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE),
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls1)# *  p = 0.0002



# WBC ~ males per female, include BM covariate

######################
#WBC ~ males per female
pgls1 = gls(scale(mean_value_White.Blood.Cells)~factor(Males_per_female_CRA) + scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1) # p = 0.20; lambda = 0.55

#Neutrophils ~ males per female
pgls1 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(Males_per_female_CRA) + scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs)) # false convergence

pgls1 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(Males_per_female_CRA) + scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) # p=0.12 * neutrophil no longer significant when control for body mass (pvalue for BM = 0.0490)


#Figure 1B
ggplot(data=traits_data, aes(x=factor(Males_per_female_CRA, labels=c("Multiple\n(n=20 spp)","Single\n(n=4 spp)")), 
                             y=mean_value_Neutrophil.Seg.Abs/1000, fill=factor(Males_per_female_CRA))) +
  geom_boxplot() + theme_bw() + ylab("Mean Neutrophils (10^3/per mm^3)") + xlab("Sexual Partners per Female") + 
  theme(legend.position="none") +  scale_fill_brewer(palette="Paired") 

#p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
#p1 = facet_plot(p, panel="Neutrophils", data=traits_data, geom=geom_point,
#                aes(x=scale(mean_value_Neutrophil.Seg.Abs), color=factor(Males_per_female_CRA))) 
#p1 + theme_tree2(panel.spacing = unit(5, "lines")) + theme(legend.position='none')

pgls1 = gls(scale(mean_value_Lymphocytes.Abs)~factor(Males_per_female_CRA) + scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls1) #BM significant p = 0.03


pgls1 = gls(scale(mean_value_Basophils.Abs)~factor(Males_per_female_CRA) + scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls1) #lambda negative

pgls1 = gls(scale(mean_value_Basophils.Abs)~factor(Males_per_female_CRA) + scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls1) 


pgls1 = gls(scale(mean_value_Eosinophils.Abs)~factor(Males_per_female_CRA)+scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls1)


pgls1 = gls(scale(mean_value_Monocytes.Abs)~factor(Males_per_female_CRA)+scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls1)


pgls1 = gls(scale(mean_value_Red.Blood.Cells)~factor(Males_per_female_CRA)+scale(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE),
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls1)# *  p = 0.0007






# WBC ~ relative testes size -- sample size gets a lot smaller

traits_data_testes = subset(traits_data, !is.na(CombinedTestesMass.in.g))
plot(CombinedTestesMass.in.g~MaleBodyMass.TestesDataset..in.g, data=traits_data_testes)


nameck = name.check(tree, traits_data_testes, data.names=traits_data_testes$Species)
nameck
tree_testes=drop.tip(tree, nameck$tree_not_data)

pgls_testes = gls(CombinedTestesMass.in.g~MaleBodyMass.TestesDataset..in.g, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE),
            method="REML")
summary(pgls_testes)

pgls_testes = gls(CombinedTestesMass.in.g~MaleBodyMass.TestesDataset..in.g, data=traits_data_testes, 
                  correlation = corPagel(0, tree_testes, form=~Species, fixed=TRUE),
                  method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls_testes)
traits_data_testes$relative_testes_size_lp = residuals(pgls_testes) #rel testes size from luepold data




pgls1 = gls(mean_value_White.Blood.Cells~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1)

traits_data_testes$preds = predict(pgls1)

# Figure 2 #
ggplot(data=traits_data_testes, aes(x=relative_testes_size_lp, y=mean_value_White.Blood.Cells)) + 
  #geom_smooth(method="lm", formula = y~x, aes(y=preds), color="#8c96c6") +
  geom_point(color="#8c96c6", size=2) +
  xlab("Relative testes mass (residuals)") +
  ylab("Mean white blood cell count (10^3/mm^3)") + theme_bw() 

pgls1 = gls(scale(mean_value_Neutrophil.Seg.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) # p = 0.1

pgls1 = gls(scale(mean_value_Lymphocytes.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Lymphocytes.Abs)) #lambda > 1 !
summary(pgls1)

pgls1 = gls(scale(mean_value_Lymphocytes.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(1, tree_testes, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls1)


pgls1 = gls(scale(mean_value_Basophils.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Monocytes.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls1) # lambda negative

pgls1 = gls(scale(mean_value_Monocytes.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0, tree_testes, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls1) # lambda negative


pgls1 = gls(scale(mean_value_Eosinophils.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Eosinophils.Abs)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls1)

pgls1 = gls(scale(mean_value_Red.Blood.Cells)~scale(relative_testes_size_lp), data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls1) # p = 0.01

plot(scale(mean_value_White.Blood.Cells)~scale(relative_testes_size_lp), data=traits_data_testes)
plot(scale(mean_value_Lymphocytes.Abs)~scale(relative_testes_size_lp), data=traits_data_testes)
plot(scale(mean_value_Monocytes.Abs)~scale(relative_testes_size_lp), data=traits_data_testes)
plot(scale(mean_value_Basophils.Abs)~scale(relative_testes_size_lp), data=traits_data_testes)
plot(scale(mean_value_Eosinophils.Abs)~scale(relative_testes_size_lp), data=traits_data_testes)
plot(scale(mean_value_Red.Blood.Cells)~scale(relative_testes_size_lp), data=traits_data_testes)
plot(scale(mean_value_Neutrophil.Seg.Abs)~scale(relative_testes_size_lp), data=traits_data_testes)


################3
# WBC ~ seasonal vs nonseasonal (in captivity)
#traits_data_sub2 = subset(traits_data, !is.na(R_Pattern_Breeding))
#nameck = name.check(tree, traits_data_sub2, data.names=traits_data_sub2$Species)
#nameck
#tree_sub2=drop.tip(tree, nameck$tree_not_data)

pgls2 = gls(scale(mean_value_White.Blood.Cells)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0002, but lambda is negative

pgls2 = gls(scale(mean_value_White.Blood.Cells)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0004 with fixed lambda = 0

#same pattern controlling for body mass too 
pgls2 = gls(scale(mean_value_White.Blood.Cells)~factor(R_Pattern_Breeding) + scale(mean_value_Body.Weight) , data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0002, but lambda is negative

pgls2 = gls(scale(mean_value_White.Blood.Cells)~factor(R_Pattern_Breeding) + scale(mean_value_Body.Weight) , data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0004


#Figure 3
p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
p1 = facet_plot(p, panel="White blood cells", data=traits_data, geom=geom_point,
                aes(x=mean_value_White.Blood.Cells, color=factor(R_Pattern_Breeding))) +
  #theme(legend.position="none") + 
  scale_color_manual(values=c("#ff7f00", "#33a02c"), name="Seasonal in\n captivity?",labels=c("No","Yes"))
p1 + theme_tree2(panel.spacing = unit(7, "lines"))


pgls2 = gls(scale(mean_value_Lymphocytes.Abs)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.7, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0004


pgls2 = gls(scale(mean_value_Lymphocytes.Abs)~factor(R_Pattern_Breeding) + scale(mean_value_Body.Weight) , data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0004 - same pattern controlling for body mass


pgls2 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls2)

pgls2 = gls(scale(mean_value_Eosinophils.Abs)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls2)

pgls2 = gls(scale(mean_value_Monocytes.Abs)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls2)

pgls2 = gls(scale(mean_value_Basophils.Abs)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs))
summary(pgls2) #lambda 0


pgls2 = gls(scale(mean_value_Red.Blood.Cells)~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls2)


# WBC ~ duration of mating season - in the WILD 
traits_data_matseas = subset(traits_data, !is.na(MatingSeasDur))
nameck = name.check(tree, traits_data_matseas, data.names=traits_data_matseas$Species)
nameck
tree_matseas = drop.tip(tree, nameck$tree_not_data)

traits_data_matseas$MatingSeasDur = as.numeric(traits_data_matseas$MatingSeasDur)
hist(traits_data_matseas$MatingSeasDur)



#
pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # *p=0.0045
plot(mean_value_White.Blood.Cells~MatingSeasDur, data=traits_data_matseas)

traits_data_matseas$preds = predict(pgls3)
ggplot(data=traits_data_matseas, aes(x=MatingSeasDur, y=mean_value_White.Blood.Cells)) + 
  geom_smooth(method="lm", formula = y~x, aes(y=preds), color="#33a02c") +
  geom_point(color="#33a02c", size=2) +
  xlab("Mating season duration in wild (months)") +
  ylab("Mean white blood cell count (10^3/mm^3)") + theme_bw() 



#p = ggtree(tree) + geom_tiplab(size=2) + coord_cartesian(clip = 'off')
#p1 = facet_plot(p, panel="White blood cells", data=traits_data, geom=geom_point,
#                aes(x=scale(mean_value_White.Blood.Cells)), pch=1) + geom_point(data=traits_data, aes(x=scale(MatingSeasDur)), pch=2)
#facet_plot(p, panel="Mating season - wild", data=traits_data, geom=geom_point,
#           aes(x=scale(MatingSeasDur)), pch=2) +
#  theme(legend.position="none") 
#p1

#p1 = facet_plot(p, panel='WBC-blue, MatSeas-green (z-transf)', data=traits_data, geom=geom_point, 
#                aes(x=scale(mean_value_White.Blood.Cells)),pch=1, color='blue3') 
#p2 = facet_plot(p1, panel='WBC-blue, MatSeas-green (z-transf)', data=traits_data, geom=geom_point, 
#                aes(x=scale(MatingSeasDur)),pch=2, color='green') 
#p2 + theme_tree2(panel.spacing = unit(5.5, "lines"))    


pgls3 = gls(mean_value_Neutrophil.Seg.Abs~MatingSeasDur+ mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls3) #p=0.038, lambda negative

pgls3 = gls(mean_value_Neutrophil.Seg.Abs~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls3) #p=0.039


pgls3 = gls(mean_value_Lymphocytes.Abs~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls3)


pgls3 = gls(mean_value_Eosinophils.Abs~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls3) 

pgls3 = gls(mean_value_Basophils.Abs ~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs))
summary(pgls3) 

pgls3 = gls(mean_value_Monocytes.Abs~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls3) # * p = 0.0448


pgls3 = gls(mean_value_Red.Blood.Cells~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Red.Blood.Cells)) 
summary(pgls3) 


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
traits_data_directpara = traits_data[!is.na(traits_data$prop_direct),]
nameck = name.check(tree, traits_data_directpara, data.names=traits_data_directpara$Species)
nameck
tree_directpara=drop.tip(tree, nameck$tree_not_data)

pgls4 = gls(mean_value_White.Blood.Cells~prop_direct, data=traits_data_directpara, 
            correlation = corPagel(0.1, tree_directpara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls4) #lambda negative

pgls4 = gls(mean_value_White.Blood.Cells~prop_direct + mean_value_Body.Weight, data=traits_data_directpara, 
            correlation = corPagel(0, tree_directpara, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls4)

pgls4 = gls(mean_value_White.Blood.Cells~prop_direct + mean_value_Body.Weight, data=traits_data_directpara, 
            correlation = corPagel(0, tree_directpara, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls4) #body mass doesn't change ... 

# Figure ? 
ggplot(data=traits_data_directpara, aes(x=prop_direct, y=mean_value_White.Blood.Cells)) + 
  geom_point(color="#ef6548", size=2) +
  xlab("Directly transmitted parasites (proportion of total)") +
  ylab("Mean white blood cell count (10^3/mm^3)") + theme_bw() 


#including BM as covariate, it doesn't change... 
pgls4 = gls(mean_value_Basophils.Abs~prop_direct + mean_value_Body.Weight, data=traits_data_directpara, 
            correlation = corPagel(0.5, tree_directpara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs)) 
summary(pgls4) 

pgls4 = gls(mean_value_Eosinophils.Abs~prop_direct + mean_value_Body.Weight , data=traits_data_directpara, 
            correlation = corPagel(0.5, tree_directpara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs)) 
summary(pgls4) 

pgls4 = gls(mean_value_Lymphocytes.Abs~prop_direct + mean_value_Body.Weight, data=traits_data_directpara, 
            correlation = corPagel(0.5, tree_directpara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs)) 
summary(pgls4) #p=0.03

pgls4 = gls(mean_value_Monocytes.Abs~prop_direct + mean_value_Body.Weight, data=traits_data_directpara, 
            correlation = corPagel(0.5, tree_directpara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls4) 

pgls4 = gls(mean_value_Neutrophil.Seg.Abs~prop_direct + mean_value_Body.Weight, data=traits_data_directpara, 
            correlation = corPagel(0.5, tree_directpara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs)) 
summary(pgls4)  #body mass is significant ... lambda negative

pgls4 = gls(mean_value_Neutrophil.Seg.Abs~prop_direct + mean_value_Body.Weight, data=traits_data_directpara, 
            correlation = corPagel(0, tree_directpara, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs)) 
summary(pgls4)  #body mass is significant ... 



#parasites ~ wbcs
plot(prop_close~mean_value_White.Blood.Cells, data=traits_data)
plot(prop_direct~mean_value_White.Blood.Cells, data=traits_data)

plot(PSR_close~mean_value_White.Blood.Cells, data=traits_data)
plot(PSR_direct~mean_value_White.Blood.Cells, data=traits_data)



### proportion close transmitted

traits_data_closepara = traits_data[!is.na(traits_data$prop_close),]
nameck = name.check(tree, traits_data_closepara, data.names=traits_data_closepara$Species)
nameck
tree_closepara=drop.tip(tree, nameck$tree_not_data)

pgls5 = gls(mean_value_White.Blood.Cells~prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0.5, tree_closepara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls5) #lambda negative

pgls5 = gls(mean_value_White.Blood.Cells ~ prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0, tree_closepara, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells)) 
summary(pgls5) 


pgls5 = gls(mean_value_Basophils.Abs ~ prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0.5, tree_closepara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs)) 
summary(pgls5)#lambda negative

pgls5 = gls(mean_value_Basophils.Abs ~ prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0, tree_closepara, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs)) 
summary(pgls5)


pgls5 = gls(mean_value_Eosinophils.Abs ~ prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0.5, tree_closepara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs)) 
summary(pgls5) #p=0.03

pgls5 = gls(mean_value_Lymphocytes.Abs ~ prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0.5, tree_closepara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs)) 
summary(pgls5)

pgls5 = gls(mean_value_Monocytes.Abs~ prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0.5, tree_closepara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls5) 

pgls5 = gls(mean_value_Neutrophil.Seg.Abs~prop_close + mean_value_Body.Weight, data=traits_data_closepara, 
            correlation = corPagel(0.5, tree_closepara, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs)) 
summary(pgls5) 









