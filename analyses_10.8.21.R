library(tidyr)
library(dplyr)
library(ape)
library(nlme)
library(geiger)
library(ggplot2)
library(ggtree)
library(PerformanceAnalytics)
library(GGally)
library(data.table)
library(phytools)


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
pgls_bm = gls(mean_value_White.Blood.Cells*1000~mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls_bm) # p = 0.72; lambda = 0.55

#neutrophils ~ BM
pgls_bm = gls(mean_value_Neutrophil.Seg.Abs/1000~mean_value_Body.Weight, data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls_bm) # p = 0.009; lambda = -0.196... but lambda is negative... so I re-fit the model fixing lambda at 0.

pgls_bm = gls(mean_value_Neutrophil.Seg.Abs~mean_value_Body.Weight, data=traits_data, 
              correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
              method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls_bm) #p= 0.008

#lymphocytes ~ BM
pgls_bm = gls(mean_value_Lymphocytes.Abs~mean_value_Body.Weight, data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls_bm) #p= 0.08

#basophils ~ BM
pgls_bm = gls(mean_value_Basophils.Abs~mean_value_Body.Weight, data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls_bm) #p= 0.70, but lambda is negative, so re-fit?

pgls_bm = gls(mean_value_Basophils.Abs~mean_value_Body.Weight, data=traits_data, 
              correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
              method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls_bm) #p= 0.35

#eosinophils ~ BM
pgls_bm = gls(mean_value_Eosinophils.Abs~mean_value_Body.Weight, data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls_bm) #p= 0.38

#monocytes ~ BM
pgls_bm = gls(mean_value_Monocytes.Abs~mean_value_Body.Weight, data=traits_data, 
              correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
              method="REML", weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls_bm) #p= 0.16

#RBC ~ BM
pgls_bm = gls(scale(mean_value_Red.Blood.Cells)~factor(mean_value_Body.Weight), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE),
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls_bm)# *  p = 0.16

pgls_bm$modelStruct[1]
pgls_bm$parAssign



###########################
# look at correlations among the wbc types 

traits_data_cor.plot = traits_data[,c(3,5:8,10)]
rownames(traits_data_cor.plot) = traits_data[,1]
colnames(traits_data_cor.plot) = c("Basophils","Eosinophils","Lymphocytes","Monocytes","Neutrophils","TotalWBC")


lower_plots <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "black") +
    geom_smooth(method="gam",...)
}  


upper_plots <- function(data, mapping) {
  # extract x and y data from mapping
  x <- data[[rlang::get_expr(mapping[["x"]])]]
  y <- data[[rlang::get_expr(mapping[["y"]])]]

  # get phylo correlation
  data.mat = matrix(c(x,y), byrow=F,dimnames=NULL)
  rownames(data.mat) = rownames(traits_data_cor.plot)

  obj<-phyl.vcv(data.mat,vcv(tree),1)
  
 #####
   r.xy<-cov2cor(obj$R)[1,2]
  #####
  
  res <- round(r.xy,2)
  # plot weighted correlation as a label in a blank plot
  ggplot() + theme_void() + geom_text(aes(0,0,label=res))
}

ggpairs(data=traits_data_cor.plot,
  diag = list(continuous = wrap("densityDiag")),
  lower = list(continuous = wrap(lower_plots, se=F)),
  upper = list(continuous = wrap(upper_plots))
) 

x="Basophils"
rlang::get_expr(x)
mat = matrix(data=c(2,5,4,3,2,6), nrow=3, dimnames=list(c(),c("Basophils","Eosinophils")))
mat[,rlang::get_expr(x)]

rlang::get_expr(mapping)

rlang::get_expr(ggpairs(data=traits_data_cor.plot))

### Supp Fig S1 - not phylogenetically controlled
ggpairs(data=traits_data_cor.plot) + theme_bw()

x=c(5,6,20)
enquo(x)
quote(x)

reverse_mapping <- function(mapping) {
  aes_args <- paste(names(mapping), stringr::str_sub(as.character(mapping), start=2), sep = "=", collapse = ", ")
  aes_text <- glue::glue("aes({aes_args})")
  aes_text
}
reverse_mapping(ex_plot)

gsub("~","",y)
gsub("~","",x)

##############
# Supp Figure S2
traits_data_gat = traits_data %>% pivot_longer(c(mean_value_Basophils.Abs, mean_value_Eosinophils.Abs, 
                                           mean_value_Lymphocytes.Abs, mean_value_Monocytes.Abs, 
                                           mean_value_Neutrophil.Seg.Abs, mean_value_White.Blood.Cells),
                                           names_to = "WBCtype", values_to = "mean_value")
traits_data_gat$WBCtype = gsub("mean_value_","",traits_data_gat$WBCtype)
traits_data_gat$WBCtype = gsub(".Abs","",traits_data_gat$WBCtype)
traits_data_gat$WBCtype[traits_data_gat$WBCtype=="Neutrophil.Seg"] = "Neutrophils"
traits_data_gat$WBCtype[traits_data_gat$WBCtype=="White.Blood.Cells"] = "Total WBC"

ggplot(data=traits_data_gat, aes(x=mean_value_Body.Weight, y=mean_value)) + 
  facet_wrap(~WBCtype, scales="free") + geom_point() + theme_bw() +
  xlab("Body mass (kg)") + ylab("Mean value")




######################
#WBC ~ males per female
pgls1 = gls(scale(mean_value_White.Blood.Cells)~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1) # p = 0.20; lambda = 0.55

#Figure 1A
ggplot(data=traits_data, aes(x=factor(Males_per_female_CRA, labels=c("Multiple\n(n=20 spp)","Single\n(n=4 spp)")), y=mean_value_White.Blood.Cells, fill=factor(Males_per_female_CRA))) +
  geom_boxplot() + theme_bw() + ylab(expression(Mean~White~Blood~Cells~(10^3~"/"~mm^3))) + xlab("Sexual Partners") + 
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
pgls1 = gls(mean_value_White.Blood.Cells~factor(Males_per_female_CRA) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1) # p = 0.11; lambda = 0.55

#Neutrophils ~ males per female
pgls1 = gls(mean_value_Neutrophil.Seg.Abs~factor(Males_per_female_CRA), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs)) # false convergence

### ***neutrophils*** ###
pgls1 = gls(mean_value_Neutrophil.Seg.Abs~factor(Males_per_female_CRA) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) # p=0.13 * neutrophil no longer significant when control for body mass (pvalue for BM = 0.0490)


#Figure 1B
ggplot(data=traits_data, aes(x=factor(Males_per_female_CRA, labels=c("Multiple\n(n=20 spp)","Single\n(n=4 spp)")), 
                             y=mean_value_Neutrophil.Seg.Abs/1000, fill=factor(Males_per_female_CRA))) +
  geom_boxplot() + theme_bw() + ylab(expression(Mean~Neutrophils~(10^3~"/"~mm^3))) + xlab("Sexual Partners") + 
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
plot(mean_value_Body.Weight~MaleBodyMass.TestesDataset..in.g, data=traits_data_testes)


nameck = name.check(tree, traits_data_testes, data.names=traits_data_testes$Species)
nameck
tree_testes=drop.tip(tree, nameck$tree_not_data)




pgls_testes = gls(CombinedTestesMass.in.g~MaleBodyMass.TestesDataset..in.g, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE),
            method="REML")
summary(pgls_testes)

pgls_testes = gls(CombinedTestesMass.in.g~MaleBodyMass.TestesDataset..in.g, data=traits_data_testes, 
                  correlation = corPagel(0, tree_testes, form=~Species, fixed=TRUE),
                  method="REML")
summary(pgls_testes)

traits_data_testes$relative_testes_size_lp = residuals(pgls_testes) #rel testes size from luepold data
traits_data_testes$pred = predict(pgls_testes)

ggplot(data=traits_data_testes, aes(x=MaleBodyMass.TestesDataset..in.g, y=CombinedTestesMass.in.g)) +
  geom_point(size=2) + theme_bw() + xlab("Body mass (g)") + ylab("Combined testes mass (g)") + 
  geom_smooth(aes(y=pred), method="lm", formula=y~x, col="black")


pgls1 = gls(mean_value_White.Blood.Cells~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls1)

traits_data_testes$preds = predict(pgls1)

# Figure 2 #
ggplot(data=traits_data_testes, aes(x=relative_testes_size_lp, y=mean_value_White.Blood.Cells)) + 
  #geom_smooth(method="lm", formula = y~x, aes(y=preds), color="#8c96c6") +
  geom_point(color="#8c96c6", size=2) +
  xlab("Relative Testes Mass (Residuals)") +
  ylab(expression(Mean~White~Blood~Cells~(10^3~"/"~mm^3))) + theme_bw() 

pgls1 = gls(mean_value_Neutrophil.Seg.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls1) 

pgls1 = gls(mean_value_Lymphocytes.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Lymphocytes.Abs)) 
summary(pgls1) #lambda > 1 !

pgls1 = gls(mean_value_Lymphocytes.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(1, tree_testes, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls1)


pgls1 = gls(mean_value_Basophils.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Basophils.Abs))
summary(pgls1)

pgls1 = gls(mean_value_Monocytes.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls1) # lambda negative

pgls1 = gls(mean_value_Monocytes.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0, tree_testes, form=~Species, fixed=TRUE), 
            method="REML", weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls1) # lambda negative


pgls1 = gls(mean_value_Eosinophils.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls1)

pgls1 = gls(mean_value_Eosinophils.Abs~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls1)

pgls1 = gls(mean_value_Red.Blood.Cells~relative_testes_size_lp, data=traits_data_testes, 
            correlation = corPagel(0.5, tree_testes, form=~Species, fixed=FALSE), 
            method="REML", weights=~I(1/n_ind_Red.Blood.Cells))
summary(pgls1) # p = 0.14

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

pgls2 = gls(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0002, but lambda is negative

pgls2 = gls(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding), data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0004 with fixed lambda = 0

#same pattern controlling for body mass too 
pgls2 = gls(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding) + mean_value_Body.Weight , data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0002



#Figure 3

alternate_names = gsub("_"," ", tree$tip.label)
name_frame = data.frame(cbind(tree$tip.label,alternate_names))
names(name_frame) = c("label","alternate_names")
name_frame$label=as.character(name_frame$label)
name_frame$alternate_names=as.character(name_frame$alternate_names)

treeplot = ggtree(tree) %<+% name_frame

traits_data$MatingSeasDur=as.numeric(traits_data$MatingSeasDur)

p = treeplot + coord_cartesian(clip = 'off')
p1 = facet_plot(p, panel="White Blood Cells", data=traits_data, geom=geom_point, size=1.5,
                aes(x=mean_value_White.Blood.Cells, fill=MatingSeasDur, shape=factor(R_Pattern_Breeding))) +
  #theme(legend.position="none") + 
  scale_fill_gradient(low="lightblue", high="darkblue", na.value="white",name="Mating season duration") +
  scale_shape_manual(values=c(24,21), name="Seasonal in\n captivity?",labels=c("No","Yes"))

p1 + theme_tree2(panel.spacing = unit(7, "lines"))  + 
  geom_tiplab(aes(label=alternate_names), 
               size=2.5)


pgls2 = gls(mean_value_Lymphocytes.Abs~factor(R_Pattern_Breeding) , data=traits_data, 
            correlation = corPagel(0.7, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0004


pgls2 = gls(mean_value_Lymphocytes.Abs/1000~factor(R_Pattern_Breeding) + mean_value_Body.Weight , data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0003 - same pattern controlling for body mass; lambda negative

pgls2 = gls(mean_value_Lymphocytes.Abs/1000~factor(R_Pattern_Breeding) + mean_value_Body.Weight , data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0004


pgls2 = gls(scale(mean_value_Neutrophil.Seg.Abs)~factor(R_Pattern_Breeding)+ mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls2)

pgls2 = gls(mean_value_Eosinophils.Abs~factor(R_Pattern_Breeding) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls2)

pgls2 = gls(mean_value_Monocytes.Abs/1000~factor(R_Pattern_Breeding) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs))
summary(pgls2) #0.009

pgls2 = gls(mean_value_Basophils.Abs~factor(R_Pattern_Breeding)+ mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0.5, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs))
summary(pgls2) #lambda neg

pgls2 = gls(mean_value_Basophils.Abs~factor(R_Pattern_Breeding)+ mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs))
summary(pgls2) #lambda 0


pgls2 = gls(mean_value_Red.Blood.Cells~factor(R_Pattern_Breeding) + mean_value_Body.Weight, data=traits_data, 
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


### Mating season duration
#
pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # *p=0.0045, lambda negative
plot(mean_value_White.Blood.Cells~MatingSeasDur, data=traits_data_matseas)

pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # *0.014
plot(mean_value_White.Blood.Cells~MatingSeasDur, data=traits_data_matseas)


traits_data_matseas$preds = predict(pgls3)

# Figure 4
ggplot(data=traits_data_matseas, aes(x=MatingSeasDur, y=mean_value_White.Blood.Cells)) + 
  geom_smooth(method="lm", formula = y~x, aes(y=preds), color="#33a02c") +
  geom_point(color="#33a02c", size=2) +
  xlab("Mating Season Duration in Wild (Months)") +
  ylab(expression(Mean~White~Blood~Cells~(10^3~"/"~mm^3))) + theme_bw() 



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


pgls3 = gls(mean_value_Neutrophil.Seg.Abs/1000~MatingSeasDur+ mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls3) #p=0.038, lambda negative

pgls3 = gls(mean_value_Neutrophil.Seg.Abs/1000~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_Neutrophil.Seg.Abs))
summary(pgls3) #p=0.039


pgls3 = gls(mean_value_Lymphocytes.Abs/1000 ~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls3)


pgls3 = gls(mean_value_Eosinophils.Abs/1000~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Eosinophils.Abs))
summary(pgls3) 

pgls3 = gls(mean_value_Basophils.Abs/1000 ~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Basophils.Abs))
summary(pgls3) 

pgls3 = gls(mean_value_Monocytes.Abs/1000~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Monocytes.Abs)) 
summary(pgls3) # * p = 0.0448


pgls3 = gls(mean_value_Red.Blood.Cells~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Red.Blood.Cells)) 
summary(pgls3) 



############################### new december 2021
### mating season duration + mating system + body mass
traits_data_matseas$Males_per_female_CRA 
pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + Males_per_female_CRA +mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0.5, tree_matseas, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # *p=0.0045, lambda negative
plot(mean_value_White.Blood.Cells~MatingSeasDur, data=traits_data_matseas)

pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + Males_per_female_CRA +mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # *0.014
plot(mean_value_White.Blood.Cells~MatingSeasDur, data=traits_data_matseas)

####################








#parasites ~ mating season duration
plot(prop_close~MatingSeasDur, data=traits_data)
plot(prop_direct~MatingSeasDur, data=traits_data)

plot(PSR_close~MatingSeasDur, data=traits_data)
plot(PSR_direct~MatingSeasDur, data=traits_data)
summary(lm(PSR_direct~MatingSeasDur, data=traits_data))







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

# Figure 5A
ggplot(data=traits_data_directpara, aes(x=prop_direct, y=mean_value_White.Blood.Cells)) + 
  geom_point(color="#ef6548", size=2) +
  xlab("Directly Transmitted Parasites") +
  ylab(expression(Mean~White~Blood~Cells~(10^3~"/"~mm^3))) + theme_bw() 


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




# Figure 5B
ggplot(data=traits_data_closepara, aes(x=prop_close, y=mean_value_White.Blood.Cells)) + 
  geom_point(color="#ef6548", size=2) +
  xlab("Close-Contact Transmitted Parasites") +
  ylab(expression(Mean~White~Blood~Cells~(10^3~"/"~mm^3))) + theme_bw() 



#Supp Figure - S5

traits_data_2 = traits_data %>% pivot_longer(cols= c(prop_close, prop_direct), names_to = "Transmission", values_to="Proportion")


ggplot(data=traits_data_2, aes(x=MatingSeasDur, y=Proportion)) +
  facet_grid(~factor(Transmission, labels=c("Close Transmission", "Direct Life Cycle"))) +
  geom_point() + theme_bw() + xlab("Mating Season Duration (Months)") + ylab("Proportion of Total Parasites") 

ggplot(data=traits_data_2, aes(x=MatingSeasDur, y=prop_direct)) +
  geom_point() + theme_bw()


traits_data_closepara_MS = subset(traits_data_closepara, !is.na(MatingSeasDur))
nameck = name.check(tree, traits_data_closepara_MS, data.names=traits_data_closepara_MS$Species)
tree_closepara_MS=drop.tip(tree, nameck$tree_not_data)


pgls.MP1 = gls(prop_close ~ MatingSeasDur, data=traits_data_closepara_MS, 
               correlation = corPagel(0.5, tree_closepara_MS, form=~Species, fixed=FALSE), method="REML") 
summary(pgls.MP1) 



traits_data_directpara_MS = subset(traits_data_directpara, !is.na(MatingSeasDur))
nameck = name.check(tree, traits_data_directpara_MS, data.names=traits_data_directpara_MS$Species)
tree_directpara_MS=drop.tip(tree, nameck$tree_not_data)

pgls.MP2 = gls(prop_direct ~ MatingSeasDur, data=traits_data_directpara_MS, 
               correlation = corPagel(0.5, tree_directpara_MS, form=~Species, fixed=FALSE), method="REML") 
summary(pgls.MP2) 







## ***** Supplement ***** #####
traits_data_testes$Males_per_female_CRA

#Figure S4
pgls0 = gls(relative_testes_size_lp~Males_per_female_CRA, data=traits_data_testes, 
            correlation = corPagel(0, tree_testes, form=~Species, fixed=TRUE), method="REML")
summary(pgls0) #not significant, but a trend toward promiscuous spp having larger testes

ggplot(data=traits_data_testes, aes(x=Males_per_female_CRA, y=relative_testes_size_lp)) +
  geom_boxplot(fill="gray") + theme_bw() + xlab("Mating partners") + ylab("Relative testes size (residual)")


traits_data_testes2 = subset(traits_data_testes, !is.na(MatingSeasDur))

nameck = name.check(tree_testes, traits_data_testes2, data.names=traits_data_testes2$Species)
nameck
tree_testes2=drop.tip(tree_testes, nameck$tree_not_data)

pgls0 = gls(relative_testes_size_lp~MatingSeasDur, data=traits_data_testes2, 
            correlation = corPagel(0.5, tree_testes2, form=~Species, fixed=FALSE), method="REML")
summary(pgls0) 

traits_data_testes$MatingSeasDur




#breeding seasonality data from Heldstab et al. paper
held_data = read.csv("Heldstab et al. data.csv", header=T)

held_data$Species=gsub(" ","_",held_data$Species)
held_data$Species=gsub("\\*","",held_data$Species)

traits_data_1 = traits_data %>% left_join(held_data) %>% filter(!is.na(Seasonality.natural.habitat) & !is.na(MatingSeasDur))

plot(MatingSeasDur~jitter(Seasonality.natural.habitat), traits_data_1)
traits_data_1$Seasonality.natural.habitat


season=lm(MatingSeasDur~Seasonality.natural.habitat, data=traits_data_1)
summary(season)
traits_data_1$predictions = NA
traits_data_1$predictions = predict(season)
length(predict(season))

ggplot(data=traits_data_1, aes(x=Seasonality.natural.habitat, y=MatingSeasDur))+
  theme_bw() + geom_jitter(width=0.1, height=0.1) + xlab("Seasonality classification - Heldstab et al. 2021") +
  ylab("Mating Season Duration (months) - this study") + geom_smooth(aes(y=predictions),method="lm", formula=y~x, col="black")


max_age_pad = c(29.6, 22.39, 34.56, 28.91, 29.67, 27.38, 29.81, 34.45, 33.66, 29.24, 33.38,
                32.92, 14.07, 23.49, 32.7, 21.52, 17.96, 18.04, 20.5, 19.25, 30.57, 25.41, 35.08, 35.73)

plot(max_age_pad~L_Median_All_Longevity_gt30d_y, data=traits_data)

