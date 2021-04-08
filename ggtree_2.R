#ok wow the installation for ggtree is weird and intense
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ggtree")

library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(nlme)

#read in data and phylogeny

data1 = read.csv("cumulative_data_final.csv", header=T)
str(data1)
colnames(data1) = c("X", "species", "Scientific_Name","Mating_Season_Duration","Mating_System","All_WBC","Body_Mass","Lymphocytes","Monocytes",
                    "Basophils","Eosinophils","Neutrophils","Monocytes_Nind","Lymphocytes_Nind",
                    "Basophils_Nind","Eosinophils_Nind","Neutrophils_Nind", "Body_Mass_Nind", "All_WBC_Nind")
data = data1[,-c(1:2)]


tree = read.nexus("consensusTree_10kTrees_Primates_allstreps.nex")
plot(tree)
str(tree$tip.label)

tree_plot = ggtree(tree) + geom_tiplab() 

tree_plot + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6))


# Make the original tree plot
p = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')


# BODY MASS 
p2 = facet_plot(p, panel="Body mass (kg)", data=data, geom=geom_point, 
                 aes(x=Body_Mass, color = Mating_System)) 

p2 + theme_tree2(panel.spacing = unit(5, "lines"))

BM = setNames(data$Body_Mass,
              data$Scientific_Name)

phylosig(tree, BM, method="lambda")
phylosig(tree, BM, method="K")


# ALL WBC
p3 = facet_plot(p, panel='All WBCs', data=data, geom=geom_segment, 
                 aes(x=0, xend=All_WBC, y=y, yend=y), size=3, color='blue4') 

p3 + theme_tree2(panel.spacing = unit(5, "lines"))


AllWBCs = setNames(data$All_WBC,
                   data$Scientific_Name)

phylosig(tree, AllWBCs, method="lambda")
phylosig(tree, AllWBCs, method="K")


# Lymphocytes
p4 = facet_plot(p, panel='Lymphocytes', data=data, geom=geom_segment, 
                aes(x=0, xend=Lymphocytes, y=y, yend=y), size=3, color='blue4') 

p4 + theme_tree2(panel.spacing = unit(5, "lines"))


Lymphs = setNames(data$Lymphocytes,
                   data$Scientific_Name)

phylosig(tree, Lymphs, method="lambda")
phylosig(tree, Lymphs, method="K")


# look at correlations 

rownames(data) = data$Scientific_Name
fancyTree(tree,type=c("scattergram"),X=data[,c(4:6)], label="off")

####
# start with body mass-lymphocytes
par(mfrow=c(1,1))
phylomorphospace(tree,data[,c(5,6)], label="off", node.size=c(0,1.2))

ols.model = lm(Lymphocytes~BM, data)
summary(ols.model)
abline(ols.model, col="red")
data$ols_resid= resid(ols.model)

pgls.model.lam = gls(Lymphocytes~BM, data, correlation=corPagel(1, tree, form = ~Scientific_Name, fixed=F), method="ML")
summary(pgls.model.lam)
abline(pgls.model.lam, col="blue")

#pgls.model.brown = gls(Lymphocytes~BM, data, correlation=corBrownian(1, tree, form = ~Scientific_Name), method="ML")
#summary(pgls.model.brown)
#abline(pgls.model.brown, col="green")



# Make the original tree plot
myplot = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')


# TRAITS AND RESIDUALS

myplot1 = facet_plot(myplot, panel='Trait', data=data, geom=geom_segment, 
                aes(x=0, xend=Lymphocytes, y=y, yend=y), size=3, color='blue3') 
myplot2 = facet_plot(myplot1, panel="Trait", data=data, geom=geom_point, 
                aes(x=Body_Mass*20 ), col="orange")

myplot2 + theme_tree2(panel.spacing = unit(5, "lines"))

myplot3 = facet_plot(myplot2, panel="OLS residuals", data=data, geom=geom_point, 
                     aes(x=ols_resid), col="black")
myplot3 + theme_tree2(panel.spacing = unit(c(5,1), "lines"))

myplot4 = facet_plot(myplot, panel="OLS residuals", data=data, geom=geom_point, 
                     aes(x=ols_resid), col="black")
myplot4 + theme_tree2(panel.spacing = unit(5, "lines"))


phylomorphospace(tree, data[,c(5,6)], label="off", xlab="Body Mass (kg)", ylab="Lymphocytes",
                 node.size=c(0,1.2))


OLS_Resids = setNames(data$ols_resid,
                      data$Scientific_Name)

phylosig(tree, OLS_Resids, method="lambda")
phylosig(tree, OLS_Resids, method="K")



####
# now look at All WBCs - BM
par(mfrow=c(1,1))
phylomorphospace(tree,data[,c(5,4)], label="off", node.size=c(0,1.2))

ols.model = lm(All_WBC~BM, data)
summary(ols.model)
abline(ols.model, col="red")
data$ols_resid= resid(ols.model)

pgls.model.lam = gls(All_WBC~BM, data, correlation=corPagel(1, tree, form = ~Scientific_Name, fixed=F), method="ML")
summary(pgls.model.lam)
abline(pgls.model.lam, col="blue")

#pgls.model.brown = gls(All_WBC~BM, data, correlation=corBrownian(1, tree, form = ~Scientific_Name), method="ML")
#summary(pgls.model.brown)
#abline(pgls.model.brown, col="green")



# Make the original tree plot
myplot = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')


# TRAITS AND RESIDUALS

myplot1 = facet_plot(myplot, panel='Trait', data=data, geom=geom_segment, 
                     aes(x=0, xend=All_WBC, y=y, yend=y), size=3, color='blue3') 
myplot2 = facet_plot(myplot1, panel="Trait", data=data, geom=geom_point, 
                     aes(x=Body_Mass*2), col="orange")

myplot2 + theme_tree2(panel.spacing = unit(5, "lines"))

myplot3 = facet_plot(myplot2, panel="OLS residuals", data=data, geom=geom_point, 
                     aes(x=ols_resid), col="black")
myplot3 + theme_tree2(panel.spacing = unit(c(5,1), "lines"))

myplot4 = facet_plot(myplot, panel="OLS residuals", data=data, geom=geom_point, 
                     aes(x=ols_resid), col="black")
myplot4 + theme_tree2(panel.spacing = unit(5, "lines"))


phylomorphospace(tree, data[,c(5,6)], label="off", xlab="Body Mass (kg)", ylab="All WBC",
                 node.size=c(0,1.2))


OLS_Resids = setNames(data$ols_resid,
                  data$Scientific_Name)

phylosig(tree, OLS_Resids, method="lambda")
phylosig(tree, OLS_Resids, method="K")

