
library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(nlme)
library(caper)

data= read.csv("cumulative_data_finalll.csv", header=T)
str(data)
data=data[,-1]
tree=read.nexus("consensusTree_10kTrees_Primates_allstreps.nex")

#Ordinary least squares regressions - if you scale y and x variables, the outcome is
#the same slope, regardless of which var is dependent vs independent
ols1 = lm(scale(Lymphocytes)~scale(BodyWeight), data=data)
summary(ols1)
plot(scale(Lymphocytes)~scale(BodyWeight), data=data,  xlim=c(-2.1,2.1), ylim=c(-2.1,2.1))
abline(ols1)

ols2 = lm(scale(BodyWeight)~scale(Lymphocytes), data=data)
summary(ols2)
plot(scale(BodyWeight)~scale(Lymphocytes), data=data, xlim=c(-2.1,2.1), ylim=c(-2.1,2.1))
abline(ols2)

comp_dat = comparative.data(tree, data, names="Scientific_Name")

#but that doesn't hold up when you scale the variables in a phylogenetic analysis
#(when you allow the model to estimate lambda; fixed = FALSE)
#this is using ape

pgls1 = gls(scale(BodyWeight)~scale(Lymphocytes), data=data, correlation = corPagel(value=0.7, tree, form=~Scientific_Name, fixed=FALSE), method="REML")
summary(pgls1) #lambda = 0.913

pgls2 = gls(scale(Lymphocytes)~scale(BodyWeight), data=data, correlation = corPagel(value=0.7, tree, form=~Scientific_Name, fixed=FALSE), method="REML")
summary(pgls2) #using maximum likelihood, I get that lambda=-3.35. I THINK negative lambda means closely related species are much more different than expected by brownian motion 
#using method="REML" gives lambda around 0.6, I think this prevents lambda from going below zero? 

#this is just a fancier way of running an ord least sq reg; this is the same as above. Note: value=0, fixed=TRUE
ols1 = gls(scale(BodyWeight)~scale(Lymphocytes), data=data, correlation = corPagel(value=0, tree, form=~Scientific_Name, fixed=TRUE))
summary(ols1)
plot(resid(ols1))

ols2 = gls(scale(Lymphocytes)~scale(BodyWeight), data=data, correlation = corPagel(value=0, tree, form=~Scientific_Name, fixed=TRUE))
summary(ols2)
plot(resid(ols2))

#look at phylogenetic signal in the residuals for the two "versions" of the model
data$resids_ybm = resid(ols1)[,1]
resids_ybm=data$resids_ybm
names(resids_ybm) = data$Scientific_Name

data$resids_ylymph = resid(ols2)[,1]
resids_ylymph=data$resids_ylymph
names(resids_ylymph) = data$Scientific_Name

phylosig(tree, x=resids_ybm, method="lambda") #lambda=0.714
phylosig(tree, x=resids_ylymph, method="lambda") #lambda = 0


p = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')
p2 = facet_plot(p, panel="BM~Lymph Resid", data=data, geom=geom_point, 
                aes(x=resids_ybm)) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel="Lymph~BM Resid", data=data, geom=geom_point, 
                aes(x=resids_ylymph)) 

p2 + theme_tree2(panel.spacing = unit(5, "lines"))


#does playing around with lambda give more similar results across the two models? Not entirely sure what to think of this... 
pgls5 = gls(scale(Lymphocytes)~scale(BodyWeight), data=data, correlation = corPagel(0.7, tree, form=~Scientific_Name, fixed=TRUE))
summary(pgls5)

pgls6 = gls(scale(BodyWeight)~scale(Lymphocytes), data=data, correlation = corPagel(0.7, tree, form=~Scientific_Name, fixed=TRUE))
summary(pgls6) 





