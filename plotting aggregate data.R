

library(ggplot2)
library(ggtree)
library(ape)
library(phytools)
library(nlme)
library(caper)

ageslope2 = read.csv("findataabslopeallwbctypes.csv", header=T)
ageslope=ageslope2[,2:8]
tree = read.nexus("consensusTree_10kTrees_Primates_allstreps.nex")
traits = read.csv("cumulative_data_finalll.csv", header=T)

traits2 = traits[,c(2:4,6)]
age_data = merge(traits2, ageslope, by="Scientific_Name")


p = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')

p2 = facet_plot(p, panel='Lymphocytes slope with age', data=age_data, geom=geom_point, 
                aes(x=ABSlopelymph, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Basophils slope with age', data=age_data, geom=geom_point, 
                aes(x=ABSlopebaso, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Monocytes slope with age', data=age_data, geom=geom_point, 
                aes(x=ABSlopemono, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Neutophils slope with age', data=age_data, geom=geom_point, 
                aes(x=ABSlopeneutro, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Eosinophils slope with age', data=age_data, geom=geom_point, 
                aes(x=ABSlopeeos, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Total WBC slope with age', data=age_data, geom=geom_point, 
                aes(x=ABSlopewbc, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


sexdim2 = read.csv("sexual dimorphism for all species and all wbc types.csv", header=T)
sexdim=sexdim2[,c(3,2,4:9)]
colnames(sexdim)[1] = "Scientific_Name"

traits2 = traits[,c(2:4,6)]
dim_data = merge(traits2, sexdim, by="Scientific_Name")

p2 = facet_plot(p, panel='Lymphocytes sexual dimorphism', data=dim_data, geom=geom_point, 
                aes(x=LYMPHOCYTEMaledivAverage_Ratio, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Basophils sexual dimorphism', data=dim_data, geom=geom_point, 
                aes(x=BASOPHILMaledivAverage_Ratio, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Monocytes sexual dimorphism', data=dim_data, geom=geom_point, 
                aes(x=MONOCYTEMaledivAverage_Ratio, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Neutophils sexual dimorphism', data=dim_data, geom=geom_point, 
                aes(x=NEUTROPHILMaledivAverage_Ratio, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Eosinophils sexual dimorphism', data=dim_data, geom=geom_point, 
                aes(x=EOSINOPHILMaledivAverage_Ratio, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


p2 = facet_plot(p, panel='Total WBC sexual dimorphism', data=dim_data, geom=geom_point, 
                aes(x=ALLWBCMaledivAverage_Ratio, color=Mating_System), size=3) 
p2 + theme_tree2(panel.spacing = unit(5, "lines"))


#testes mass
#categorized mating system
#seasonality

##### replicate nunn study
#get papers for leila & read one together?
#analyses with wbc types - absolute counts
#dimorphism in wbc types - leila look at literature
#correlations lymphocytes neutrophils
#age??? 

pgls()




