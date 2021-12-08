traits_data = read.csv("trait_data_strepsirrhines.csv", header=T)
traits_data$Species

traits_data$Nocturnal_CRA = c(1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,1,0,0,0,0)



# WBC ~ duration of mating season - in the WILD 
traits_data_matseas = subset(traits_data, !is.na(MatingSeasDur))
nameck = name.check(tree, traits_data_matseas, data.names=traits_data_matseas$Species)
nameck
tree_matseas = drop.tip(tree, nameck$tree_not_data)

pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=TRUE ), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # *0.014

pgls3 = gls(mean_value_Lymphocytes.Abs~MatingSeasDur + factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data_matseas, 
            correlation = corPagel(0, tree_matseas, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls3) # *0.014


traits_data_matseas2 = subset(traits_data_matseas, Nocturnal_CRA==0)
nameck = name.check(tree, traits_data_matseas2, data.names=traits_data_matseas2$Species)
nameck
tree_matseas2 = drop.tip(tree, nameck$tree_not_data)

pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + mean_value_Body.Weight, data=traits_data_matseas2, 
            correlation = corPagel(0, tree_matseas2, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # 


pgls2 = gls(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding)+ factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # 

pgls2 = gls(mean_value_White.Blood.Cells~factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # *p=0.0004 with fixed lambda = 0


pgls2 = gls(mean_value_Lymphocytes.Abs~factor(R_Pattern_Breeding)+ factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0004 with fixed lambda = 0



pgls2 = gls(mean_value_Lymphocytes.Abs~ factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data, 
            correlation = corPagel(0, tree, form=~Species, fixed=FALSE), method="REML",  
            weights=~I(1/n_ind_Lymphocytes.Abs))
summary(pgls2) # *p=0.0004 with fixed lambda = 0






#restrict to just madagascar
traits_data$Mada = c(1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,0,0,1,1,1,1)
cbind(traits_data$Species,traits_data$Mada)

traits_data_md = subset(traits_data_matseas, Mada=="1")
nameck = name.check(tree, traits_data_md, data.names=traits_data_md$Species)
nameck
tree_md = drop.tip(tree, nameck$tree_not_data)


pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + mean_value_Body.Weight, data=traits_data_md, 
            correlation = corPagel(0, tree_md, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # 

pgls3 = gls(mean_value_White.Blood.Cells~MatingSeasDur + factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data_md, 
            correlation = corPagel(0, tree_md, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls3) # 




traits_data_md2 = subset(traits_data, Mada=="1")
nameck = name.check(tree, traits_data_md2, data.names=traits_data_md2$Species)
nameck
tree_md2 = drop.tip(tree, nameck$tree_not_data)

pgls2 = gls(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding) + mean_value_Body.Weight, data=traits_data_md2, 
            correlation = corPagel(0, tree_md2, form=~Species, fixed=TRUE), method="REML",  
            weights=~I(1/n_ind_White.Blood.Cells))
summary(pgls2) # 




#Type III error - to address collinearity issue
#I haven't weighted the data here by sample size... not entirely sure how to do it analogously to above
#library(car)

#############
#just madagascar lemurs
length(traits_data_md$Species)
plot(MatingSeasDur~mean_value_White.Blood.Cells,data=traits_data_md)
plot(Nocturnal_CRA~mean_value_White.Blood.Cells,data=traits_data_md)
plot(mean_value_Body.Weight~mean_value_White.Blood.Cells,data=traits_data_md)


#only mating season and nocturnality
model = lm(mean_value_White.Blood.Cells~MatingSeasDur + factor(Nocturnal_CRA), data=traits_data)
anova(model) #type 1
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~factor(Nocturnal_CRA) + MatingSeasDur, data=traits_data)
anova(model) #type 1
Anova(model, type="III") 

#control for body mass 
model = lm(Nocturnal_CRA ~ mean_value_Body.Weight, data=traits_data_md)
model = lm(mean_value_White.Blood.Cells ~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_md)
model = lm(mean_value_White.Blood.Cells ~ factor(Nocturnal_CRA), data=traits_data)


model = lm(mean_value_White.Blood.Cells~MatingSeasDur + factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data_md)
anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~factor(Nocturnal_CRA) + MatingSeasDur + mean_value_Body.Weight, data=traits_data_md)
#anova(model) 
Anova(model, type="III") 

#enter body mass first?
model = lm(mean_value_White.Blood.Cells~ mean_value_Body.Weight + factor(Nocturnal_CRA) + MatingSeasDur , data=traits_data_md)
#anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~ mean_value_Body.Weight + MatingSeasDur + factor(Nocturnal_CRA) , data=traits_data_md)
#anova(model) 
Anova(model, type="III") 

#just mat seas dur and body mass
model = lm(mean_value_White.Blood.Cells~ MatingSeasDur + mean_value_Body.Weight, data=traits_data_md)
#anova(model)        
Anova(model, type="III")

#just nocturnal and body mass
model = lm(mean_value_White.Blood.Cells~ factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data_md)
#anova(model)        
Anova(model, type="III")


########
#all strepsirrhines - mating season duration
length(traits_data$Species)
plot(MatingSeasDur~mean_value_White.Blood.Cells,data=traits_data)
plot(Nocturnal_CRA~mean_value_White.Blood.Cells,data=traits_data)
plot(mean_value_Body.Weight~mean_value_White.Blood.Cells,data=traits_data)


#only mating season and nocturnality
model = lm(mean_value_White.Blood.Cells~MatingSeasDur + factor(Nocturnal_CRA), data=traits_data)
#anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~factor(Nocturnal_CRA) + MatingSeasDur, data=traits_data)
#anova(model) 
Anova(model, type="III") 

#control for body mass 
model = lm(mean_value_White.Blood.Cells~MatingSeasDur + factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data)
#anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~factor(Nocturnal_CRA) + MatingSeasDur + mean_value_Body.Weight, data=traits_data)
#anova(model) 
Anova(model, type="III") 

#enter body mass first?
model = lm(mean_value_White.Blood.Cells~ mean_value_Body.Weight + factor(Nocturnal_CRA) + MatingSeasDur , data=traits_data)
#anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~ mean_value_Body.Weight + MatingSeasDur + factor(Nocturnal_CRA) , data=traits_data)
#anova(model) 
Anova(model, type="III") 

#just mat seas dur and body mass
model = lm(mean_value_White.Blood.Cells~ MatingSeasDur + mean_value_Body.Weight, data=traits_data)
#anova(model)        
Anova(model, type="III")

#just nocturnal and body mass
model = lm(mean_value_White.Blood.Cells~ factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data)
#anova(model)        
Anova(model, type="III")


#########
#all strepsirrhines - seasonal breeding in captivity (binary)
length(traits_data$Species)


#only mating season and nocturnality
model = lm(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding) + factor(Nocturnal_CRA), data=traits_data)
#anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~factor(Nocturnal_CRA) + factor(R_Pattern_Breeding), data=traits_data)
#anova(model) 
Anova(model, type="III") 

#control for body mass 
model = lm(mean_value_White.Blood.Cells~factor(R_Pattern_Breeding) + factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data)
#anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~factor(Nocturnal_CRA) + factor(R_Pattern_Breeding) + mean_value_Body.Weight, data=traits_data)
#anova(model) 
Anova(model, type="III") 

#enter body mass first?
model = lm(mean_value_White.Blood.Cells~ mean_value_Body.Weight + factor(Nocturnal_CRA) + factor(R_Pattern_Breeding) , data=traits_data)
#anova(model) 
Anova(model, type="III") 

model = lm(mean_value_White.Blood.Cells~ mean_value_Body.Weight + factor(R_Pattern_Breeding) + factor(Nocturnal_CRA) , data=traits_data)
#anova(model) 
Anova(model, type="III") 

#just mat seas dur and body mass
model = lm(mean_value_White.Blood.Cells~ factor(R_Pattern_Breeding) + mean_value_Body.Weight, data=traits_data)
#anova(model)        
Anova(model, type="III")

#just nocturnal and body mass
model = lm(mean_value_White.Blood.Cells~ factor(Nocturnal_CRA) + mean_value_Body.Weight, data=traits_data)
#anova(model)        
Anova(model, type="III")




