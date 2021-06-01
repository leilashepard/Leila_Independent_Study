# EVERYTHING I ACTUALLY USED IN THE PAPER SIMPLIFIED


# Figure 1- chart.correlation
cumulative_data_finalll = read.csv("cumulative_data_finalll.csv", header=TRUE)
chart.Correlation(cumulative_data_finalll[, c(3, 5, 6, 7, 8, 9, 10, 11)], histogram=TRUE, pch=19)

# Figures 2, 3, 4, and 5
treeallstreps = read.nexus("consensusTree_10kTrees_Primates_allstreps.nex") #read in the tree file
lemurcomp = comparative.data(treeallstreps, cumulative_data_finalll, 'Scientific_Name')


    # did this for all trait correlations - Figure 2
pgls.model1 = pgls(AllWBC~Mating_System, data=lemurcomp, lambda="ML")
summary(pgls.model1)

    # did this for all traits - Figure 3
p = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')
p7 = facet_plot(p, panel="Body Weight", data=cumulative_data_finalll, geom=geom_point,
                aes(x=BodyWeight, size=body_weight_n_ind)) 

p7 + theme_tree2(panel.spacing = unit(5, "lines"))

BM = setNames(cumulative_data_finalll$BodyWeight,
              cumulative_data_finalll$Scientific_Name)

phylosig(tree, BM, method="lambda", test=TRUE)
phylosig(tree, BM, method="K")




# Figures 5 and 6 - lymphocyte conc vs. age and wbc conc vs. age
wbc_age = read.csv("wbc all species.csv")
wbc_age1 = subset(wbc_age, select = -c(birthdate, deathdate, age_at_death))
wbc_age_pivot_longer1 = pivot_longer(data=wbc_age1, cols=starts_with("y"), names_to="age", values_to = "bloodcellconcentration")
wbc_age_pivot_longer = wbc_age_pivot_longer1[-which(is.na(wbc_age_pivot_longer1$bloodcellconcentration)),]
unique_species = unique(wbc_age_pivot_longer$species)
wbc_age_pivot_longer$age = sub("y", "", wbc_age_pivot_longer$age)
wbc_age_pivot_longer3 = merge(wbc_age_pivot_longer, cumulative_data_finall[, c("species", "Mating_System")], by="species", all.x=FALSE)
wbc_age_pivot_longer2 = merge(wbc_age_pivot_longer3, cumulative_data_finall[,c("species", "Mating_Season_Duration")], by="species", all.x=FALSE)
wbc_age_pivot_longer1 = merge(wbc_age_pivot_longer2, cumulative_data_finall[,c("species", "wbc_n_ind")], by="species", all.x=FALSE)
wbc_age_pivot_longer = merge(wbc_age_pivot_longer1, cumulative_data_finall[,c("species", "Scientific_Name")], by="species", all.x=FALSE)
count = ddply(wbc_age_pivot_longer,.(age,species),nrow)
length(count$V1)
wbc_age_pivot_longer = wbc_age_pivot_longer[!duplicated(wbc_age_pivot_longer), ]
wbc_age_pivot_longer = wbc_age_pivot_longer %>% group_by(species, age) %>% mutate(wbc_n_ind = n())
final_table = wbc_age_pivot_longer %>% group_by(age, species) %>%
  summarize(Mating_System, Scientific_Name, Mating_Season_Duration, wbc_n_ind, mean_bcc = mean(bloodcellconcentration, na.rm=TRUE))
final_table = final_table[!duplicated(final_table), ]
final_table$age = strtoi(final_table$age)
final_table = final_table[order(final_table[,1] ),]
reg = lm(mean_bcc~age*Mating_System, data=final_table, weights=wbc_n_ind)
summary(reg)
final_table$regpredict = predict.lm(reg, final_table, weights=wbc_n_ind)
ggplot(data=final_table, aes(x=age, y=mean_bcc, col=Mating_System))+
  ggtitle("WBC Concentration vs. Age")+
  scale_color_brewer(palette="Spectral")+
  geom_smooth(method=lm, aes(y=regpredict), se= FALSE, fullrange=FALSE)+geom_point(aes(size=final_table$wbc_n_ind))+facet_wrap(~Scientific_Name, nrow=3)+
  xlim(0,40)


lymphocyte_age = read.csv("all species lymphocytes.csv")
lymphocyte_age1 = subset(lymphocyte_age, select = -c(birthdate, deathdate, age_at_death))
lymphocyte_age_pivot_longer1 = pivot_longer(data=lymphocyte_age1, cols=starts_with("y"), names_to="age", values_to = "bloodcellconcentration")
lymphocyte_age_pivot_longer = lymphocyte_age_pivot_longer1[-which(is.na(lymphocyte_age_pivot_longer1$bloodcellconcentration)),]
unique_species = unique(lymphocyte_age_pivot_longer$species)
lymphocyte_age_pivot_longer$age = sub("y", "", lymphocyte_age_pivot_longer$age)
lymphocyte_age_pivot_longer3 = merge(lymphocyte_age_pivot_longer, cumulative_data_finall[, c("species", "Mating_System")], by="species", all.x=FALSE)
lymphocyte_age_pivot_longer2 = merge(lymphocyte_age_pivot_longer3, cumulative_data_finall[,c("species", "Mating_Season_Duration")], by="species", all.x=FALSE)
lymphocyte_age_pivot_longer1 = merge(lymphocyte_age_pivot_longer2, cumulative_data_finall[,c("species", "lymphocytes_n_ind")], by="species", all.x=FALSE)
lymphocyte_age_pivot_longer = merge(lymphocyte_age_pivot_longer1, cumulative_data_finall[,c("species", "Scientific_Name")], by="species", all.x=FALSE)
count = ddply(lymphocyte_age_pivot_longer,.(age,species),nrow)
length(count$V1)
lymphocyte_age_pivot_longer = lymphocyte_age_pivot_longer[!duplicated(lymphocyte_age_pivot_longer), ]
lymphocyte_age_pivot_longer = lymphocyte_age_pivot_longer %>% group_by(species, age) %>% mutate(lymphocytes_n_ind = n())
final_table = lymphocyte_age_pivot_longer %>% group_by(age, species) %>%
  summarize(Mating_System, Scientific_Name, Mating_Season_Duration, lymphocytes_n_ind, mean_bcc = mean(bloodcellconcentration, na.rm=TRUE))
final_table = final_table[!duplicated(final_table), ]
final_table$age = strtoi(final_table$age)
final_table = final_table[order(final_table[,1] ),]
reg = lm(mean_bcc~age*Mating_System, data=final_table, weights=lymphocytes_n_ind)
summary(reg)
final_table$regpredict = predict.lm(reg, final_table, weights=lymphocytes_n_ind)
ggplot(data=final_table, aes(x=age, y=mean_bcc, col=Mating_System))+
  ggtitle("Lymphocyte Concentration vs. Age")+
  scale_color_brewer(palette="Spectral")+
  geom_smooth(method=lm, aes(y=regpredict), se= FALSE, fullrange=FALSE)+geom_point(aes(size=final_table$lymphocytes_n_ind))+facet_wrap(~Scientific_Name, nrow=3)+
  xlim(0,40)


# Figure 9 - Blood Cell Concentrations by Sex
clean_data= function(file_name){
  primate_data = read.csv(paste(file_name, ".csv", sep=""), header = T)
  
  
  
  # get mean for each subject
  for(i in 1:length(primate_data$subject)){
    primate_data$n_years[i] = length(which(!is.na(primate_data[i,12:82])))
    primate_data$mean_measurement[i]=mean(as.numeric(primate_data[i,12:82]), na.rm=T)
    
  }
  #primate_data_nona = primate_data[-which(is.na(primate_data$mean_measurement)),] # BODY WEIGHT
  primate_data_nona=primate_data
  
  primate_data_sp = data.frame(species=unique(primate_data_nona$species), sex="all", n_ind=NA, mean_measurement=NA)
  # creates function that tells how many data points there are for each species when species=x
  count_sp = function(x){length(primate_data_nona$subject[which(primate_data_nona$species==x)])}
  #gives the datapoint count for every species
  primate_data_sp$n_ind=as.numeric(lapply(paste(primate_data_sp$species), FUN=count_sp))
  #adds together all the datapoints
  
  #take the mean WBC % for every species as a whole
  mean_sp = function(x){mean(as.numeric(primate_data_nona$mean_measurement[which(primate_data_nona$species==x)]))}
  primate_data_sp$mean_measurement=as.numeric(lapply(paste(primate_data_sp$species), FUN=mean_sp))
  
  
  # get a mean for each sex within each species
  
  primate_data_sex = data.frame(species=rep(unique(primate_data_nona$species), 2),
                                sex=c(rep("Male", length(unique(primate_data_nona$species))),
                                      rep("Female", length(unique(primate_data_nona$species)))),
                                n_ind=NA, mean_measurement=NA)
  
  
  for (i in 1:length(primate_data_sex$species)){
    
    
    primate_data_sex$mean_measurement[i]=mean(primate_data_nona$mean_measurement[which(primate_data_nona$species==primate_data_sex$species[i] & primate_data_nona$sex==primate_data_sex$sex[i])])
    primate_data_sex$n_ind[i]=length(primate_data_nona$mean_measurement[which(primate_data_nona$species==primate_data_sex$species[i] & primate_data_nona$sex==primate_data_sex$sex[i])])
    
  }
  
  # final table with all data (male and female together and separate)- save
  combined_table = rbind(primate_data_sp, primate_data_sex)
  write.csv(combined_table, paste(file_name, "clean_data.csv"))
}

clean_data("monocytes all species")
clean_data("neutrophil seg all species")
clean_data("All Species Downloaded Data for Eosinophils")
clean_data("all species lymphocytes")
clean_data("basophils all species")
clean_data("wbc all species")
clean_data("neutrophil banded all species")
clean_data("body weight all species downloaded")

monocytes = read.csv("monocytes all species clean_data.csv", header=T)
neutrophils_seg = read.csv("neutrophil seg all species clean_data.csv", header=T)
eosinophils = read.csv("All Species Downloaded Data for Eosinophils clean_data.csv", header=T)
lymphocytes = read.csv("all species lymphocytes clean_data.csv", header=T)
basophils = read.csv("basophils all species clean_data.csv", header=T)
wbc = read.csv("wbc all species clean_data.csv", header=T)
neutrophils_band = read.csv("neutrophil banded all species clean_data.csv", header=T)
body_weight = read.csv("all species downloaded for body weight clean_data.csv", header=T)


monocytes = monocytes[,-1] #*******************
colnames(monocytes)[c(3, 4)] = c("monocytes_n_ind", "monocytes_mean_measurement")
colnames(monocytes)

neutrophils_seg = neutrophils_seg[,-1]
colnames(neutrophils_seg)[c(3, 4)] = c("neutrophils_n_ind", "neutrophils_mean_measurement")
colnames(neutrophils_seg)

eosinophils = eosinophils[,-1]
colnames(eosinophils)[c(3, 4)] = c("eosinophils_n_ind", "eosinophils_mean_measurement")

lymphocytes = lymphocytes[,-1]
colnames(lymphocytes)[c(3, 4)] = c("lymphocytes_n_ind", "lymphocytes_mean_measurement")

basophils = basophils[, -1]
colnames(basophils)[c(3, 4)] = c("basophils_n_ind", "basophils_mean_measurement")
colnames(basophils)

wbc = wbc[,-1]
colnames(wbc)[c(3, 4)] = c("wbc_n_ind", "wbc_mean_measurement")
colnames(wbc)

neutrophils_band = neutrophils_band[,-1]
colnames(neutrophils_band)[c(3, 4)] = c("neutrophils_n_ind", "neutrophils_mean_measurement")
colnames(neutrophils_band)

body_weight = body_weight[,-1]
colnames(body_weight)[c(3, 4)] = c("body_weight_n_ind", "body_weight_mean_measurement")
colnames(body_weight)

merge_dataset1 = merge(monocytes, neutrophils_seg, by=c("species", "sex"), all=TRUE)
merge_dataset1
merge_dataset2 = merge(merge_dataset1, eosinophils, by=c("species", "sex"), all=TRUE)
merge_dataset2
merge_dataset3 = merge(merge_dataset2, lymphocytes, by=c("species", "sex"), all=TRUE)
merge_dataset3
merge_dataset4 = merge(merge_dataset3, basophils, by=c("species", "sex"), all=TRUE)
merge_dataset4
merge_dataset5 = merge(merge_dataset4, wbc, by=c("species", "sex"), all=TRUE)
merge_dataset5
merge_dataset6 = merge(merge_dataset5, neutrophils_band, by=c("species", "sex"), all=TRUE)
merge_dataset6
merge_dataset7 = merge(merge_dataset6, body_weight, by=c("species", "sex"), all=TRUE)
merge_dataset7


duration_phy = read.csv("mating season duration for species phylogeny3.csv", header=T)
colnames(duration_phy)[1]="species"

str(duration_phy)


duration_phy = duration_phy[,c(1, 4, 6, 9, 10)] # 2 and move over
final_dataset = merge(merge_dataset7, duration_phy, by=c("species"), all=TRUE)
colnames(final_dataset)[c(5, 6, 15, 16)]= c("neutrophils_seg_n_ind", "neutrophils_seg_mean_measurement", "neutrophils_band_n_ind", "neutrophils_band_mean_measurement")
write.csv(final_dataset, file="finaldatasetwbcandweightmeanssexphylogeny.csv")

final_data = read.csv("finaldatasetwithinfo.csv")


final_data_lemur = final_data[final_data$Taxonomy == "strep",]
final_data_lemur
final_data_lemur_nosex = subset(final_data_lemur, sex=="all")

final_data_lemur_noall = subset(final_data_lemur, sex!="all")
#------------------------------------------------------------------------------------------
#install.packages("ggplot")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("tidyr")
library(tidyr)
#install.packages("ggpmisc")
library(ggpmisc)

colnames(final_data_lemur_noall)
final_data_lemur_noall$Scientific_Name2 = gsub(" ", "_", final_data_lemur_noall$Scientific_Name2)


pivot_longer_data = data.frame(final_data_lemur_noall[, c(2, 3, 5, 7, 9, 11, 13, 15, 17, 23)])
pivot_longer_data = merge(pivot_longer_data, cumulative_data_finall[, c("species", "Scientific_Name", "Mating_System")], by="species")

pivot_longer_official_data = pivot_longer(data=pivot_longer_data, cols = ends_with("_mean_measurement"), names_to="blood_cell_type", values_to= "mean_blood_cell")
gather(pivot_longer_official_data, key=sex, value="")


pivot_spread_data = pivot_longer_official_data %>%
  spread(key=sex, value=mean_blood_cell) %>%
  filter(Mating_System.y=="Monogamous"|Mating_System.y=="Polygynous"|Mating_System.y=="Polyandrous"|Mating_System.y=="Polygamous") %>%
  filter(blood_cell_type!="neutrophils_band_mean_measurement")

pivot_spread_data$blood_cell_type=factor(pivot_spread_data$blood_cell_type, labels=c("Basophils", "Eosinophils", "Lymphocytes", "Monocytes", "Neutrophils", "All WBC"))
#write.csv(pivot_spread_data, "male and female counts.csv")
ggplot(data=pivot_spread_data, aes(x=Male,
                                   y=Female))+
  ggtitle("Blood Cell Concentrations by Sex")+
  geom_text(aes(label=Scientific_Name), check_overlap = TRUE, size=2)+
  geom_abline(slope=1, intercept=0)+facet_wrap(~blood_cell_type, nrow=2, scales="free")+
  geom_point(aes(col=Mating_System.y)) + scale_color_brewer(palette="Spectral") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)





# Figure 10- sexual dimorphism phyogenetic signal- NOTHING CONCLUSIVE
# have not added to paper yet- should i?
pivot_spread_data
data_subset_lymph = subset(x=pivot_spread_data, blood_cell_type=="Lymphocytes")
data_subset_lymph = merge(data_subset_lymph, cumulative_data_finalll[,c(2,7)], by="Scientific_Name", all.x=FALSE)
data_subset_lymph["MaledivAverage_Ratio"] = (data_subset_lymph$Male)/(data_subset_lymph$Lymphocytes)
p = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')
p10 = facet_plot(p, panel="Lymphocyte Ratio", data=data_subset_lymph, geom=geom_point,
                aes(x=MaledivAverage_Ratio)) 
p10 + theme_tree2(panel.spacing = unit(5, "lines"))
BM = setNames(data_subset_lymph$MaledivAverage_Ratio,
              data_subset_lymph$Scientific_Name)
phylosig(tree, BM, method="lambda", test=TRUE)
phylosig(tree, BM, method="K")

data_subset_allwbc = subset(x=pivot_spread_data, blood_cell_type=="All WBC")
data_subset_allwbc = merge(data_subset_allwbc, cumulative_data_finalll[,c(2,5)], by = "Scientific_Name", all.x=FALSE)
data_subset_allwbc["MaledivAverage_Ratio"] = (data_subset_allwbc$Male)/(data_subset_allwbc$AllWBC)
p = ggtree(tree) + geom_tiplab(size=2)+ coord_cartesian(clip = 'off')
p10 = facet_plot(p, panel="WBC Ratio", data=data_subset_allwbc, geom=geom_point,
                 aes(x=MaledivAverage_Ratio)) 
p10 + theme_tree2(panel.spacing = unit(5, "lines"))
BM = setNames(data_subset_allwbc$MaledivAverage_Ratio,
              data_subset_allwbc$Scientific_Name)
phylosig(tree, BM, method="lambda", test=TRUE)
phylosig(tree, BM, method="K")


# Figure 11- the phylogenetic signal of the slope of wbc count over age
# not put into paper yet

# seeing if slope between WBC and age has a phylogenetic signal
FINDSLOPE = function(Scientific_Name, data1)
{
  subsetdm = data1[data1$Scientific_Name==Scientific_Name,]
  mod1 = lm(mean_age~age, data=subsetdm)
  sum1 = summary(mod1)
  b = sum1$coefficients[2,1]
  error = sum1$coefficients[2,2]
  slope = c(b, error)
  return(slope)
}

speciesnames = unique(wbc_age_lemurphy$Scientific_Name)
findata = matrix(NA, 23,3, byrow = TRUE)
# results = FINDSLOPE(Scientific_Name = "Eulemur_coronatus", data=wbc_age_lemurphy)


for (i in 1:23)
{
  data = as.data.frame(FINDSLOPE(Scientific_Name = speciesnames[i], data=wbc_age_lemurphy))
  findata[i,] = c(speciesnames[i], data1[1], data1[2])
  
}

findata1 = as.data.frame(findata)
colnames(findata1) = c("Scientific_Name", "ABSlope", "Error")
findata1$ABSlope = as.numeric(findata1$ABSlope)
findata1$Error = as.numeric(findata1$Error)

p5 = facet_plot(p, panel="ABSlope", data=findata1, geom=geom_point, #add weight by 1/error
                aes(x=findata1$ABSlope, size=1/sqrt(Error))) 

p5 + theme_tree2(panel.spacing = unit(5, "lines"))

BM = setNames(findata1$ABSlope,
              findata1$Scientific_Name)

phylosig(tree, BM, method="lambda")
phylosig(tree, BM, method="K")

1/sqrt(findata1$Error)


#age slope phy
# quarterly