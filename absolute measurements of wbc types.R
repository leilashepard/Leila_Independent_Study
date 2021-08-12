library(tidyr)
setwd("./PAD data")

#downloaded the WBC data from PAD on 5 Aug 2021

basophils = read.csv("Basophils-Abs_subjects_measurements_yearly_50th_percentile.csv", header=T)
eosinophils = read.csv("Eosinophils-Abs_subjects_measurements_yearly_50th_percentile.csv", header=T)
lymphocytes = read.csv("Lymphocytes-Abs_subjects_measurements_yearly_50th_percentile.csv", header=T)
monocytes = read.csv("Monocytes-Abs_subjects_measurements_yearly_50th_percentile.csv", header=T)
neutrophils = read.csv("Neutrophil_Seg-Abs_subjects_measurements_yearly_50th_percentile.csv", header=T) #for banded, most measurements are 0
rbcs = read.csv("RBC_subjects_measurements_yearly_50th_percentile.csv", header=T)
wbcs = read.csv("WBC_subjects_measurements_yearly_50th_percentile.csv", header=T)
bodymass = read.csv("Body Mass_subjects_measurements_yearly_50th_percentile.csv", header=T)

#make one big dataframe with all the wbc types
leuks = rbind(basophils,eosinophils,lymphocytes,monocytes,neutrophils,rbcs,wbcs,bodymass)
length(unique(leuks$subject)) #1411
#separate by age
leuks_longer = leuks %>% pivot_longer(cols=y0:y70, names_to = "age", values_to = "value", values_drop_na = "TRUE")

subyr = leuks_longer %>% group_by(subject) %>% summarize(nyears=length(unique(age)))
length(subyr$subject)
View(subyr)
sum(subyr$nyears)

#take mean for each sex-age
leuks_sum_sex_age = leuks_longer %>% group_by(species, sex, age, measurement) %>% summarize(mean_value = mean(value), n_ind=length(unique(subject)))
leuks_sum_sex_agew = leuks_sum_sex_age %>% pivot_wider(id_cols=c(species, sex, age, measurement), names_from=measurement, values_from = c(mean_value, n_ind))
write.csv(leuks_sum_sexw, "PAD_leukocytes_strepsirrhines_sepsexage.csv")


#take mean for each age from male and female
leuks_sum_age = leuks_sum_sex_age %>% group_by(species, age, measurement) %>% summarize(mean_value = mean(mean_value), n_ind=sum(n_ind))
leuks_sum_agew = leuks_sum_age %>% pivot_wider(id_cols=c(species, age, measurement), names_from=measurement, values_from = c(mean_value, n_ind))
write.csv(leuks_sum_agew, "PAD_leukocytes_strepsirrhines_sepage.csv")

#take mean for each sex across the ages
leuks_sum = leuks_sum_sex_age %>% group_by(species, sex, measurement) %>% summarize(mean_value = mean(mean_value))
leuks_nind = leuks_longer %>% group_by(species, sex, measurement) %>% summarize(n_ind = length(unique(subject)))

#separate estimates for males and females
leuks_sum_sex = leuks_sum %>% left_join(leuks_nind)
leuks_sum_sexw = leuks_sum_sex %>% pivot_wider(id_cols=c(species, sex, measurement), names_from=measurement, values_from = c(mean_value, n_ind))
write.csv(leuks_sum_sexw, "PAD_leukocytes_strepsirrhines_sepsex.csv")

#mean of males and females for each species
leuks_sum_all = leuks_sum_sex %>% group_by(species, measurement) %>% summarize(mean_value=mean(mean_value), n_ind=sum(n_ind))
leuks_sum_allw = leuks_sum_all %>% pivot_wider(id_cols=c(species,measurement), names_from=measurement, values_from = c(mean_value, n_ind))

write.csv(leuks_sum_allw, "PAD_leukocytes_strepsirrhines_all.csv")

  
