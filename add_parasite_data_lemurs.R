library(dplyr)
library(tidyr)
library(stringr)
library(fossil)


#get big primate parasite dataset
gmpd_data = read.csv("~/Documents/Google Drive/FILES/A-G/Flowering Time in Primates/gmpd database query mysql wkbnch 12.1.2020.csv", header=T)

gmpd_data$HostCorrectedName_Corbet = gsub(" ","_",gmpd_data$HostCorrectedName_Corbet)
gmpd_data$HostCorrectedName_MSW05 = gsub(" ","_",gmpd_data$HostCorrectedName_MSW05)
#gmpd_data$HostReportedName = gsub(" ","_",gmpd_data$HostReportedName)

#traits data made in another sheet- analyses 5 Aug 2021
traits_data$Species[traits_data$Species %in% gmpd_data$HostCorrectedName_Corbet]
traits_data$Species[!traits_data$Species %in% gmpd_data$HostCorrectedName_MSW05]
#all species but 4 represented in the data! better than I thought! - missing Eulemur_flavifrons, Nycticebus_pygameus, Eulemur_sanfordi, Loris_tardigradus
#traits_data$Species[!traits_data$Species %in% gmpd_data$HostReportedName]

gmpd_data_lemursub = gmpd_data[gmpd_data$HostCorrectedName_MSW05 %in% traits_data$Species,]
summary(gmpd_data_lemursub)

#two studies that sampled for sexually transmitted diseases, prevalence was zero - 
#Irwin et al. 2010 (Propithecus diadema), Dutton et al. 2008 (Varecia rubra)

length(which(gmpd_data_lemursub$Direct==1)) #279 records
length(which(gmpd_data_lemursub$CloseT==1)) #119 records


##omit parasites labelled 'sp.', if there is already a named species 
#for the same genus recorded in the same host

#group by host and eliminate sp.s if a named species exists for the parasite genus
#also remove 0 prevalence (host has been found to NOT have that parasite)

byhostpar <- gmpd_data_lemursub %>%
  group_by(HostCorrectedName_MSW05)%>%   #Group by host
  filter(TotalPrevalence != 0 )   %>%   ## only take columns that don't have 0 prev, but retain NA
  mutate(parasitegenus=word(ParasiteCorrectedName, 1))%>% ## make a column for the name of the parasite genus
  mutate(parasitespecies=word(ParasiteCorrectedName,-1)) %>% #make a column for the parasite name
  ungroup() %>%
  group_by(HostCorrectedName_MSW05, parasitegenus)%>%   #Group by host AND parasite genus
  mutate(nonsp_in_gen = length(parasitegenus[parasitespecies != "sp."])>0) %>% #T/F column for whether there is another record of identified species in the same parasite genus
  ungroup() %>% 
  mutate(colsp = parasitespecies=="sp.") %>% #T/F for whether the record is a "sp."
  arrange(HostCorrectedName_MSW05, parasitegenus, parasitespecies) %>%
  filter(nonsp_in_gen == FALSE | colsp == FALSE) #throw out sp's for which there is another identified species in the genus (i.e. TRUE in both columns)

str(byhostpar) #5207 records

#write.csv(byhostpar, file = "checking sp issue.csv") # -- looks good!

##how many hosts and parasites?
length(unique(byhostpar$HostCorrectedName_MSW05)) #20 host species with at least one parasite with prev > 0
length(unique(byhostpar$ParasiteCorrectedName)) #105 parasite species

byhostpar1=byhostpar
byhostpar=as.data.frame(byhostpar[,c("ID","ParasiteCorrectedName","HostCorrectedName_MSW05","Citation")])

createMatrix = function(data) {
  data$Citation = factor(data$Citation)
  data$ParasiteCorrectedName = factor(data$ParasiteCorrectedName)
  m=matrix(0,nrow=length(levels(data$ParasiteCorrectedName)), 
           ncol=length(levels(data$Citation)), 
           dimnames=list(levels(data$ParasiteCorrectedName),levels(data$Citation)))
  
  for(i in 1:ncol(m)){
    for(j in 1:nrow(m)){
      data2 = data[which(data$Citation==colnames(m)[i]),]
      if(rownames(m)[j] %in% data2$ParasiteCorrectedName) m[j,i]=1
    }
  }
  m
}

combinations = split(byhostpar, list(byhostpar$HostCorrectedName_MSW05))
head(combinations)
length(combinations)
allmatrices= lapply(combinations, createMatrix)

#how many citations
num_cit = lapply(allmatrices, ncol)
length(num_cit)

#make citations back into a dataframe through some coaxing
num_cit2 = do.call("rbind", num_cit)
num_cit2 = as.matrix(num_cit2)
num_cit2 = as.data.frame(num_cit2)
num_cit2$Host=rownames(num_cit2)
rownames(num_cit2) = NULL
colnames(num_cit2) = c("Citations", "host")

#how many recorded parasites
PSR = lapply(allmatrices, nrow)
length(PSR)

#make true PSR estimates back into a dataframe through some coaxing
PSR_total = do.call("rbind", PSR)
PSR_total = as.matrix(PSR_total)
PSR_total = as.data.frame(PSR_total)
PSR_total$Host=rownames(PSR_total)
rownames(PSR_total) = NULL

colnames(PSR_total) = c("PSR", "host")
head(PSR_total)


#do chao2 estimate on all host matrices
chao_est_PSR = lapply(allmatrices, FUN=chao2)

#make chao estimates back into a dataframe through some coaxing
chao_est = do.call("rbind", chao_est_PSR)
chao_est = as.matrix(chao_est)
chao_est = as.data.frame(chao_est)
chao_est$Host=rownames(chao_est)
rownames(chao_est) = NULL

colnames(chao_est) = c("chao_est_PSR", "host")
head(chao_est)

psrtot_lemurs1=merge(PSR_total,num_cit2, by="host", all=T)
psrtot_lemurs=merge(psrtot_lemurs1,chao_est, by="host", all=T)

plot(chao_est_PSR~Citations, psrtot_lemurs)


# just look at direct, close transmitted diseases
byhostpar_direct = byhostpar1[byhostpar1$Direct==1,c("ID","ParasiteCorrectedName","HostCorrectedName_MSW05","Citation")]

combinations = split(byhostpar_direct, list(byhostpar_direct$HostCorrectedName_MSW05))
head(combinations)
length(combinations)
allmatrices= lapply(combinations, createMatrix)
PSR_direct = lapply(allmatrices, nrow)

#make true PSR estimates back into a dataframe through some coaxing
PSR_direct = do.call("rbind", PSR_direct)
PSR_direct = as.matrix(PSR_direct)
PSR_direct = as.data.frame(PSR_direct)
PSR_direct$Host=rownames(PSR_direct)
rownames(PSR_direct) = NULL

colnames(PSR_direct) = c("PSR_direct", "host")
head(PSR_direct)



# just look at direct, close transmitted diseases
byhostpar_close = byhostpar1[byhostpar1$CloseT==1,c("ID","ParasiteCorrectedName","HostCorrectedName_MSW05","Citation")]

combinations = split(byhostpar_close, list(byhostpar_close$HostCorrectedName_MSW05))
head(combinations)
length(combinations)
allmatrices= lapply(combinations, createMatrix)
PSR_close = lapply(allmatrices, nrow)

#make true PSR estimates back into a dataframe through some coaxing
PSR_close = do.call("rbind", PSR_close)
PSR_close = as.matrix(PSR_close)
PSR_close = as.data.frame(PSR_close)
PSR_close$Host=rownames(PSR_close)
rownames(PSR_close) = NULL

colnames(PSR_close) = c("PSR_close", "host")
head(PSR_close)

psrtot_lemurs1 = merge(psrtot_lemurs, PSR_direct, by="host",all=T)
psrtot_lemurs2 = merge(psrtot_lemurs1, PSR_close, by="host",all=T)
psrtot_lemurs2

psrtot_lemurs2$PSR_close[is.na(psrtot_lemurs2$PSR_close)] = 0

plot(PSR_close/PSR~Citations, psrtot_lemurs2)
plot(PSR_direct/PSR~Citations, psrtot_lemurs2)


psrtot_lemurs2$prop_close = psrtot_lemurs2$PSR_close/psrtot_lemurs2$PSR
psrtot_lemurs2$prop_direct = psrtot_lemurs2$PSR_direct/psrtot_lemurs2$PSR

write.csv(psrtot_lemurs2, "PSR data for lemurs.csv")

