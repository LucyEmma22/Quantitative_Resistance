library(gtools) 
library(stringr)
library(dplyr)
library(tidyr)
library(diptest)
library(multimode)
library(mixR)

################################################################################################################################################
# READ IN DATA
setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance/CSVs")
file_list<-list.files() # Creates a vector where each element is the name of a CSV in your folder e.g. [1] "Amikacin.csv" [2] "Amoxicillin-clavulanic acid (fixed).csv"
all_data<-data.frame() # Create an empty data frame called "all_data"
for (i in 1:length(file_list)){ # Set up a loop - let i be 1, then 2, then 3 etc and stop when i is equal to the length of the vector file_list
  data<-read.csv(file_list[i]) # Read in the CSV which is number i in the file_list
  data<-data[1:24] # Select the columns with data in, removing any extra columns containing only NAs. Yours will be [1:24] as you've added the extra columns
  data$CI<-ifelse(is.na(data$CI)==TRUE,"",data$CI)
  data<-na.omit(data) # Remove any rows at the bottom that contain only NAs
  data$Antibiotic<-file_list[i] # Add a column called "Antibiotic" where each row is the name of the antibiotic for this CSV (element i in file_list)
  data$Antibiotic<-gsub(".csv","", data$Antibiotic) # Remove ".csv" from the end of all the antibiotic names
  names(data)[1] <- "Bacteria" # Name the first column "Bacteria"
  all_data<-smartbind(all_data,data) # Add the data from the current CSV onto a dataframe containing the data from all the CSVs
}

################################################################################################################################################
# TIDY DATAFRAME

setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance")
data_ECOFF<-all_data %>%
  filter(Distributions>=5) %>% # Only include combinations made from 5 or more independent MIC distributions
  dplyr::rename(ECOFF=X.T.ECOFF) %>% 
  filter(is.na(as.numeric(ECOFF))==F) %>%
  gather("MIC","freq",2:20) %>% 
  full_join(read.csv("Drug Classes.csv"),by="Antibiotic") %>% 
  separate(Bacteria,c("genus","species"),sep=" ",extra="merge",remove=FALSE) %>% 
  dplyr::select(-CI) %>% 
  mutate(ECOFF=as.numeric(ECOFF)) %>%
  mutate(MIC=as.numeric(sub("X", "", MIC))) %>% 
  mutate(genus=sub(",", "", genus)) %>%
  filter(is.na(Bacteria)==FALSE) %>% 
  filter(ATC_code!="Antifungal" & ATC_code!="Antiprotozoan") %>% # Remove fungi and protozoa
  #filter((grepl('J01', ATC_code))) %>% # Include only antibiotics for systemic use
  mutate(Bacteria_Antibiotic=paste0(Bacteria,"_",Antibiotic))

data_ECOFF<-data_ECOFF %>% full_join(data_ECOFF %>% filter(MIC>ECOFF) %>% group_by(Bacteria_Antibiotic) %>% summarise_at(vars(freq),list(resistant = sum)))
bins<-data.frame(min=c(0.001,sort(unique(data_ECOFF$MIC))[-length(sort(unique(data_ECOFF$MIC)))]),MIC=sort(unique(data_ECOFF$MIC)))
data_ECOFF<-full_join(data_ECOFF,bins,by="MIC") %>% filter(resistant!=0)
write.csv(data_ECOFF,file="data_ECOFF.csv",row.names = FALSE)

################################################################################################################################################
# MIXTURE MODELS + TESTS

bacteria_antibiotic_list<-unique(data_ECOFF$Bacteria_Antibiotic)
bimodality_results<-data.frame()

for (i in 1:length(bacteria_antibiotic_list)){
  
  data<-filter(data_ECOFF,Bacteria_Antibiotic==bacteria_antibiotic_list[i])
  X<-dplyr::select(data,min,MIC,freq) %>% mutate(min=log2(min)) %>% mutate(MIC=log2(MIC)) 
  X_reinstate<-reinstate(as.matrix(X))
  
  ############### MIXTURE MODELS ###############
  
  # Bimodal Model
  tryCatch({
    b<-mixfit(X_reinstate,ncomp=2,ev=FALSE)
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Bimodal not calculated"), "\n")})
  
  # Unimodal Model
  tryCatch({
    u<-mixfit(X_reinstate,ncomp=1)
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Unimodal not calculated"), "\n")})
  
  # Are there 2 peaks in the bimodal curve?
  if (is.null(b)==F){
    df2<-data.frame(x=density(b)$x,y=density(b)$y) %>%  mutate(Diff = y - lag(y))
    df2$Diff[1]<-0
    for (j in 1:(nrow(df2)-1)){if(df2$Diff[j]<0 & df2$Diff[j+1]>0){peaks<-"2 Peaks"
    break}else{peaks<-"1 Peak"}}
  }else{peaks<-NA}
  
  # Overlap
  prop_overlap<-ifelse(is.null(b)==F,max( sum(pmin(density(b)$comp[,1],density(b)$comp[,2]))/sum(density(b)$comp[,1]) , sum(pmin(density(b)$comp[,1],density(b)$comp[,2]))/sum(density(b)$comp[,2]) ),NA)
  
  # Probability Best
  uBIC<-ifelse(is.null(u)==F,u$bic,NA)
  bBIC<-ifelse(is.null(b)==F,b$bic,NA)
  ri<-exp((min(uBIC,bBIC)-max(uBIC,bBIC))/2)
  prob_best<-ri/(ri+1)
  bimodality_results$prob_bimodal<-ifelse(bimodality_results$best_model=="Bimodal",1-bimodality_results$prob_best,bimodality_results$prob_best)
  
  ############### BIMODALITY TESTS ###############
  
  # Dip Test
  diptest<-NA
  tryCatch({
    diptest<-dip.test(X_reinstate)
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Dip Test not calculated"), "\n")})
  
  # Critical Bandwidth
  multimode_criticalbandwidth_si<-NA
  tryCatch({
    multimode_criticalbandwidth_si<-modetest(X_reinstate,mod0=1,method="SI")
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Critical Bandwidth not calculated"), "\n")})
  
  # Excess Mass
  multimode_excessmass_acr<-NA
  tryCatch({
    multimode_excessmass_acr<-modetest(X_reinstate,mod0=1,method="ACR")
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Excess Mass not calculated"), "\n")})
  
  ############### RESULTS ###############
  
  bimodality_results<-rbind(bimodality_results,data.frame(Bacteria_Antibiotic=bacteria_antibiotic_list[i],
                                                          ECOFF,
                                                          distributions=unique(data$Distributions),
                                                          observations=unique(data$Observations),
                                                          resistance_freq=sum(filter(X,MIC>log2(ECOFF))$freq)/sum(X$freq),
                                                          best_model=ifelse(uBIC<bBIC,"Unimodal","Bimodal"),
                                                          unimodalBIC=ifelse(is.null(u)==F,u$bic,NA),
                                                          bimodalBIC=ifelse(is.null(b)==F,b$bic,NA),
                                                          prob_best,
                                                          peaks,
                                                          prop_overlap,
                                                          diptest=as.numeric(diptest[1]),diptest_pval=as.numeric(diptest[2]),
                                                          criticalbandwidth=as.numeric(multimode_criticalbandwidth_si[1]),criticalbandwidth_pval=as.numeric(multimode_criticalbandwidth_si[2]),
                                                          excessmass=as.numeric(multimode_excessmass_acr[1]),excessmass_pval=as.numeric(multimode_excessmass_acr[2])))
  print(i) # keep track of how the loop is running
} 

bimodality_results$prob_bimodal<-ifelse(bimodality_results$best_model=="Bimodal",1-bimodality_results$prob_best,bimodality_results$prob_best)
bimodality_results$prop_overlap<-ifelse(bimodality_results$best_model=="Unimodal",1,bimodality_results$prop_overlap)
write.csv(bimodality_results,file="bimodality_results.csv",row.names = FALSE)