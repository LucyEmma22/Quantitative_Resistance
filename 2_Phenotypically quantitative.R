################################################################################################################################################
# IMPORT CSV FILES
setwd("~/OneDrive - University of Edinburgh/Bhavya/CSVs")
library(gtools) # Library for 'smartbind' function
library(stringr) # Library for 'gsub' function
################################################################################################################################################

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
setwd("~/OneDrive - University of Edinburgh/Bhavya")
library(dplyr)
library(tidyr)
################################################################################################################################################

data_ECOFF<-all_data %>%
  filter(Distributions>=5) %>% 
  dplyr::rename(ECOFF=X.T.ECOFF) %>% 
  filter(is.na(as.numeric(ECOFF))==F) %>%
  gather("MIC","freq",2:20) %>% 
  full_join(read.csv("Drug Classes Bhavya.csv"),by="Antibiotic") %>% 
  separate(Bacteria,c("genus","species"),sep=" ",extra="merge",remove=FALSE) %>% 
  dplyr::select(-CI) %>% 
  mutate(ECOFF=as.numeric(ECOFF)) %>%
  mutate(MIC=as.numeric(sub("X", "", MIC))) %>% 
  mutate(genus=sub(",", "", genus)) %>%
  filter(is.na(Bacteria)==FALSE) %>% 
  filter(ATC_code!="Antifungal" & ATC_code!="Antiprotozoan") %>% 
  #filter(data_ECOFF,(grepl('J01', ATC_code))) %>% 
  mutate(Bacteria_Antibiotic=paste0(Bacteria,"_",Antibiotic))

data_ECOFF<-data_ECOFF %>% full_join(data_ECOFF %>% filter(MIC>ECOFF) %>% group_by(Bacteria_Antibiotic) %>% summarise_at(vars(freq),list(resistant = sum)))
bins<-data.frame(min=c(0.001,sort(unique(data_ECOFF$MIC))[-length(sort(unique(data_ECOFF$MIC)))]),MIC=sort(unique(data_ECOFF$MIC)))
data_ECOFF<-full_join(data_ECOFF,bins,by="MIC") 

################################################################################################################################################
# MIXTURE MODELS + TESTS
library(diptest)
library(Rfolding)
library(multimode)
library(mixR)
library(gridExtra)
library(ggplot2)
################################################################################################################################################
data_ECOFF<-data_ECOFF %>% filter(resistant!=0)
unimodal_plots<-list()
bimodal_plots<-list()
bacteria_antibiotic_list<-unique(data_ECOFF$Bacteria_Antibiotic)
bimodality_results<-data.frame()
for (i in 1:length(bacteria_antibiotic_list)){
  
  data<-filter(data_ECOFF,Bacteria_Antibiotic==bacteria_antibiotic_list[i])
  X<-dplyr::select(data,min,MIC,freq) %>% mutate(min=log2(min)) %>% mutate(MIC=log2(MIC)) 
  X_reinstate<-reinstate(as.matrix(X))
  
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
  
  
  ################################################################################################################################################
  
  # Dip Test
  diptest<-NA
  tryCatch({
    diptest<-dip.test(X_reinstate)
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Dip Test not calculated"), "\n")})
  
  # Folding Test
  foldingtest<-NA
  tryCatch({
    foldingtest<-folding.test(as.matrix(X_reinstate))
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Folding Test not calculated"), "\n")})
  
  # Multimode
  multimode_criticalbandwidth_si<-NA
  tryCatch({
    multimode_criticalbandwidth_si<-modetest(X_reinstate,mod0=1,method="SI")
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Critical Bandwidth not calculated"), "\n")})
  
  multimode_excessmass_acr<-NA
  tryCatch({
    multimode_excessmass_acr<-modetest(X_reinstate,mod0=1,method="ACR")
  }, error=function(e){cat(paste0("ERROR : ",bacteria_antibiotic_list[i], " Excess Mass not calculated"), "\n")})
  
  ################################################################################################################################################
  # Plots
  
  if (is.null(b)==F){
    
    ECOFF<-unique(data$ECOFF)
    X<-filter(X,freq!=0)
    X_reinstate_df<-data.frame(X_reinstate)
    X_reinstate_df$colour<-ifelse(X_reinstate_df$X_reinstate<=log2(ECOFF),"dodgerblue","firebrick")
    density_bimodal<-data.frame(x=density(b)$x,y_all=density(b)$y,y_s=density(b)$comp[,1],y_r=density(b)$comp[,2])
    density_unimodal<-data.frame(x=density(u)$x,y=density(u)$y)
    
    bimodal_plots[[i]]<-ggplot(data=X_reinstate_df,aes(x=X_reinstate)) +
      geom_histogram(breaks=unique(c(X$min,X$MIC)),aes(fill=colour,y=stat(count)/sum(stat(count))),colour="white",size=0.1) + 
      scale_fill_manual(values=c("dodgerblue","firebrick"))+
      geom_polygon(data=density_bimodal,aes(x=x,y=y_s),fill="dodgerblue",alpha=0.5) + 
      geom_polygon(data=density_bimodal,aes(x=x,y=y_r),fill="firebrick",alpha=0.5) + 
      geom_line(data=density_bimodal,aes(x=x,y=y_s),colour="dodgerblue3",size=0.5) + 
      geom_line(data=density_bimodal,aes(x=x,y=y_r),colour="firebrick3",size=0.5) + 
      geom_line(data=density_bimodal,aes(x=x,y=y_all),linetype="dashed",size=0.5) + 
      theme_void() + 
      #labs(title=paste0(unique(data$Bacteria),"\n",unique(data$Antibiotic)),x=NULL, y=NULL) +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(legend.position="none")
    
    unimodal_plots[[i]]<-ggplot(data=X_reinstate_df,aes(x=X_reinstate)) +
      geom_histogram(breaks=unique(c(X$min,X$MIC)),aes(fill=colour,y=stat(count)/sum(stat(count))),colour="white",size=0.1) + 
      scale_fill_manual(values=c("dodgerblue","firebrick"))+
      geom_polygon(data=density_unimodal,aes(x=x,y=y),fill="dodgerblue",alpha=0.2) + 
      geom_polygon(data=density_unimodal,aes(x=x,y=y),fill="firebrick",alpha=0.2) + 
      geom_line(data=density_unimodal,aes(x=x,y=y),linetype="dashed",size=0.5) + 
      theme_void() + 
      #labs(title=paste0(unique(data$Bacteria),"\n",unique(data$Antibiotic)),x=NULL, y=NULL) +
      theme(plot.title = element_text(hjust = 0.5))+
      theme(legend.position="none")
  }
  ################################################################################################################################################
  
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
                                                          foldingtest=as.numeric(foldingtest[3]),foldingtest_pval=as.numeric(foldingtest[2]),
                                                          criticalbandwidth=as.numeric(multimode_criticalbandwidth_si[1]),criticalbandwidth_pval=as.numeric(multimode_criticalbandwidth_si[2]),
                                                          excessmass=as.numeric(multimode_excessmass_acr[1]),excessmass_pval=as.numeric(multimode_excessmass_acr[2])))
  
  print(i)
} 
bimodality_results$prob_bimodal<-ifelse(bimodality_results$best_model=="Bimodal",1-bimodality_results$prob_best,bimodality_results$prob_best)

################################################################################################################################################

#rownames(bimodality_results)<-NULL
bimodality_results$number<-1:nrow(bimodality_results)

bimodality_results$foldingtest_prediction<-ifelse(bimodality_results$foldingtest<1,2,1)
bimodality_results$diptest_prediction<-ifelse(bimodality_results$diptest_pval<0.05,2,1)
bimodality_results$cb_prediction<-ifelse(bimodality_results$criticalbandwidth_pval<0.05,2,1)
bimodality_results$em_prediction<-ifelse(bimodality_results$excessmass_pval<0.05,2,1)
bimodality_results$test_prediction<-(bimodality_results$diptest_prediction + bimodality_results$cb_prediction + bimodality_results$em_prediction + bimodality_results$foldingtest_prediction)/4
bimodality_results$model_prediction<-ifelse(bimodality_results$best_model=="Unimodal"|bimodality_results$peaks=="1 Peak",1,2)

bimodality_results$overlap<-ifelse(bimodality_results$best_model=="Unimodal",1,bimodality_results$prop_overlap)

bimodality_results<-bimodality_results[, c(20, 1:19,21:27)]

################################################################################################################################################
# PLOTS
################################################################################################################################################

# Dot and bar plots for tests and models
ggplot(data.frame(table(bimodality_results$model_prediction,bimodality_results$test_prediction)),aes(Var1,Var2))+
  geom_point(aes(size=Freq,colour=Freq))+
  scale_size_continuous(range = c(-1, 30))+
  theme_classic()+
  labs(x="Model",y="Tests")+
  theme(legend.position="none")+
  scale_colour_gradient(low = "white", high = "mediumpurple")+
  geom_text(aes(label=Freq))

ggplot(na.omit(bimodality_results),aes(x=as.factor(test_prediction)))+
  geom_bar(fill="mediumpurple")+
  facet_wrap(.~model_prediction)+
  theme_light()+
  labs(x="Test Prediction",y="Count")

# Distribution of Sample Sizes
sample_sizes<-distinct(dplyr::select(data_ECOFF,Observations,Bacteria_Antibiotic))
ggplot(sample_sizes,aes(log(Observations)))+
  geom_histogram(bins=20,colour="black",boundary=0,fill="lightgrey")+
  theme_light()+
  labs(title="Sample Sizes in EUCAST Data",x="Log (Sample Size)",y="Frequency")+
  theme(plot.title = element_text(hjust = 0.5))


# Sample Size Vs Number Resistant
ggplot(data=bimodality_results,aes(log10(observations),log10(observations*resistance_freq)))+
  geom_point(colour="black")+
  theme_light()+
  labs(x="Log10 (Sample Size)", y="Log10 (Number Resistant)")

# Number of Tests/Models that weren't run
data.frame(Test=c("Unimodal Model","Bimodal Model","Dip Test","Folding Test","Critical Bandwidth","Excess Mass"),
           Missing=c(nrow(subset(bimodality_results,is.na(bimodality_results$unimodalBIC))),
                     nrow(subset(bimodality_results,is.na(bimodality_results$bimodalBIC))),
                     nrow(subset(bimodality_results,is.na(bimodality_results$diptest))),
                     nrow(subset(bimodality_results,is.na(bimodality_results$foldingtest))),
                     nrow(subset(bimodality_results,is.na(bimodality_results$criticalbandwidth))),
                     nrow(subset(bimodality_results,is.na(bimodality_results$excessmass)))
           ))

################################################################################################################################################
# TEST P-VALUE PLOTS

# Calculations
test_results<-bimodality_results %>% dplyr::select(criticalbandwidth_pval,excessmass_pval,diptest_pval) %>% gather("test","pval",1:3) %>% na.omit()

em_mean<-mean(filter(test_results,test=="excessmass_pval")$pval)
dt_mean<-mean(filter(test_results,test=="diptest_pval")$pval)
cb_mean<-mean(filter(test_results,test=="criticalbandwidth_pval")$pval)

em_bimodal<-nrow(filter(test_results,test=="excessmass_pval" & pval<0.05))/nrow(filter(test_results,test=="excessmass_pval"))
dt_bimodal<-nrow(filter(test_results,test=="diptest_pval" & pval<0.05))/nrow(filter(test_results,test=="diptest_pval"))
cb_bimodal<-nrow(filter(test_results,test=="criticalbandwidth_pval" & pval<0.05))/nrow(filter(test_results,test=="criticalbandwidth_pval"))
print(paste0("Excess Mass:",em_bimodal," ___ Dip Test:",dt_bimodal," ___ Critical Bandwidth: ",cb_bimodal))

# Violin Plots - Pval
all_means<-data.frame(mean=c(cb_mean,dt_mean,em_mean),test=c("criticalbandwidth_pval","diptest_pval","excessmass_pval"),y=c(1,2,3),bimodal=c(cb_bimodal,dt_bimodal,em_bimodal))
pval_violin_plot<-ggplot(test_results,aes(pval,test))+
  geom_jitter(aes(colour=test),size=0.3,alpha=0.3,height=0.45)+
  geom_violin(aes(fill=test),colour=NA,alpha=0.6)+
  scale_fill_manual(values=c("mediumpurple","mediumseagreen","goldenrod"))+
  scale_colour_manual(values=c("mediumpurple","mediumseagreen","goldenrod"))+
  theme_light()+
  labs(x="P-Value",y=NULL)+
  theme(legend.position="none")+
  geom_segment(aes(x = mean, y = y-0.45, xend = mean, yend = y+0.45,colour=test),data=all_means)+
  geom_vline(xintercept=0.05,linetype="dashed")+
  #geom_text(data=all_means,aes(x=0.07,y=y,label=round(unimodal,2)),size=3,hjust=0)+
  #geom_text(data=all_means,aes(x=0.03,y=y,label=round(bimodal,2)),size=3,hjust=1)+
  scale_y_discrete(labels=c("excessmass_pval" = str_wrap("Excess Mass",5), "diptest_pval" = "Dip Test","criticalbandwidth_pval" = str_wrap("Critical Bandwidth",5)))

prop_bimodal_tests<-data.frame(table((filter(test_results,pval<0.05))$test)) %>% full_join(data.frame(table(test_results$test)),by="Var1") %>% mutate(prop_bimodal=Freq.x/Freq.y)
test_bar_plot<-ggplot(prop_bimodal_tests,aes(x=prop_bimodal,y=Var1))+
  geom_bar(stat="identity",aes(fill=Var1),alpha=0.8,colour="black")+
  theme_light()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("mediumpurple","mediumseagreen","goldenrod"))+
  labs(x="Proportion Bimodal (P-Value < 0.05)",y=NULL)+
  geom_label(aes(label=round(prop_bimodal,2)))+
  scale_y_discrete(labels=c("excessmass_pval" = str_wrap("Excess Mass",5), "diptest_pval" = "Dip Test","criticalbandwidth_pval" = str_wrap("Critical Bandwidth",5)))

grid.arrange(pval_violin_plot,test_bar_plot)
################################################################################################################################################
# EXAMPLE PLOTS

library(ggpubr)

label<-c("I: Low Overlap (10^-13%)", "I: High Overlap (10%)","II: 1 Peak","III: 1 Component")
title<-c("S.agalactiae + CIP", "K.pneumoniae + CST", "E.coli + C/T", "C.jejuni + SXT")
title_col<-c("mediumpurple","mediumseagreen","goldenrod","mediumvioletred")
examples<-c(
  "Streptococcus agalactiae_Ciprofloxacin",
  "Klebsiella pneumoniae_Colistin",
  "Escherichia coli_Ceftolozane-tazobactam",
  "Campylobacter jejuni_Trimethoprim-sulfamethoxazole "
  )

example_results<-data.frame()
example_plots<-list()
for (i in 1:length(examples)){
  data<-filter(data_ECOFF,Bacteria_Antibiotic==examples[i])
  X<-dplyr::select(data,min,MIC,freq) %>% mutate(min=log2(min)) %>% mutate(MIC=log2(MIC)) 
  X_reinstate<-reinstate(as.matrix(X))
  b<-mixfit(X_reinstate,ncomp=2,ev=FALSE)
  u<-mixfit(X_reinstate,ncomp=1)
  # Are there 2 peaks in the bimodal curve?
  df2<-data.frame(x=density(b)$x,y=density(b)$y) %>%  mutate(Diff = y - lag(y))
  df2$Diff[1]<-0
  for (j in 1:(nrow(df2)-1)){if(df2$Diff[j]<0 & df2$Diff[j+1]>0){peaks<-"2 Peaks"
  break}else{peaks<-"1 Peak"}}
  # Overlap
  prop_overlap<-ifelse(is.null(b)==F,max( sum(pmin(density(b)$comp[,1],density(b)$comp[,2]))/sum(density(b)$comp[,1]) , sum(pmin(density(b)$comp[,1],density(b)$comp[,2]))/sum(density(b)$comp[,2]) ),NA)
  # Probability Best
  uBIC<-ifelse(is.null(u)==F,u$bic,NA)
  bBIC<-ifelse(is.null(b)==F,b$bic,NA)
  ri<-exp((min(uBIC,bBIC)-max(uBIC,bBIC))/2)
  prob_best<-ri/(ri+1)
  
  ECOFF<-unique(data$ECOFF)
  X<-filter(X,freq!=0)
  X_reinstate_df<-data.frame(X_reinstate)
  X_reinstate_df$colour<-ifelse(X_reinstate_df$X_reinstate<=log2(ECOFF),"dodgerblue","firebrick")
  if(b$bic<u$bic){df_density<-data.frame(x=density(b)$x,y_all=density(b)$y,y_s=density(b)$comp[,1],y_r=density(b)$comp[,2])
  }else{df_density<-data.frame(x=density(u)$x,y_all=density(u)$y,y_s=density(u)$y,y_r=density(u)$y)}
  example_plots[[i]]<-ggplot(data=X_reinstate_df,aes(x=X_reinstate)) +
    geom_histogram(breaks=unique(c(X$min,X$MIC)),aes(fill=colour,y=stat(count)/sum(stat(count))),colour="white",size=0.1) + 
    scale_fill_manual(values=c("dodgerblue","firebrick"))+
    geom_polygon(data=df_density,aes(x=x,y=y_s),fill="dodgerblue",alpha=0.3) + 
    geom_polygon(data=df_density,aes(x=x,y=y_r),fill="firebrick",alpha=0.3) + 
    geom_line(data=df_density,aes(x=x,y=y_s),colour="dodgerblue3",size=0.5) + 
    geom_line(data=df_density,aes(x=x,y=y_r),colour="firebrick3",size=0.5) + 
    geom_line(data=df_density,aes(x=x,y=y_all),linetype="dashed",size=0.5) + 
    theme_void() + 
    labs(title=label[i],subtitle=title[i],x=NULL, y=NULL) +
    theme(plot.title = element_text(hjust = 0.5,size=10,colour = title_col[i],face="bold"),plot.subtitle = element_text(hjust = 0.5,size=8),legend.position="none")
  
  example_results<-rbind(example_results,data.frame(Bacteria_Antibiotic=bacteria_antibiotic_list[i],
                                                    ECOFF,
                                                    distributions=unique(data$Distributions),
                                                    observations=unique(data$Observations),
                                                    resistance_freq=sum(filter(X,MIC>log2(ECOFF))$freq)/sum(X$freq),
                                                    best_model=ifelse(uBIC<bBIC,"Unimodal","Bimodal"),
                                                    unimodalBIC=ifelse(is.null(u)==F,u$bic,NA),
                                                    bimodalBIC=ifelse(is.null(b)==F,b$bic,NA),
                                                    prob_best,
                                                    peaks,
                                                    prop_overlap))
}

all_examples_plots<-grid.arrange(grobs=example_plots,nrow=1)

################################################################################################################################################
# BAR PLOT

plot_data<-data.frame(Group=c("I","II","III"),Count=c(nrow(bimodality_results %>% filter(best_model=="Bimodal", peaks=="2 Peaks")),
                                                   nrow(bimodality_results %>% filter(best_model=="Bimodal", peaks=="1 Peak")),
                                                   nrow(bimodality_results %>% filter(best_model=="Unimodal"))))

library(ggpattern)
bar_plot<-ggplot(plot_data,aes(x=Count,y=Group))+
  geom_bar_pattern(stat="identity",alpha=0.8,colour="black",aes(fill=Group,pattern_fill = Group, pattern_colour=Group), pattern = 'stripe',pattern_density=0.4,pattern_size=0,pattern_spacing=1,pattern_angle=20,pattern_alpha=0.8)+
  scale_pattern_fill_manual(values = c("mediumpurple",NA,NA)) +
  scale_pattern_colour_manual(values = c("mediumpurple",NA,NA)) +
  scale_fill_manual(values=c("mediumseagreen","goldenrod","mediumvioletred"))+
  theme_light()+
  theme(legend.position="none")+
  geom_label(aes(label=Count))

################################################################################################################################################
# OVERLAP PLOT

# Violin Plot - all
bimodality_results$col<-ifelse(bimodality_results$Bacteria_Antibiotic==examples[1],"mediumpurple",
                               ifelse(bimodality_results$Bacteria_Antibiotic==examples[2],"mediumseagreen",
                               ifelse(bimodality_results$Bacteria_Antibiotic==examples[3],"goldenrod",
                               ifelse(bimodality_results$Bacteria_Antibiotic==examples[4],"mediumvioletred","black"))))
overlap_mean_all<-mean(bimodality_results$overlap,na.rm=T)
overlap_plot_all<-ggplot(bimodality_results,aes(overlap,1))+
  #geom_jitter(data=filter(bimodality_results,!Bacteria_Antibiotic %in% examples),aes(x=overlap,y=1),size=1,shape=16,height=0.4,colour="mediumpurple",stroke=0)+
  geom_jitter(size=1,shape=16,height=0.4,colour="black",alpha=0.3)+
  geom_violin(alpha=0.6,fill="grey",color=NA)+
  #geom_boxplot(width=0.1,fill="darkgrey",alpha=0.7,colour="black",outlier.size = -1)+
  #geom_vline(xintercept=overlap_mean_all,colour="darkgrey")+
  geom_segment(x = overlap_mean_all, y = 0.6, xend = overlap_mean_all, yend = 1.4,colour="grey")+
  geom_point(data=filter(bimodality_results,Bacteria_Antibiotic %in% examples),aes(x=overlap,y=1),colour=filter(bimodality_results,Bacteria_Antibiotic %in% examples)$col,size=4,shape=4,stroke=1.5)+
  theme_light()+
  labs(x="Overlap",y="Frequency")+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

grid.arrange(bar_plot,overlap_plot_all,all_examples_plots,nrow=3,heights=c(2/5,2/5,1/5))


pval_violin_plot<-ggplot(test_results,aes(pval,test))+
  geom_jitter(aes(colour=test),size=0.3,alpha=0.3,height=0.45)+
  geom_violin(aes(fill=test),colour=NA,alpha=0.6)+
  scale_fill_manual(values=c("goldenrod","mediumseagreen","mediumpurple"))+
  scale_colour_manual(values=c("goldenrod","mediumseagreen","mediumpurple"))+
  theme_light()+
  labs(x="P-Value",y=NULL)+
  theme(legend.position="none")+
  geom_segment(aes(x = mean, y = y-0.45, xend = mean, yend = y+0.45,colour=test),data=all_means)+
  geom_vline(xintercept=0.05,linetype="dashed")+
  #geom_text(data=all_means,aes(x=0.07,y=y,label=round(unimodal,2)),size=3,hjust=0)+
  #geom_text(data=all_means,aes(x=0.03,y=y,label=round(bimodal,2)),size=3,hjust=1)+
  scale_y_discrete(labels=c("excessmass_pval" = str_wrap("Excess Mass",5), "diptest_pval" = "Dip Test","criticalbandwidth_pval" = str_wrap("Critical Bandwidth",5)))


# Violin Plot - 2 Peak
overlap_mean_twopeak<-mean(filter(bimodality_results,peaks=="2 Peaks" & best_model=="Bimodal")$prop_overlap,na.rm=T)
overlap_plot_twopeak<-ggplot(filter(bimodality_results,peaks=="2 Peaks" & best_model=="Bimodal"),aes(prop_overlap,1))+
  geom_jitter(size=0.6,height=0.4,colour="mediumpurple4")+
  geom_violin(alpha=0.5,fill="mediumpurple",colour="mediumpurple4")+
  geom_boxplot(width=0.1,fill="mediumpurple",alpha=0.7,colour="mediumpurple4",outlier.size = -1)+
  geom_vline(xintercept=overlap_mean_twopeak,colour="mediumpurple4",linetype="dashed")+
  theme_light()+
  labs(x="Overlap",y="Frequency")+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

################################################################################################################################################
# LINEAR MODELS
################################################################################################################################################

library(lme4)
library(MCMCglmm)
model_data<-bimodality_results %>% 
  left_join((data_ECOFF %>% dplyr::select(Bacteria_Antibiotic,class,class1,class2) %>% distinct()),by="Bacteria_Antibiotic") %>% 
  separate(Bacteria_Antibiotic, c("Bacteria", "Antibiotic"), sep="_") %>% 
  separate (Bacteria,c("genus","species"),sep=" ",extra="merge") %>% 
  mutate(genus=gsub(",", "",genus)) %>% 
  dplyr::select(class,Antibiotic,genus,species,prop_overlap) %>% 
  na.omit() %>% 
  mutate(logit_prop_overlap=ifelse(prop_overlap!=1,logit(prop_overlap),6))

model<-lmer(logit_prop_overlap~1+(1|species)+(1|Antibiotic),data=model_data)
model<-lmer(logit_prop_overlap~1+(1|genus)+(1|class),data=model_data)

nitt = 500000
burnin = 100000
thin = 200
mcmc_model <- MCMCglmm(logit_prop_overlap ~ 1, random = ~species + Antibiotic, data = model_data,nitt=nitt,burnin=burnin,thin=thin)
summary(mcmc_model)

library(tibble)
fixed_posterior<-as.data.frame(summary(mcmc_model)$solutions)
random_variance<-rbind(as.data.frame(summary(mcmc_model)$Rcovariances),as.data.frame(summary(mcmc_model)$Gcovariances)) %>% rownames_to_column()
random_variance$rowname <- factor(random_variance$rowname, levels = random_variance$rowname)
#ggplot(data=random_variance,aes(x=rowname,y=post.mean))+geom_bar(stat="identity")+theme_light()+labs(x="Predictor",y="Variance")+geom_errorbar(aes(ymin=`l-95% CI`, ymax=`u-95% CI`,group=rowname), width=0.2,position=position_dodge(0.9))
grid.arrange(
  ggplot(data=random_variance,aes(x=1,y=post.mean/sum(post.mean),fill=rowname))+
    geom_bar(stat="identity",position="stack",alpha=0.8)+
    theme_light()+
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    labs(x=NULL,y="Proportion Variance Explained",fill="Predictor")+
    scale_fill_manual(values=c("grey","mediumpurple","mediumseagreen"),labels=c("Residual","Bacteria","Antibiotic")),
ggplot(na.omit(blah),aes(prop_overlap,genus))+
  geom_boxplot(fill="mediumpurple",colour="black",alpha=0.8)+
  theme_light()+
  labs(title="Bacterial Genus",x=NULL,y=NULL)+
  theme(plot.title = element_text(hjust = 0.5)),
ggplot(na.omit(blah),aes(prop_overlap,class))+
  geom_boxplot(fill="mediumseagreen",colour="black",alpha=0.8)+
  theme_light()+
  labs(title="Antibiotic Class",x=NULL,y=NULL)+
  theme(plot.title = element_text(hjust = 0.5)),
nrow=1,bottom="Overlap Proportion")

################################################################################################################################################
# AMR BURDEN
################################################################################################################################################

bimodality_results<-bimodality_results %>% 
  separate(Bacteria_Antibiotic,c("Bacteria","Antibiotic"),sep="_",remove=F) %>% 
  full_join(read.csv("AMR Burden_bacteria.csv"),by="Bacteria") %>% 
  full_join(read.csv("AMR Burden_antibiotics.csv"),by="Antibiotic") %>% 
  full_join(read.csv("AMR Burden_data.csv") %>% gather("Antibiotic2","burden",2:16) %>% na.omit()) %>% 
  drop_na(overlap) %>% 
  mutate(group=ifelse(is.na(burden)==F,"Attributable Deaths Estimated","Attributable Deaths Not Estimated"))

wtest<-wilcox.test(overlap~ group,data = bimodality_results,exact = FALSE)
label<-paste0("Wilcoxon Test = ", round(wtest$statistic), ",  P-Value = ",round(wtest$p.value,4))
lab<-data.frame(x=1,y=2.55,label=label)
burden_violin<-ggplot(data=bimodality_results,aes(x=overlap,y=str_wrap(group,2)))+
  geom_violin(fill="mediumpurple",alpha=0.5,colour=NA)+
  geom_boxplot(fill="mediumpurple",alpha=0.5,width=0.15,colour="mediumpurple4")+
  geom_jitter(size=0.5,colour="mediumpurple4")+
  theme_light()+
  geom_label(data=lab,aes(x=x,y=y,label=label),hjust=1,vjust=1)+
  labs(x="Overlap Proportion",y=NULL)

burden<-bimodality_results %>% dplyr::select(Bacteria2,Antibiotic2,overlap,burden) %>% na.omit() %>%  mutate(Bacteria_Antibiotic=paste0(Bacteria2,"_",Antibiotic2)) %>% 
  group_by(Bacteria_Antibiotic) %>% summarise_at(vars(overlap), list(mean_overlap = mean)) %>% 
  full_join(bimodality_results %>% dplyr::select(Bacteria2,Antibiotic2,overlap,burden) %>% na.omit() %>%  mutate(Bacteria_Antibiotic=paste0(Bacteria2,"_",Antibiotic2))) %>% 
  dplyr::select(mean_overlap,burden,Bacteria_Antibiotic) %>% distinct()

stest<-cor.test(burden$mean_overlap, burden$burden, method = "spearman")
label<-paste0("Spearman Correlation = ", round(stest$estimate), ",  P-Value = ",round(stest$p.value,4))
lab<-data.frame(x=max(burden$mean_overlap),y=max(log10(burden$burden))+0.1,label=label)
burden_scatter<-ggplot(data=burden,aes(mean_overlap,log10(burden)))+
  geom_point(colour="mediumpurple4")+
  theme_light()+
  geom_label(data=lab,aes(x=x,y=y,label=label),hjust=1,vjust=0)+
  labs(x="Overlap Proportion",y="log10 (Attributable Deaths)")+
  ylim(2.5,5.3)

grid.arrange(burden_scatter,burden_violin,nrow=1)

################################################################################################################################################
# CSV OUTPUT
################################################################################################################################################
setwd("~/Downloads")
write.csv(bimodality_results %>% dplyr::select(Bacteria_Antibiotic,ECOFF,distributions,observations,resistance_freq,best_model,prob_bimodal,unimodalBIC,bimodalBIC,peaks,prop_overlap,diptest_pval,criticalbandwidth_pval,excessmass_pval), file="Tests for unimodality.csv")

setwd("~/OneDrive - University of Edinburgh/PhD Year 2/Quantitative Resistance Project_Final")
write.csv(bimodality_results,file="bimodality_results.csv",row.names = FALSE)
write.csv(data_ECOFF,file="data_ECOFF.csv",row.names = FALSE)
save(mcmc_model,file="mcmc_model.Rdata")
