################################################################################################################################################
# IMPORT CSV FILES
setwd("~/OneDrive - University of Edinburgh/Quantitative_Resistance/Quantitative_Resistance/CSVs")
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
setwd("~/OneDrive - University of Edinburgh/Quantitative_Resistance/Quantitative_Resistance")
library(dplyr)
library(tidyr)
################################################################################################################################################

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
library(diptest)
library(multimode)
library(mixR)
################################################################################################################################################

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

################################################################################################################################################
# PLOTS
library(gridExtra)
library(ggplot2)
library(patchwork)
library(ggpubr)
################################################################################################################################################

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


############### TEST P-VALUE PLOTS ###############

# Calculations
test_results<-bimodality_results %>% dplyr::select(criticalbandwidth_pval,excessmass_pval,diptest_pval) %>% gather("test","pval",1:3) %>% na.omit()

em_mean<-mean(filter(test_results,test=="excessmass_pval")$pval)
dt_mean<-mean(filter(test_results,test=="diptest_pval")$pval)
cb_mean<-mean(filter(test_results,test=="criticalbandwidth_pval")$pval)

em_bimodal<-nrow(filter(test_results,test=="excessmass_pval" & pval<0.05))/nrow(filter(test_results,test=="excessmass_pval"))
dt_bimodal<-nrow(filter(test_results,test=="diptest_pval" & pval<0.05))/nrow(filter(test_results,test=="diptest_pval"))
cb_bimodal<-nrow(filter(test_results,test=="criticalbandwidth_pval" & pval<0.05))/nrow(filter(test_results,test=="criticalbandwidth_pval"))

print(data.frame(Test=c("Excess Mass", "Dip Test","Critical Bandwidth"),Mean_Pvalue=c(em_mean,dt_mean,cb_mean),Proportion_Bimodal=c(em_bimodal,dt_bimodal,cb_bimodal)))

# Violin Plots - Pval
all_means<-data.frame(mean=c(cb_mean,dt_mean,em_mean),test=c("criticalbandwidth_pval","diptest_pval","excessmass_pval"),y=c(1,2,3),bimodal=c(cb_bimodal,dt_bimodal,em_bimodal))
pval_violin_plot<-ggplot(test_results,aes(pval,test))+
  geom_jitter(aes(colour=test),size=0.3,alpha=0.3,height=0.45)+
  geom_violin(aes(fill=test),colour=NA,alpha=0.6)+
  scale_fill_manual(values=c("mediumpurple","mediumseagreen","goldenrod"))+
  scale_colour_manual(values=c("mediumpurple","mediumseagreen","goldenrod"))+
  theme_light()+
  labs(title="Unimodality Tests",x="P-Value",y=NULL)+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
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


############### MIXTURE MODEL EXAMPLE PLOTS ###############

label<-c("Low Overlap", "High Overlap","1 Peak","1 Component")
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
    theme(plot.title = element_text(hjust = 0.5,size=10,colour = title_col[i]),plot.subtitle = element_text(hjust = 0.5,size=8),legend.position="none")
  
  example_results<-rbind(example_results,data.frame(Bacteria_Antibiotic=examples[i],
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


############### MIXTURE MODEL BAR PLOT ###############

plot_data<-data.frame(Group=c("Low Overlap","High Overlap","1 Peak","1 Component"),Count=c(
  nrow(bimodality_results %>% filter(best_model=="Bimodal", peaks=="2 Peaks",prop_overlap<0.01)),
  nrow(bimodality_results %>% filter(best_model=="Bimodal", peaks=="2 Peaks",prop_overlap>0.01)),
  nrow(bimodality_results %>% filter(best_model=="Bimodal", peaks=="1 Peak")),
  nrow(bimodality_results %>% filter(best_model=="Unimodal"))))

plot_data$Group <- factor(plot_data$Group, levels = plot_data$Group)

bar_plot<-ggplot(plot_data,aes(x=Group,y=Count))+
  geom_bar(stat="identity",alpha=0.8,colour="black",aes(fill=Group))+
  scale_fill_manual(values=c("mediumpurple","mediumseagreen","goldenrod","mediumvioletred"))+
  theme_light()+
  labs(title="Mixture Models")+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
  geom_label(aes(label=Count))

plot_data<-bimodality_results %>% mutate(Group=ifelse(best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap<0.01,"Low Overlap",
                                                      ifelse(best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap>0.01,"High Overlap",
                                                             ifelse(best_model=="Bimodal"& peaks=="1 Peak","1 Peak","1 Component"))))
ggplot(plot_data,aes(x=prop_overlap,y=Group))+geom_boxplot()

mean(filter(bimodality_results,best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap<0.01)$prop_overlap)

mean(filter(bimodality_results,best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap>0.01)$prop_overlap)
min(filter(bimodality_results,best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap>0.01)$prop_overlap)
max(filter(bimodality_results,best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap>0.01)$prop_overlap)

mean(filter(bimodality_results,best_model=="Bimodal"& peaks=="1 Peak")$prop_overlap)
min(filter(bimodality_results,best_model=="Bimodal"& peaks=="1 Peak")$prop_overlap)
max(filter(bimodality_results,best_model=="Bimodal"& peaks=="1 Peak")$prop_overlap)

mean(filter(bimodality_results,best_model=="Unimodal")$prop_overlap)


############### MIXTURE MODEL OVERLAP PLOT ###############

# Violin Plot - all
overlap_mean_all<-mean(bimodality_results$prop_overlap,na.rm=T)
excluding_examples_overlap<-bimodality_results %>% filter(!Bacteria_Antibiotic %in% examples)
examples_overlap<-bimodality_results %>% filter(Bacteria_Antibiotic %in% examples) %>% mutate(col=ifelse(Bacteria_Antibiotic==examples[1],"mediumpurple",ifelse(Bacteria_Antibiotic==examples[2],"mediumseagreen",ifelse(Bacteria_Antibiotic==examples[3],"goldenrod","mediumvioletred"))))

overlap_plot_all<-ggplot(excluding_examples_overlap,aes(prop_overlap,1))+
  geom_jitter(size=1,shape=16,height=0.4,colour="black",alpha=0.3)+
  geom_violin(data=bimodality_results,aes(x=prop_overlap,y=1),alpha=0.6,fill="grey",color=NA)+
  geom_segment(x = overlap_mean_all, y = 0.6, xend = overlap_mean_all, yend = 1.4,colour="grey")+
  geom_point(data=examples_overlap,aes(x=prop_overlap,y=1),colour=examples_overlap$col,size=4,shape=4,stroke=1.5)+
  theme_light()+
  labs(x="Overlap",y="Frequency")+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())

############### ARRANGING PLOTS ###############

examples<-ggarrange(example_plots[[1]],example_plots[[2]],example_plots[[3]],example_plots[[4]],nrow=1)
patchwork <- (pval_violin_plot / test_bar_plot) | (overlap_plot_all / examples / bar_plot) 
patchwork + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 1))

wrap_elements(grid::textGrob('Unimodality Tests')) / pval_violin_plot / test_bar_plot + plot_layout(heights = c(1,6,6)) | (wrap_elements(grid::textGrob('Mixture Models')) / examples /bar_plot / overlap_plot_all + plot_layout(heights = c(1,4,4,4)))
(pval_violin_plot / test_bar_plot) | (bar_plot / (example_plots[[1]]|example_plots[[2]]|example_plots[[3]]|example_plots[[4]]) / overlap_plot_all) 

################################################################################################################################################
# AMR BURDEN
################################################################################################################################################

bimodality_results2<-bimodality_results %>% 
  separate(Bacteria_Antibiotic,c("Bacteria","Antibiotic"),sep="_",remove=F) %>% 
  full_join(read.csv("AMR Burden_bacteria.csv"),by="Bacteria") %>% 
  full_join(read.csv("AMR Burden_antibiotics.csv"),by="Antibiotic") %>% 
  full_join(read.csv("AMR Burden_data.csv") %>% gather("Antibiotic2","burden",2:16) %>% na.omit()) %>% 
  drop_na(prop_overlap) %>% 
  mutate(group=ifelse(is.na(burden)==F,"Attributable Deaths Estimated","Attributable Deaths Not Estimated"))

wtest<-wilcox.test(prop_overlap ~ group,data = bimodality_results2,exact = FALSE)
label<-paste0("Wilcoxon Test = ", round(wtest$statistic), ",  P-Value = ",round(wtest$p.value,4))
lab<-data.frame(x=1,y=2.55,label=label)
burden_violin<-ggplot(data=bimodality_results2,aes(x=prop_overlap,y=str_wrap(group,2)))+
  geom_violin(fill="mediumpurple",alpha=0.5,colour=NA)+
  geom_boxplot(fill="mediumpurple",alpha=0.5,width=0.15,colour="mediumpurple4")+
  geom_jitter(size=0.5,colour="mediumpurple4")+
  theme_light()+
  #geom_label(data=lab,aes(x=x,y=y,label=label),hjust=1,vjust=1)+
  labs(x="Overlap Proportion",y=NULL)

burden<-bimodality_results2 %>% dplyr::select(Bacteria2,Antibiotic2,prop_overlap,burden) %>% na.omit() %>%  mutate(Bacteria_Antibiotic=paste0(Bacteria2,"_",Antibiotic2)) %>% 
  group_by(Bacteria_Antibiotic) %>% summarise_at(vars(prop_overlap), list(mean_overlap = mean)) %>% 
  full_join(bimodality_results2 %>% dplyr::select(Bacteria2,Antibiotic2,prop_overlap,burden) %>% na.omit() %>%  mutate(Bacteria_Antibiotic=paste0(Bacteria2,"_",Antibiotic2))) %>% 
  dplyr::select(mean_overlap,burden,Bacteria_Antibiotic) %>% distinct()

stest<-cor.test(burden$mean_overlap, burden$burden, method = "spearman")
label<-paste0("Spearman Correlation = ", round(stest$estimate), ",  P-Value = ",round(stest$p.value,4))
lab<-data.frame(x=max(burden$mean_overlap),y=max(log10(burden$burden))+0.1,label=label)
burden_scatter<-ggplot(data=burden,aes(mean_overlap,log10(burden)))+
  geom_point(colour="mediumpurple4")+
  theme_light()+
  #geom_label(data=lab,aes(x=x,y=y,label=label),hjust=1,vjust=0)+
  labs(x="Overlap Proportion",y="log10 (Attributable Deaths)")+
  ylim(2.5,5.3)

################################################################################################################################################
# LINEAR MODELS
library(MCMCglmm)
################################################################################################################################################

bimodality_results3<-bimodality_results %>% 
  separate(Bacteria_Antibiotic, c("Bacteria", "Antibiotic"), sep="_") %>% 
  mutate(logit_prop_overlap=ifelse(prop_overlap!=1,logit(prop_overlap),6))

nitt = 500000
burnin = 100000
thin = 200
mcmc_model <- MCMCglmm(logit_prop_overlap ~ 1, random = ~Bacteria + Antibiotic, data = bimodality_results3,nitt=nitt,burnin=burnin,thin=thin)
save(mcmc_model,file="mcmc_model.Rdata")

library(tibble)
fixed_posterior<-as.data.frame(summary(mcmc_model)$solutions)
random_variance<-rbind(as.data.frame(summary(mcmc_model)$Rcovariances),as.data.frame(summary(mcmc_model)$Gcovariances)) %>% rownames_to_column()
random_variance$rowname <- factor(random_variance$rowname, levels = random_variance$rowname)

mcmc_model_plot<-ggplot(data=random_variance,aes(x=1,y=post.mean/sum(post.mean),fill=rowname))+
  geom_bar(stat="identity",position="stack",alpha=0.8)+
  theme_light()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  labs(x=NULL,y="Proportion Variance Explained",fill="Predictor")+
  scale_fill_manual(values=c("grey","mediumpurple","mediumseagreen"),labels=c("Residual","Bacteria","Antibiotic"))

patchwork<-mcmc_model_plot|burden_violin |burden_scatter
patchwork + plot_layout(widths = c(1, 2,2)) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 1))
