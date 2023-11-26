setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance")
library(gridExtra)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyr)
library(dplyr)
library(stringr)
library(mixR)
library(MCMCglmm)
################################################################################################################################################
bimodality_results<-read.csv("bimodality_results.csv")
data_ECOFF<-read.csv("data_ECOFF.csv")

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
  scale_fill_manual(values=c("darkslategray2","darkslategray3","darkslategray4"))+
  scale_colour_manual(values=c("darkslategray2","darkslategray3","darkslategray4"))+
  theme_light()+
  labs(title="Unimodality Tests",x="P-Value",y=NULL)+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))+
  geom_segment(aes(x = mean, y = y-0.45, xend = mean, yend = y+0.45,colour=test),data=all_means)+
  geom_vline(xintercept=0.05,linetype="dashed")+
  scale_y_discrete(labels=c("excessmass_pval" = str_wrap("Excess Mass",5), "diptest_pval" = "Dip Test","criticalbandwidth_pval" = str_wrap("Critical Bandwidth",5)))

prop_bimodal_tests<-data.frame(table((filter(test_results,pval<0.05))$test)) %>% full_join(data.frame(table(test_results$test)),by="Var1") %>% mutate(prop_bimodal=Freq.x/Freq.y)
test_bar_plot<-ggplot(prop_bimodal_tests,aes(x=prop_bimodal,y=Var1))+
  geom_bar(stat="identity",aes(fill=Var1),alpha=0.8,colour="black")+
  theme_light()+
  theme(legend.position="none")+
  scale_fill_manual(values=c("darkslategray2","darkslategray3","darkslategray4"))+
  labs(x="Proportion Bimodal (P-Value < 0.05)",y=NULL)+
  scale_y_discrete(labels=c("excessmass_pval" = str_wrap("Excess Mass",5), "diptest_pval" = "Dip Test","criticalbandwidth_pval" = str_wrap("Critical Bandwidth",5)))

############### MIXTURE MODEL EXAMPLE PLOTS ###############

label<-c("Low Overlap", "High Overlap","1 Peak","1 Component")
title<-c("S.agalactiae \n+ CIP", "K.pneumoniae \n+ CST", "E.coli \n+ C/T", "C.jejuni \n+ SXT")
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
    geom_line(data=df_density,aes(x=x,y=y_s),colour="dodgerblue4",size=0.5) + 
    geom_line(data=df_density,aes(x=x,y=y_r),colour="firebrick4",size=0.5) + 
    #geom_line(data=df_density,aes(x=x,y=y_all),linetype="dashed",size=0.5) + 
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
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))

plot_data<-bimodality_results %>% mutate(Group=ifelse(best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap<0.01,"Low Overlap",
                                                      ifelse(best_model=="Bimodal"& peaks=="2 Peaks"&prop_overlap>0.01,"High Overlap",
                                                             ifelse(best_model=="Bimodal"& peaks=="1 Peak","1 Peak","1 Component"))))

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
  labs(title="Mixture Models",x="Overlap",y="Frequency")+
  theme(plot.title = element_text(hjust = 0.5),axis.text.y = element_blank(),axis.ticks.y = element_blank())

############### ARRANGING PLOTS ###############

examples<-ggarrange(example_plots[[1]],example_plots[[2]],example_plots[[3]],example_plots[[4]],nrow=1)
patchwork <- (pval_violin_plot / test_bar_plot) | (overlap_plot_all / examples / bar_plot + plot_layout(heights = c(2,1,2))) 
patchwork + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 1), text = element_text(size = 14))

################################################################################################################################################
# PROPORTION VARIANCE IN OVERLAP EXPLAINED BY ANTIBIOTIC AND BACTERIA

bimodality_results3<-bimodality_results %>% 
  separate(Bacteria_Antibiotic, c("Bacteria", "Antibiotic"), sep="_") %>% 
  mutate(logit_prop_overlap=ifelse(prop_overlap!=1,logit(prop_overlap),6))

nitt = 500000
burnin = 100000
thin = 200
mcmc_model <- MCMCglmm(logit_prop_overlap ~ 1, random = ~Bacteria + Antibiotic, data = bimodality_results3,nitt=nitt,burnin=burnin,thin=thin)
save(mcmc_model,file="mcmc_model.Rdata")
load("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance/mcmc_model.Rdata")

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

mcmc_model_plot & theme(text = element_text(size = 18))

################################################################################################################################################
# AMR BURDEN

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
  geom_violin(fill="darkslategray3",alpha=0.5,colour=NA)+
  geom_boxplot(fill="darkslategray3",alpha=0.5,width=0.15,colour="darkslategray4")+
  geom_jitter(size=0.5,colour="darkslategray4")+
  theme_light()+
  labs(x="Overlap Proportion",y=NULL)

burden<-bimodality_results2 %>% dplyr::select(Bacteria2,Antibiotic2,prop_overlap,burden) %>% na.omit() %>%  mutate(Bacteria_Antibiotic=paste0(Bacteria2,"_",Antibiotic2)) %>% 
  group_by(Bacteria_Antibiotic) %>% summarise_at(vars(prop_overlap), list(mean_overlap = mean)) %>% 
  full_join(bimodality_results2 %>% dplyr::select(Bacteria2,Antibiotic2,prop_overlap,burden) %>% na.omit() %>%  mutate(Bacteria_Antibiotic=paste0(Bacteria2,"_",Antibiotic2))) %>% 
  dplyr::select(mean_overlap,burden,Bacteria_Antibiotic) %>% distinct()

stest<-cor.test(burden$mean_overlap, burden$burden, method = "spearman")
label<-paste0("Spearman Correlation = ", round(stest$estimate), ",  P-Value = ",round(stest$p.value,4))
lab<-data.frame(x=max(burden$mean_overlap),y=max(log10(burden$burden))+0.1,label=label)
burden_scatter<-ggplot(data=burden,aes(mean_overlap,log10(burden)))+
  geom_point(colour="darkslategray4")+
  theme_light()+
  labs(x="Overlap Proportion",y="log10 (Attributable Deaths)")+
  ylim(2.5,5.3)

(burden_violin |burden_scatter) + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 1), text = element_text(size = 16))