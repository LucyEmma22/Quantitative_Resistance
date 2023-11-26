setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)

all_results<-read.csv("within_host_results.csv")
examples_data<-read.csv("within_host_results_examples.csv")

#######################################################################################################################
# Example dynamics

cols <- c("dodgerblue","firebrick")
colGradient <- colorRampPalette(cols)
examples_data$facet <- factor(examples_data$facet, levels = unique(examples_data$facet))
example_dynamics<-ggplot(examples_data,aes(x=time,y=frequency,colour=as.factor(round(mic,1))))+
  geom_line(key_glyph = "rect", size=0.2)+
  theme_light()+
  scale_y_log10(limits = c(1,10^10))+
  facet_wrap(~facet,nrow=1)+
  scale_colour_manual(name="MIC",values=colGradient(mic))+
  theme(legend.position="bottom")+
  guides(colour = guide_legend(nrow = 1))+
  labs(x="Time (Days)",y="Density of Bacteria")

#######################################################################################################################
# Resistant cell density for n=2 and n=10

total<-all_results %>% group_by(n,dose,run) %>% summarise_at(vars(freq),list(total = sum)) %>% group_by(n,dose) %>% summarise_at(vars(total),list(mean_total = mean))
resistant<-all_results %>% filter(mic>2) %>% group_by(n,dose,run) %>% summarise_at(vars(freq),list(resistant = sum)) %>% group_by(n,dose) %>% summarise_at(vars(resistant),list(mean_resistant = mean))
total_and_resistant<-full_join(total,resistant,by=c("n","dose")) %>% mutate(mean_prop_resistant=mean_resistant/mean_total)

twoten<-filter(total_and_resistant,n==2|n==10)
twoten %>% filter(mean_resistant==max(mean_resistant))
twoten_plot<-ggplot(twoten,aes(x=dose,y=mean_resistant))+
  geom_line(size=0.8)+
  geom_line(data=twoten,aes(x=dose,y=mean_prop_resistant*max(twoten$mean_total)),colour="mediumpurple",linetype="dashed",size=0.5)+
  geom_line(data=twoten,aes(x=dose,y=mean_total),colour="mediumseagreen",linetype="dashed",size=0.5)+
  facet_wrap(.~n)+
  labs(x="Antibiotic Dose",y="Total Resistant Cell Density")+
  theme_light()+
  geom_vline(xintercept=c(20,35),linetype="solid",size=0.6,colour="grey")

#######################################################################################################################
# Best Dose for different n and clinical ranges

df<-data.frame(min=rep(seq(5,35,5),each=7),max=rep(seq(10,40,5),7)) %>% filter(max>min)
best_dose<-data.frame()
for(i in 1:nrow(df)){
  best_dose<-rbind(best_dose,filter(total_and_resistant,dose==df$min[i]|dose==df$max[i]) %>% group_by(n) %>%slice(which.min(mean_resistant)) %>% dplyr::select(n,dose) %>% mutate(min=df$min[i],max=df$max[i]))
}
best_dose <- best_dose %>% mutate(best=ifelse(dose==min,"Low","High"))
adjust<-1
best_dose_adj <- best_dose %>% mutate(min=ifelse(n==2|n==3|n==4,min+adjust,min)) %>% 
  mutate(min=ifelse(n==8|n==9|n==10,min-adjust,min)) %>% 
  mutate(max=ifelse(n==2|n==5|n==8,max-adjust,max)) %>% 
  mutate(max=ifelse(n==4|n==7|n==10,max+adjust,max))

best_dose_plot<-ggplot(best_dose_adj,aes(x=max,y=min))+
  geom_point(aes(fill=best),shape=22,size=6,colour="black",stroke=0.2)+
  geom_point(data=filter(best_dose_adj,n==2 & min==20+adjust & max==35-adjust),aes(x=max,y=min),shape=22,size=6,colour="red",stroke=0.9)+
  geom_point(data=filter(best_dose_adj,n==10 & min==20-adjust & max==35+adjust),aes(x=max,y=min),shape=22,size=6,colour="red",stroke=0.9)+
  theme_light()+
  labs(fill="Best Dose", x="Maximum Clinical Dose",y="Minimum Clinical Dose")+
  scale_fill_manual(values=c("dodgerblue","goldenrod"))+
  geom_text(aes(label=n),colour="black",size=3)

#######################################################################################################################
# Arrange plots

patchwork<-example_dynamics/twoten_plot / best_dose_plot  + plot_layout(heights = c(1,2,4))
patchwork + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 1), text = element_text(size = 12))