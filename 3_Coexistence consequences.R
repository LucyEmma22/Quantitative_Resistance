setwd("~/OneDrive - University of Edinburgh/PhD Year 2/Toy Model")

library(dplyr)
library(stringr)

#####################################################################################################################################
# CALCULATING FITTEST MIC AND PROPORTION OF RESISTANCE FOR DIFFERENT ANTIBIOTIC USAGE RATES AND AMOUNTS OF STANDING GENETIC VARIANCE
#####################################################################################################################################

sd_rs_list<-c(0.001, 0.5, 1) # Standing genetic variance 

ts_list<-seq(0,20,length.out=1001) # Antibiotic usage x antibiotic clearance rate

# Data frame for standing genetic variance and treatment rate
sd_rs_rep<-data.frame(sd_rs=rep(sd_rs_list,each=length(ts_list)))
ts_rep<-data.frame(ts=rep(ts_list,times=length(sd_rs_list))) 
SGV_and_treatment<-cbind(sd_rs_rep,ts_rep) 

c<-0.3 # cost of resistance 
b<-1 # transmission rate
g<-0.5 # immune clearance rate
n<-1 # host population density
cutoff<-2 # MIC threshold for resistance

df<-data.frame()
for (i in 1:nrow(SGV_and_treatment)){
  t<-SGV_and_treatment$ts[i] # Treatment rate
  sd_r<-SGV_and_treatment$sd_rs[i] # Standing genetic variance
  
  # R0 = stuff_in / stuff_out
  # R0 = (transmission x host_pop)  /  ( cost x MIC  +  immune_clearance  +  AB_clearance x treatment_rate x exp(-MIC) )
  # R0 = (beta*n)/(c*MIC + gamma + t*exp(-MIC))
  # Differentiate R0 with respect to MIC (how does R0 with change with change in MIC)
  # dR0/dMIC = -(n*beta*(c-t*exp(-MIC))) / (c*MIC + gamma + t*exp(-MIC))^2
  # Solve for MIC when R0 is maximum (dR0/dMIC=0)
  # MIC = -log(c/t) = log(t/c)
  fittest_mic<-ifelse(log(t/c)>0, log(t/c), 0) # Fittest MIC (model with clearance cost)

  prop_resistant<-pnorm(log(cutoff),mean=log(fittest_mic),sd=sd_r,lower.tail=F) # Proportion of strains above threshold for resistance
  # Normal distribution, mean = log(MIC), SD = standing genetic variance (0, 0.1, 1 or 10)
  # we are saying that the strains in circulation are normally distributed around the fittest MIC, with standing genetic variance
  # given our distribution of MIC (mean = fittest MIC, SD = standing genetic variance), what is the probability of drawing a number above the cutoff (threshold for resistance)

  x<-data.frame(t,sd_r,fittest_mic,prop_resistant)
  df<-rbind(df,x)
}

######################################################################################################################################################
# GRAPHS
library(gridExtra)
library(ggplot2)
library(patchwork)
library(ggpubr)
######################################################################################################################################################

tr_examlpes<-c(1,3,5)
tr_titles<-c("I","II","III")
lab<-filter(df, t %in% tr_examlpes) %>%full_join(data.frame(t=tr_examlpes,label=tr_titles),by="t")
  
fittest_mic_plot<-ggplot(df,aes(x=t,y=fittest_mic))+
  geom_line(lwd=0.8)+
  geom_label(data=lab,aes(x=t,y=fittest_mic,label=label))+
  labs(x="Treatment Rate",y="Fittest MIC")+
  theme_light()

prop_resistance_plot<-ggplot(df,aes(x=t,y=prop_resistant,colour=as.factor(sd_r)))+
  geom_line(lwd=0.8)+
  theme_light()+
  labs(x="Treatment Rate",y="Proportion Resistant")+
  scale_colour_manual(name=str_wrap("Standing Genetic Variation", width=30),values=c("mediumpurple","mediumseagreen","goldenrod"))+
  geom_label(data=lab,aes(x=t,y=prop_resistant,label=label), show.legend = F, label.padding = unit(0.1, "lines"))+
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid", colour ="black"),legend.position="bottom")
  
# Examples 
fittest_mic<-ifelse(log(tr_examlpes/c)>0, log(tr_examlpes/c), 0)
to_plot<-data.frame(fittest_mic=rep(fittest_mic,length(sd_rs_list)),variance=rep(sd_rs_list,each=length(tr_examlpes)))
title<-rep(tr_titles,length(sd_rs_list))
title_col<-rep(c("mediumpurple","mediumseagreen","goldenrod"),each=length(sd_rs_list))

plots<-list()
for(i in 1:nrow(to_plot)){
mean<-to_plot$fittest_mic[i]
variance<-to_plot$variance[i]
df2<-data.frame(x=seq(-4,8,0.01),y=dnorm(seq(-4,8,0.01), mean, variance))
plots[[i]]<-ggplot(data=df2,aes(x,y))+
  geom_ribbon(data=subset(df2,x>2),aes(ymax=y),ymin=0,fill="firebrick3",colour="firebrick3")+
  geom_ribbon(data=subset(df2,x<2),aes(ymax=y),ymin=0,fill="dodgerblue",colour="dodgerblue")+
  #geom_line(data=df2,aes(y=y,x=x),fill="white",colour="black")+
  theme_void()+
  theme(axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(colour = title_col[i],hjust=0.5,face="bold"))+
  geom_vline(xintercept = 2, linetype="dashed", color = "black", size=0.5)+
  labs(title=title[i],x=NULL,y=NULL)
}

example_plots<-ggarrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],plots[[9]])
example_plots<-annotate_figure(example_plots,bottom = text_grob("MIC"),left = text_grob("Frequency",rot = 90))
patchwork<-(fittest_mic_plot|prop_resistance_plot|example_plots) + plot_layout(widths = c(3,3,3))
patchwork + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 1))
