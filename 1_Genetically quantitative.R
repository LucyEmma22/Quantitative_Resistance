setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(stringr)
library(seqinr)
library(easyPubMed)
library(patchwork)
library(gridExtra)
library(ggpubr)
library(ggtext)
#########################################################################################################
# GWAS
#########################################################################################################

GWAS_data<-read.csv("GWAS_data.csv") %>% filter(value!=0)

#########################################################################################################
# GWAS GENES AND MUTATIONS

GWAS_mutations<-filter(GWAS_data,gene_mutation=="mutation",value!=0)
GWAS_genes_plot<-ggplot(GWAS_genes,aes(x=log10(value),y=method))+
  geom_boxplot(fill="mediumseagreen",alpha=0.8)+ 
  theme_light()+ 
  labs(x="Log10 (Number of Genes)",y="bGWAS Method")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept=log10(mean(GWAS_genes$value,na.rm=T)),colour="mediumseagreen",linetype="dashed")

GWAS_genes<-filter(GWAS_data,gene_mutation=="gene",value!=0)
GWAS_mutations_plot<-ggplot(GWAS_mutations,aes(x=log10(value),y=method))+
  geom_boxplot(fill="mediumpurple",alpha=0.8)+ 
  theme_light()+ 
  labs(x="Log10 (Number of Mutations)",y="bGWAS Method")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept=log10(mean(GWAS_mutations$value,na.rm=T)),colour="mediumpurple",linetype="dashed")

#########################################################################################################
# HERITABILITY

heritability<-read.csv("heritability.csv")
heritability_plot<-ggplot(heritability,aes(x=value,y=antibiotic,fill=factor(Heritability, levels=c("h2_missing","h2_explained"))))+
  geom_bar(stat="identity",colour="black",position="stack")+
  scale_fill_manual(name=NULL,values=c("white","grey"),labels=c("Missing","Explained"))+
  theme_light()+
  labs(y="Antibiotic",x="Heritability")+
  facet_wrap(~bacteria,scales="free")+
  theme(legend.position="bottom")

#########################################################################################################
# Distribution of Sample Sizes

ggplot(GWAS_data,aes(log10(n)))+
  geom_histogram(bins=20,boundary=0,fill="grey")+
  theme_light()+
  labs(title="GWAS (Distribution of Sample Sizes)",x="Log10 (Sample Size)",y="Frequency")+
  theme(plot.title = element_text(hjust = 0.5))

#########################################################################################################
# RESFINDER
#########################################################################################################

# RESFINDER GENES

setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance/ResFinder_ARGs")
files<-list.files()
files<-files[grepl( ".fsa", files, fixed = TRUE)]
drug_list<-gsub('.fsa','',files)
gene_nos<-data.frame()
for (i in 1:length(files)){
  data<-read.fasta(file = files[i], as.string = TRUE)
  gene_nos<-rbind(gene_nos,data.frame(antibiotic=str_to_sentence(drug_list[i]),genes_mutated=length(data)))
}

gene_resfinder<-ggplot(gene_nos,aes(x=log10(genes_mutated),y=antibiotic))+
  geom_bar(stat="identity",colour="black",fill="mediumseagreen",alpha=0.8)+ 
  theme_light()+ 
  labs(x="Log10 (Number of ARGs)",y="Antibiotic")+ 
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept=log10(mean(gene_nos$genes_mutated)),colour="mediumseagreen",linetype="dashed")

# RESFINDER GENES VS PUBMED SEARCH

pubmed_results1<-data.frame()
for (i in 1:nrow(gene_nos)){
  pubmed_results1<-pubmed_results1 %>% 
    rbind(data.frame(antibiotic=gene_nos$antibiotic[i],pubmed=as.numeric(get_pubmed_ids(paste0("(drug resistance, bacterial[MeSH Terms]) AND ",gene_nos$antibiotic[i]))$Count)))
}
gene_nos<-left_join(gene_nos,pubmed_results1)

model<-lm(log10(genes_mutated)~log10(pubmed),data=gene_nos)
c<-model$coef[1]
m<-model$coef[2]
x=log10(seq(min(gene_nos$pubmed),max(gene_nos$pubmed),length=100))
df<-data.frame(x=x,y=m*x+c)
title<-paste0("Estimate = ", round(m,2), ",  P-Value = ",round(summary(model)$coefficients[2,4],3))

gene_resfinder_pubmed<-ggplot(gene_nos,aes(x=log10(pubmed),y=log10(genes_mutated)))+
  geom_point(colour="mediumseagreen")+
  theme_light()+
  labs(x="Log10 (Number of PubMed Search Results)", y="Log10 (Number of ARGs)")+
  geom_line(data=df,aes(x=x,y=y),colour="mediumseagreen")
  #geom_text(data=data.frame(x=2.5,y=3,lab=title),aes(x=x,y=y,label=title),hjust=0,size=3)

setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance")
write.csv(gene_nos,file="ResFinder_ARGs_and_PubMed_Search_Results.csv")

# RESFINDER MUTATIONS

setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance/PointFinder_mutations")
folders<-list.files()
mut_nos<-data.frame()

for (i in 1:length(folders)){
  setwd(paste("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance/PointFinder_mutations",folders[i],sep="/"))
  data<-read.delim2("resistens-overview.txt")[2:7] %>% 
    mutate(Resistance = strsplit(as.character(Resistance), ",")) %>% unnest(Resistance) %>% mutate(Resistance=trimws(Resistance)) %>% mutate(Resistance=str_to_sentence(Resistance)) %>% 
    filter(is.na(as.numeric(Codon_pos))==F) %>% filter(is.na(as.numeric(Ref_nuc))==T)
  indel<-data %>% filter(Ref_codon=="ins"|Ref_codon=="del") %>% mutate(Res_codon = strsplit(as.character(Res_codon), ",")) %>% unnest(Res_codon) %>% mutate(Res_codon=trimws(Res_codon))
  snp<-data %>%  filter(Ref_codon!="ins"&Ref_codon!="del") %>% mutate(Res_codon = strsplit(as.character(Res_codon), ",")) %>% unnest(Res_codon) %>% mutate(Res_codon=trimws(Res_codon)) 
  data<-data.frame(rbind(indel,snp))
  mut_nos<-rbind(mut_nos,
                 data.frame(table(data$Resistance)) %>% rename(codons=Freq) %>%
                   full_join(data.frame(table(distinct(dplyr::select(data,Gene_name,Codon_pos,Resistance))$Resistance)),by="Var1") %>% rename (mutations=Freq) %>%
                   full_join(data.frame(table(distinct(dplyr::select(data,Gene_name,Resistance))$Resistance)),by="Var1") %>% rename(genes=Freq) %>% rename(antibiotic=Var1) %>% mutate(bacteria=gsub("_"," ",str_to_title(folders[i])))
  )
}
mut_nos$antibiotic<-gsub("Gentamicin c","Gentamicin",mut_nos$antibiotic)
mut_nos$antibiotic<-gsub("Para-aminosalicyclic acid","Para-aminosalicylic acid",mut_nos$antibiotic)

mut_resfinder<-ggplot(mut_nos,aes(x=log10(mutations),y=str_wrap(bacteria,5)))+
  geom_boxplot(fill="mediumpurple",alpha=0.8)+ 
  theme_light()+ 
  labs(x="Log10 (Number of Chromosomal Mutations)",y="Bacterial Species")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept=log10(mean(mut_nos$mutations)),colour="mediumpurple",linetype="dashed")

# RESFINDER MUTATIONS VS PUBMED SEARCH 

pubmed_results2<-data.frame()
for (i in 1:nrow(mut_nos)){
  pubmed_results2<-pubmed_results2 %>% 
    rbind(data.frame(antibiotic=mut_nos$antibiotic[i],bacteria=mut_nos$bacteria[i],pubmed=as.numeric(get_pubmed_ids(paste0("(drug resistance, bacterial[MeSH Terms]) AND ",mut_nos$antibiotic[i]," AND ",mut_nos$bacteria[i]))$Count)))
}
mut_nos<-left_join(mut_nos,pubmed_results2)

model<-lm(log10(mutations)~log10(pubmed),data=mut_nos)
c<-model$coef[1]
m<-model$coef[2]
x=log10(seq(min(mut_nos$pubmed),max(mut_nos$pubmed),length=100))
df<-data.frame(x=x,y=m*x+c)
title<-paste0("Estimate = ", round(m,2), ",  nP-Value = ",round(summary(model)$coefficients[2,4],3))

mut_resfinder_pubmed<-ggplot(mut_nos,aes(log10(pubmed),log10(mutations)))+
  geom_point(colour="mediumpurple")+
  theme_light()+
  labs(x="Log10 (Number of PubMed Search Results)", y="Log10 (Number of Chromosomal Mutations)")+
  geom_line(data=df,aes(x=x,y=y),colour="mediumpurple")
  #geom_text(data=data.frame(x=2.5,y=3,lab=title),aes(x=x,y=y,label=title),hjust=0,size=3)

setwd("~/Library/CloudStorage/OneDrive-UniversityofEdinburgh/Quantitative_Resistance/Quantitative_Resistance")
write.csv(mut_nos,file="PointFinder_Mutations_and_PubMed_Search_Results.csv")

#########################################################################################################
# ARRANGE FIGURE
#########################################################################################################

patchwork <- ((gene_resfinder|gene_resfinder_pubmed) / (mut_resfinder|mut_resfinder_pubmed) / (GWAS_genes_plot|GWAS_mutations_plot) / heritability_plot) + plot_layout(widths = c(1,1,2))
patchwork + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(face = 1), text = element_text(size = 8)) 

#########################################################################################################
# DISPLAY MEAN/MIN/MAX 
#########################################################################################################

heritability_2 <- heritability %>% spread(Heritability,value) %>% mutate(h2_total=h2_missing+h2_explained) %>% mutate(prop_explained=h2_explained/h2_total)

# Rounded
print(data.frame(Name=c("GWAS Mutations","GWAS Genes","ResFinder Genes","ResFinder Mutations","Total Heritability","Proprotion Heritability Explained"),
                 Mean=c(round(mean(GWAS_mutations$value,na.rm=T),2),round(mean(GWAS_genes$value,na.rm=T),2),round(mean(gene_nos$genes_mutated),2),round(mean(mut_nos$mutations),2),round(mean(heritability_2$h2_total),3),round(mean(heritability_2$prop_explained),3)),
                 Min=c(min(GWAS_mutations$value,na.rm=T),min(GWAS_genes$value,na.rm=T),min(gene_nos$genes_mutated),min(mut_nos$mutations),min(heritability_2$h2_total),round(min(heritability_2$prop_explained),3)),
                 Max=c(max(GWAS_mutations$value,na.rm=T),max(GWAS_genes$value,na.rm=T),max(gene_nos$genes_mutated),max(mut_nos$mutations),max(heritability_2$h2_total),round(max(heritability_2$prop_explained),3))
))

# Unrounded
print(data.frame(Name=c("GWAS Mutations","GWAS Genes","ResFinder Genes","ResFinder Mutations","Total Heritability","Proprotion Heritability Explained"),
                 Mean=c(mean(GWAS_mutations$value,na.rm=T),mean(GWAS_genes$value,na.rm=T),mean(gene_nos$genes_mutated),mean(mut_nos$mutations),mean(heritability_2$h2_total),mean(heritability_2$prop_explained)),
                 Min=c(min(GWAS_mutations$value,na.rm=T),min(GWAS_genes$value,na.rm=T),min(gene_nos$genes_mutated),min(mut_nos$mutations),min(heritability_2$h2_total),min(heritability_2$prop_explained)),
                 Max=c(max(GWAS_mutations$value,na.rm=T),max(GWAS_genes$value,na.rm=T),max(gene_nos$genes_mutated),max(mut_nos$mutations),max(heritability_2$h2_total),max(heritability_2$prop_explained))
))              