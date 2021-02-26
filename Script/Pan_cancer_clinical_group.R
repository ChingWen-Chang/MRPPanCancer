
setwd("/Users/ching-wenchang/Documents/Sur")
clinical_stage<-read.delim("/Users/ching-wenchang/Documents/Sur/clinical_stage_MRP_group.txt",header = T)
clinical_stage[1:10,1:10]

library("gplots")
library(tidytext)
library(janeaustenr)
library(dplyr)
library(ggplot2)

####AJCC_stage
result<-clinical_stage %>% 
  group_by(Group) %>% count(stage)
result<-result[-which(is.na(result$stage)),]
result<-result[-which(result$stage==0),]
result$n <- ave(result$n, result$stage, FUN=function(x) x/max(x)) 
stage_plot<-ggplot(result, aes(fill=as.factor(stage), y=n, x=Group)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette = "PuRd")+ labs(fill = "Factor_Factor_Factor")+ ggtitle("Stage") 
pdf(file="TCGA_MRP_stage_plot_nor.pdf", width=3.5, height=2.5, compress=TRUE, useDingbats=F)
stage_plot
dev.off()
write.table(result,"clinical_stage_MRP_group_n.txt",quote = F,sep = "\t")

clinical_stage_MRP_group_n<-read.delim("/Users/ching-wenchang/Documents/Sur/clinical_stage_MRP_group_n_1.txt",header = T)
stage <- as.table(as.matrix(clinical_stage_MRP_group_n))
rownames(stage)<-stage[ ,1]
stage<-stage[ ,-1]
stage <- as.table(as.matrix(clinical_stage_MRP_group_n))
rownames(stage)<-stage[ ,1]
stage<-stage[ ,-1]
stage<-stage[,-c(2,3)]
stage_1<-apply(matrix(stage,ncol=2,nrow=4),1,as.numeric)
dimnames(stage_1) <- list(Group = c("WT", "LD"),Stage = c("Stage1", "Stage2", "Stage3", "Stage4"))
chisq<- chisq.test(unlist(t(stage_1)))

####Age
result_Age<-clinical_stage %>% 
  group_by(Group) %>% count(Age)
result_Age<-result_Age[-which(is.na(result_Age$Age)),]
result_Age$n <- ave(result_Age$n, result_Age$Age, FUN=function(x) x/max(x)) 
Age_plot<-ggplot(result_Age, aes(fill=as.factor(Age), y=n, x=Group)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette = "RdPu")+ labs(fill = "Factor_Factor_Factor")+ ggtitle("Age") 
pdf(file="TCGA_MRP_Age_plot_nor.pdf", width=3.5, height=2.5, compress=TRUE, useDingbats=F)
Age_plot
dev.off()


clinical_Age_MRP_group_n<-read.delim("/Users/ching-wenchang/Documents/Sur/clinical_Age_MRP_group_n.txt",header = T)
Age <- as.table(as.matrix(clinical_Age_MRP_group_n))
rownames(Age)<-Age[ ,1]
Age<-Age[ ,-1]
Age <- as.table(as.matrix(clinical_Age_MRP_group_n))
rownames(Age)<-Age[ ,1]
Age<-Age[ ,-1]
Age<-Age[,-c(2,3)]
Age_1<-apply(matrix(Age,ncol=2,nrow=2),1,as.numeric)
dimnames(Age_1) <- list(Group = c("WT", "LD"),Age = c("<65", ">65"))
chisq<- chisq.test(unlist(t(Age_1)))


####gender
result_gender<-clinical_stage %>% 
  group_by(Group) %>% count(gender)
result_gender$n <- ave(result_gender$n, result_gender$gender, FUN=function(x) x/max(x)) 
gender_plot<-ggplot(result_gender, aes(fill=as.factor(gender), y=n, x=Group)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette = "YlOrBr")+  labs(fill = "Factor_Factor_Factor")+ ggtitle("Gender") 
pdf(file="TCGA_MRP_gender_plot_nor.pdf", width=3.5, height=2.5, compress=TRUE, useDingbats=F)
gender_plot
dev.off()


clinical_gender_MRP_group_n<-read.delim("/Users/ching-wenchang/Documents/Sur/clinical_gender_MRP_group_n.txt",header = T)
gender <- as.table(as.matrix(clinical_gender_MRP_group_n))
rownames(gender)<-gender[ ,1]
gender<-gender[ ,-1]
gender<-gender[,-c(2,3)]
gender_1<-apply(matrix(gender,ncol=2,nrow=2),1,as.numeric)
dimnames(gender_1) <- list(Group = c("WT", "LD"),gender = c("Female", "Male"))
chisq<- chisq.test(unlist(t(gender_1)))

####race
result_race<-clinical_stage %>% 
  group_by(Group) %>% count(race)
result_race<-result_race[-which(is.na(result_race$race)),]
result_race$n <- ave(result_race$n, result_race$race, FUN=function(x) x/max(x)) 
race_plot<-ggplot(result_race, aes(fill=as.factor(race), y=n, x=Group)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette = "Pastel1")+ labs(fill = "Factor_Factor_Factor")+ ggtitle("Race") 

pdf(file="TCGA_MRP_race_plot_nor.pdf", width=5.3, height=2.5, compress=TRUE, useDingbats=F)
race_plot
dev.off()


clinical_race_MRP_group_n<-read.delim("/Users/ching-wenchang/Documents/Sur/clinical_race_MRP_group_n.txt",header = T)
race <- as.table(as.matrix(clinical_race_MRP_group_n))
rownames(race)<-race[ ,1]
race<-race[ ,-1]
race<-race[,-c(2,3)]
race_1<-apply(matrix(race,ncol=2,nrow=3),1,as.numeric)
dimnames(race_1) <- list(Group = c("WT", "LD"),race = c("ASIAN", "BLACK OR AFRICAN AMERICAN", "WHITE"))
chisq<- chisq.test(unlist(t(race_1)))



####treatment_outcome
result_treatment_outcome<-clinical_stage %>% 
  group_by(Group) %>% count(treatment_outcome)
result_treatment_outcome<-result_treatment_outcome[-which(is.na(result_treatment_outcome$treatment_outcome)),]
result_treatment_outcome$treatment_outcome1<-factor(result_treatment_outcome$treatment_outcome, levels=c("Complete Remission/Response", "Partial Remission/Response", "Stable Disease", "Progressive Disease"))
result_treatment_outcome$n <- ave(result_treatment_outcome$n, result_treatment_outcome$treatment_outcome, FUN=function(x) x/max(x)) 
treatment_outcome_plot<-ggplot(result_treatment_outcome, aes(fill=as.factor(treatment_outcome1), y=n, x=Group)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette = "Purples")+ labs(fill = "Factor_Factor_Factor")+ ggtitle("Treatment_outcome") 

pdf(file="TCGA_MRP_treatment_outcome_plot_nor.pdf", width=4, height=2.5, compress=TRUE, useDingbats=F)
treatment_outcome_plot
dev.off()

clinical_treatment_outcome_MRP_group_n<-read.delim("/Users/ching-wenchang/Documents/Sur/clinical_treatment_outcome_MRP_group_n.txt",header = T)
treatment_outcome <- as.table(as.matrix(clinical_treatment_outcome_MRP_group_n))
rownames(treatment_outcome)<-treatment_outcome[ ,1]
treatment_outcome<-treatment_outcome[ ,-1]
treatment_outcome<-treatment_outcome[,-c(2,3)]
treatment_outcome_1<-apply(matrix(treatment_outcome,ncol=2,nrow=4),1,as.numeric)
dimnames(treatment_outcome_1) <- list(Group = c("WT", "LD"),treatment_outcome = c("Complete Remission/Response", "Partial Remission/Response", "Progressive Disease", "Stable Disease"))
chisq<- chisq.test(unlist(t(treatment_outcome_1)))


####tumor_type
result_tumor_type<-clinical_stage %>% 
  group_by(Group) %>% count(tumor_type)
result_tumor_type<-result_tumor_type[-which(is.na(result_tumor_type$tumor_type)),]
result_tumor_type$tumor_type1<-factor(result_tumor_type$tumor_type, levels=c("Biochemical evidence of disease", "Regional lymph node", "Locoregional tumor", "New primary tumor", "Metastasis", "Recurrence"))
result_tumor_type$n <- ave(result_tumor_type$n, result_tumor_type$tumor_type, FUN=function(x) x/max(x)) 

tumor_type_plot<-ggplot(result_tumor_type, aes(fill=as.factor(tumor_type1), y=n, x=Group)) + 
  geom_bar(position="fill", stat="identity")+scale_fill_brewer(palette = "GnBu")+ labs(fill = "Factor_Factor_Factor")+ ggtitle("tumor_type") 

pdf(file="TCGA_MRP_tumor_type_plot_nor.pdf", width=4, height=2.5, compress=TRUE, useDingbats=F)
tumor_type_plot
dev.off()



clinical_tumor_type_MRP_group_n<-read.delim("/Users/ching-wenchang/Documents/Sur/clinical_tumor_type_MRP_group_n.txt",header = T)
tumor_type <- as.table(as.matrix(clinical_tumor_type_MRP_group_n))
rownames(tumor_type)<-tumor_type[ ,1]
tumor_type<-tumor_type[ ,-1]
tumor_type<-tumor_type[,-c(2,3)]
tumor_type_1<-apply(matrix(tumor_type,ncol=2,nrow=6),1,as.numeric)
dimnames(tumor_type_1) <- list(Group = c("WT", "LD"),tumor_type = c("Biochemical evidence of disease", "Locoregional tumor", "Metastasis", "New primary tumor", "Recurrence", "Regional lymph node"))
chisq<- chisq.test(unlist(t(tumor_type_1)))

