setwd("/Users/changc11/Documents/mutation_sing/mut_files")

Pan_cancer_TP53_mutation_group_all<- read.delim("/Users/changc11/Documents/Pan_cancer_mutation2/Mutation_TP53/Pan_cancer_TP53_mutation_group_all.txt",header = T)
Pan_cancer_WT_Sing_Patient_100<- read.delim("/Users/changc11/Documents/mutation_sing/mut_files/Pan_cancer_WT_Sing_Patient_100.txt",header = T)
Pan_cancer_MD.LD_Sing_Patient_100 <- read.delim("~/Documents/mutation_sing/mut_files/Pan_cancer_MD.LD_Sing_Patient_100.txt", header=T, row.names=1)


library(ggplot2)
library(dplyr)
library(ggpubr)
Pan_cancer_WT_Sing_Patient_100_P<-rownames(Pan_cancer_WT_Sing_Patient_100)
Pan_cancer_WT_Sing_Patient_100_P<-gsub(".weights","", Pan_cancer_WT_Sing_Patient_100_P)
rownames(Pan_cancer_WT_Sing_Patient_100)<-Pan_cancer_WT_Sing_Patient_100_P
Pan_cancer_WT_Sing_Patient_100$Tumor_Sample_Barcode<-Pan_cancer_WT_Sing_Patient_100_P
tp53_mutation_Pan_cancer_WT_Sing_Patient_100<-merge(Pan_cancer_TP53_mutation_group_all, Pan_cancer_WT_Sing_Patient_100, by="Tumor_Sample_Barcode")
Pan_cancer_MD.LD_Sing_Patient_100_P<-rownames(Pan_cancer_MD.LD_Sing_Patient_100)
Pan_cancer_MD.LD_Sing_Patient_100_P<-gsub(".weights","", Pan_cancer_MD.LD_Sing_Patient_100_P)
rownames(Pan_cancer_MD.LD_Sing_Patient_100)<-Pan_cancer_MD.LD_Sing_Patient_100_P
Pan_cancer_MD.LD_Sing_Patient_100$Tumor_Sample_Barcode<-Pan_cancer_MD.LD_Sing_Patient_100_P
tp53_mutation_Pan_cancer_MD.LD_Sing_Patient_100<-merge(Pan_cancer_TP53_mutation_group_all, Pan_cancer_MD.LD_Sing_Patient_100, by="Tumor_Sample_Barcode")
Pan_cancer_Sing<-rbind(tp53_mutation_Pan_cancer_WT_Sing_Patient_100, tp53_mutation_Pan_cancer_MD.LD_Sing_Patient_100)
Pan_cancer_Sing_all<-cbind(group_new=Pan_cancer_Sing$group_new, Pan_cancer_Sing[,10:40])
Pan_cancer_Sing_all$group_new
Pan_cancer_Sing_all$group_new<-as.factor(Pan_cancer_Sing_all$group_new)
#write.table(Pan_cancer_Sing,file="Pan_cancer_Sing.txt",quote = F,sep = "\t")

library(dplyr)
library('tidyverse')
library('broom')
Pan_cancer_Sing_all_sum<-Pan_cancer_Sing_all %>%
  group_by(group_new) %>%
  summarise_all(tibble::lst(mean))
#write.table(Pan_cancer_Sing_all_sum,file="Pan_cancer_Sing_all_sum.txt",quote = F,sep = "\t")


all.colors1 <-c("#67001f","#980043", "#ce1256", "#e7298a", "#df65b0", "#c994c7", "#d4b9da", "#e7e1ef",
                "#ffffcc", "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#4292c6", "#2171b5", "#08519c", "#08306b",
                "#bd0026", "#e31a1c", "#fc4e2a", "#fd8d3c", "#feb24c", "#fed976", "#ffeda0", "#ccece6", "#99d8c9",
                "#66c2a4", "#41ae76", "#238b45", "#006d2c", "#00441b")


all.sigs1 <-c("Signature.7","Signature.11", "Signature.4", "Signature.13", "Signature.24", "Signature.29", "Signature.2", "Signature.22",
              "Signature.3", "Signature.9", "Signature.1", "Signature.26", "Signature.21", "Signature.15", "Signature.20", "Signature.6", "Signature.10",
              "Signature.25", "Signature.8", "Signature.19", "Signature.17", "Signature.unknown", "Signature.18", "Signature.16", "Signature.5", "Signature.30",
              "Signature.27", "Signature.12", "Signature.23", "Signature.28", "Signature.14")

Pan_cancer_Sing_all_sum_t_1 <- read.delim("~/Documents/mutation_sing/mut_files/Pan_cancer_Sing_all_sum_t.txt")
library(ggplot2)


Pan_cancer_Sing_all_sum_t_1$group_new <- factor(Pan_cancer_Sing_all_sum_t_1$group_new, levels = all.sigs1)
levels(Pan_cancer_Sing_all_sum_t_1$group_new)
Pan_cancer_Sing_all_sum_t_1$Group<-factor(Pan_cancer_Sing_all_sum_t_1$Group, levels = c("non_deletion", "P53", "MRP_del_LD", "P53_MRP_del"))
# Barplot
bp<- ggplot(Pan_cancer_Sing_all_sum_t_1, aes(x=Group, y=Signature, fill=group_new))+
  geom_bar(width = 1, stat = "identity")
bp+ scale_fill_manual(values=all.colors1)
pdf(file="Pan_cancer_All_plotSignatures1.pdf",width=6, height=4,compress=TRUE,useDingbats=F)
bp+ scale_fill_manual(values=all.colors1)
dev.off()


########################################
library(plyr)
library(dplyr)
m1 = ddply(Pan_cancer_Sing_all_sum_t_1, .(group_new), summarize, ratio=Signature/sum(Signature))
m2 = Pan_cancer_Sing_all_sum_t_1[order(Pan_cancer_Sing_all_sum_t_1$group_new),]
# combine them:
mydf = data.frame(m2,ratio=m1$ratio)
mydf = ddply(mydf, .(group_new), transform, position = cumsum(Signature) - 0.5*Signature) 
#write.table(mydf,paste("Pan_cancer_group_Sing",patient_name,sep = "",".txt"),quote = F,sep = "\t")



mydf_Group<-mydf %>%
  group_split(Group)
non_deletion <-mydf_Group[[1]]
non_deletion<-data.frame(non_deletion)
names(non_deletion)[names(non_deletion) == "ratio"] <- "non_deletion"
P53 <-mydf_Group[[2]]
P53<-data.frame(P53)
names(P53)[names(P53) == "ratio"] <- "P53"
MRP_del <-mydf_Group[[3]]
MRP_del<-data.frame(MRP_del)
names(MRP_del)[names(MRP_del) == "ratio"] <- "MRP_del"
P53_MRP_del <-mydf_Group[[4]]
P53_MRP_del<-data.frame(P53_MRP_del)
names(P53_MRP_del)[names(P53_MRP_del) == "ratio"] <- "P53_MRP_del"
mydf_Group_4<-data.frame(non_deletion, P53, MRP_del, P53_MRP_del)
mydf_Group_4_P53_WT = ddply(mydf_Group_4, .(group_new), summarize, P53_WT=P53/non_deletion) 
mydf_Group_4_P53_WT$P53_WT_ratio_log<-log2(mydf_Group_4_P53_WT$P53_WT)
mydf_Group_4_MRP_del_WT = ddply(mydf_Group_4, .(group_new), summarize, MRP_del_WT=MRP_del/non_deletion) 
mydf_Group_4_MRP_del_WT$MRP_del_WT_ratio_log<-log2(mydf_Group_4_MRP_del_WT$MRP_del_WT)
mydf_Group_4_P53_MRP_del_WT = ddply(mydf_Group_4, .(group_new), summarize, P53_MRP_del_WT=P53_MRP_del/non_deletion) 
mydf_Group_4_P53_MRP_del_WT$P53_MRP_del_WT_ratio_log<-log2(mydf_Group_4_P53_MRP_del_WT$P53_MRP_del_WT)
mydf_Group_4_P53_MRP_del_P53 = ddply(mydf_Group_4, .(group_new), summarize, P53_MRP_del_P53=P53_MRP_del/P53) 
mydf_Group_4_P53_MRP_del_P53$P53_MRP_del_P53_ratio_log<-log2(mydf_Group_4_P53_MRP_del_P53$P53_MRP_del_P53)
mydf_Group_4_foldchange<-merge(mydf_Group_4_P53_WT, mydf_Group_4_MRP_del_WT, by="group_new")
mydf_Group_4_foldchange<-merge(mydf_Group_4_foldchange, mydf_Group_4_P53_MRP_del_WT, by="group_new")
mydf_Group_4_foldchange<-merge(mydf_Group_4_foldchange, mydf_Group_4_P53_MRP_del_P53, by="group_new")
mydf_Group_4_foldchange$group_new <- factor(mydf_Group_4_foldchange$group_new, levels = all.sigs1)


pdf(file="Pan_cancer_All_plotSignatures_foldchange_P53_MRP_del_WT.pdf",width=5, height=3,compress=TRUE,useDingbats=F)
ggplot(mydf_Group_4_foldchange, aes(group_new, P53_MRP_del_WT_ratio_log, fill=ifelse(P53_MRP_del_WT_ratio_log>0,"P53_MRP_del","WT"))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("red", "blue"), name="WT or P53_MRP_del") +
  scale_y_continuous(labels = function(y) y )+ theme(axis.text.x = element_text(face="bold", color="black", size=8, angle=45))
dev.off()

pdf(file="Pan_cancer_All_plotSignatures_foldchange_P53_MRP_del_P53.pdf",width=5, height=3,compress=TRUE,useDingbats=F)
ggplot(mydf_Group_4_foldchange, aes(group_new, P53_MRP_del_P53_ratio_log, fill=ifelse(P53_MRP_del_P53_ratio_log>0,"P53_MRP_del","P53"))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("red", "blue"), name="P53 or P53_MRP_del") +
  scale_y_continuous(labels = function(y) y )+ theme(axis.text.x = element_text(face="bold", color="black", size=8, angle=45))
dev.off()

pdf(file="Pan_cancer_All_plotSignatures_foldchange_MRP_del_WT.pdf",width=5, height=3,compress=TRUE,useDingbats=F)
ggplot(mydf_Group_4_foldchange, aes(group_new, MRP_del_WT_ratio_log, fill=ifelse(MRP_del_WT_ratio_log>0,"MRP_del","WT"))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("red", "blue"), name="WT or MRP_del") +
  scale_y_continuous(labels = function(y)y )+ theme(axis.text.x = element_text(face="bold", color="black", size=8, angle=45))
dev.off()

pdf(file="Pan_cancer_All_plotSignatures_foldchange_P53_WT.pdf",width=5, height=3,compress=TRUE,useDingbats=F)
ggplot(mydf_Group_4_foldchange, aes(group_new, P53_WT_ratio_log, fill=ifelse(P53_WT_ratio_log>0,"P53","WT"))) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("red", "blue"), name="WT or P53") +
  scale_y_continuous(labels = function(y) y )+ theme(axis.text.x = element_text(face="bold", color="black", size=8, angle=45))
dev.off()



library(dplyr)
library(ggpubr)
library(ggplot2)
Pan_cancer_Sing_t<- read.delim("/Users/changc11/Documents/mutation_sing/mut_files/Pan_cancer_Sing_t.txt",header = T)
# Box plot facetted by "dose"
p <- ggboxplot(Pan_cancer_Sing_t, x = "group_new", y = "Signature",
               color = "group_new", palette = "jco",
               facet.by = "Signature_Type", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(method = "anova", label.y = 1.2)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = ".all.", data =Pan_cancer_Sing_t)+ theme(axis.text.x = element_text(face="bold", color="black", size=8, angle=45))           
Pan_cancer_Sing_p_Wilcoxon<-compare_means(Signature~group_new, data=Pan_cancer_Sing_t, group.by = "Signature_Type")
Pan_cancer_Sing_p<-compare_means(Signature~group_new, data=Pan_cancer_Sing_t, group.by = "Signature_Type",method = "t.test")
write.table(Pan_cancer_Sing_p,file="Pan_cancer_Sing_p.txt",quote = F,sep = "\t")
write.table(Pan_cancer_Sing_p_Wilcoxon,file="Pan_cancer_Sing_p_Wilcoxon.txt",quote = F,sep = "\t")



##############################
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
require(graphics)
library(RColorBrewer)

mut.sig_p<- read.delim("/Users/changc11/Documents/mutation_sing/mut.sig_p.txt",header = T)
mut.sig_p <- mut.sig_p[,-1]
breaks = c(0, 0.001, 0.05, 0.06, 1)
colors = c("red", "#fb6a4a", "#ffffb2", "#bdd7e7")
p2<-pheatmap(mut.sig_p, color = colors, 
             show_rownames=T,
             show_colnames =T, breaks=breaks, cluster_rows=F, cluster_cols=F)
pdf(file="Pan_cancer_mut_sig_heatmapcrp.pdf",width=3, height=6)
p2
dev.off()
