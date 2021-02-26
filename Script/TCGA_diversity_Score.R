RNA_data <- read.delim("/Users/ching-wenchang/Documents/Pan_cancer_mRNA/All_TCGA_mRNA_group.txt",header = T,row.names = NULL)
RNA_data[1:10,1:10]

setwd("/Users/ching-wenchang/Documents/TCGA_diversity")

library('DescTools')
library('entropy')
RNA_data$X<-as.factor(RNA_data$X)
RNA_data$row.names<-as.factor(RNA_data$row.names)
RNA_data$PATIENT_ID<-as.factor(RNA_data$PATIENT_ID)
RNA_data$Group<-as.factor(RNA_data$Group)
RNA_data$Cancer_Type<-as.factor(RNA_data$Cancer_Type)
RNA_data$group_new<-as.factor(RNA_data$group_new)
RNA_data$Hugo_Symbol<-as.factor(RNA_data$Hugo_Symbol)
RNA_data_del<-(RNA_data)[,-c(1:8)]
RNA_data_del$X <- NULL
RNA_data_del[1:10,1:10]
RNA_data_log<-log(RNA_data_del)
RNA_data_log[1:10,1:10]
system.time(RNA_data_log_dat <- do.call(data.frame,lapply(RNA_data_log, function(x) replace(x, is.infinite(x),NA))))
system.time(is.na(dat) <- sapply(dat, is.infinite))
RNA_data_log1<-do.call(data.frame,lapply(RNA_data_log, function(x) replace(x, is.infinite(x),NA)))
RNA_data_log1[is.na(RNA_data_log1)] <- 0
RNA_data_log1[1:10,1:10]
Pan_cancer_diversity1<-cbind(RNA_data[,c(1:8)], RNA_data_log1)
Pan_cancer_diversity1[1:10,1:10]
Entropy_value1 <-c()
for (patient in 1:nrow(Pan_cancer_diversity1))
  
{
  #patient <- 1
  tem <- Entropy(as.numeric(Pan_cancer_diversity1[patient,9:ncol(Pan_cancer_diversity1),]) , y = NULL, base = 2)
  Entropy_value1 <- c(Entropy_value1,tem)
}
Pan_cancer_diversity1_diversity<- cbind(diversity = Entropy_value1,Pan_cancer_diversity1)
Pan_cancer_diversity1_diversity[1:10,1:10]
summary(Pan_cancer_diversity1_diversity$diversity)
RNA_data_diversity<-Pan_cancer_diversity1_diversity
doit <- function(x) {(x - min(x, na.rm=TRUE))/(max(x,na.rm=TRUE) - 
                                                 min(x, na.rm=TRUE))} 
RNA_data_diversity$diversity_pct_Per = doit(RNA_data_diversity$diversity)
RNA_data_diversity[1:10,1:10]
write.table(RNA_data_diversity,"Pan_cancer_RNA_data_diversity.txt",sep="\t",quote = F)
RNA_data_diversity_all<-cbind.data.frame(Tumor_Sample_Barcode=(as.character( RNA_data_diversity$PATIENT_ID)), diversity=(as.numeric(RNA_data_diversity$diversity)) )
RNA_data_diversity_all<-cbind.data.frame(RNA_data_diversity_all, diversity_pct_Per=(RNA_data_diversity$diversity_pct_Per))
write.table(RNA_data_diversity_all,"Pan_cancer_RNA_data_diversity2.txt",sep="\t",quote = F)

#########Cancer_Type
pan_cancer_MRP_per_Group28<- read.delim("/Users/ching-wenchang/Documents/Pan_cancer_CNV_mRNA/pan_cancer_MRP_per_Group28_LDHD_1.txt",header = T)
diversity_pct_Per_mean<-aggregate(RNA_data_diversity$diversity_pct_Per, list(RNA_data_diversity$Cancer_Type), mean)
names(diversity_pct_Per_mean)[names(diversity_pct_Per_mean) == "Group.1"] <- "Type"
names(diversity_pct_Per_mean)[names(diversity_pct_Per_mean) == "x"] <- "diversity_Per"
diversity_pct_Per_MRP<-merge(diversity_pct_Per_mean, pan_cancer_MRP_per_Group28, by= "Type")
diversity_pct_Per_MRP$Type1 <-factor(diversity_pct_Per_MRP$Type, levels=c("THCA",
                                                                                  "AML",
                                                                                  "UVM",
                                                                                  "THYM",
                                                                                  "PRAD",
                                                                                  "LGG",
                                                                                  "KIRP",
                                                                                  "KIRC",
                                                                                  "UCEC",
                                                                                  "MESO",
                                                                                  "GBM",
                                                                                  "CESC",
                                                                                  "COAD",
                                                                                  "PAAD",
                                                                                  "LIHC",
                                                                                  "ACC",
                                                                                  "HNSCC",
                                                                                  "READ",
                                                                                  "STAD",
                                                                                  "BRCA",
                                                                                  "SKCM",
                                                                                  "SARC",
                                                                                  "LUAD",
                                                                                  "BLCA",
                                                                                  "ESCA",
                                                                                  "LUSC",
                                                                                  "ESAD",
                                                                                  "OV"))



library(ggplot2)
gg <- ggplot(diversity_pct_Per_MRP, aes(x=diversity_Per, y=Per)) + 
  geom_point(aes(col=Type1), size=3) +  # Set color to vary based on state categories.
  geom_smooth(method="lm", col="firebrick", size=2) +  
  labs(title="MRP_deletion_Per_diversity",  y="Frequences of MRP Deletion", x="Diversity_Score", caption="R=0.43, p=0.024")
plot(gg)

pdf(file="Pan_cancer_MRP_deletion_Per_diversity.pdf",width=4.5, height=3,compress=TRUE,useDingbats=F)
gg
dev.off()

write.table(diversity_pct_Per_MRP, "diversity_pct_Per_MRP.txt", quote = F, sep = "\t")
###############


Pan_cancer_cor_Per <- read.delim("/Users/ching-wenchang/Documents/Pan_cancer_CNV_mRNA/pan_cancer_MRP_per_Group28_LDHD.txt")

library("ggpubr")
g <-ggscatter(Pan_cancer_cor_Per, x = "cor", y = "Per", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "cor", ylab = "Del_Per")

pdf(file="pan_cancer_cor__Del_Per_28.pdf",width=2.5, height=2,compress=TRUE,useDingbats=F)
g
dev.off()

Pan_cancer_cor_Per$Type1 <-factor(Pan_cancer_cor_Per$Type, levels=c("THCA",
                                                                          "AML",
                                                                          "UVM",
                                                                          "THYM",
                                                                          "PRAD",
                                                                          "LGG",
                                                                          "KIRP",
                                                                          "KIRC",
                                                                          "UCEC",
                                                                          "MESO",
                                                                          "GBM",
                                                                          "CESC",
                                                                          "COAD",
                                                                          "PAAD",
                                                                          "LIHC",
                                                                          "ACC",
                                                                          "HNSCC",
                                                                          "READ",
                                                                          "STAD",
                                                                          "BRCA",
                                                                          "SKCM",
                                                                          "SARC",
                                                                          "LUAD",
                                                                          "BLCA",
                                                                          "ESCA",
                                                                          "LUSC",
                                                                          "ESAD",
                                                                          "OV"))

library(ggplot2)
gg <- ggplot(Pan_cancer_cor_Per, aes(x=cor, y=Per)) + 
  geom_point(aes(col=Type1), size=3) +  # Set color to vary based on state categories.
  geom_smooth(method="lm", col="firebrick", size=2) +  
  labs(title="MRP_deletion_Per_diversity",  y="Frequences of MRP Deletion", x="CNV-mRNA Corr", caption="R=0.88, p=4.1e-10")
plot(gg)

pdf(file="pan_cancer_cor__Del_Per_28_1.pdf",width=4.5, height=3,compress=TRUE,useDingbats=F)
gg
dev.off()

###############
Pan_cancer_HR_Per <- read.delim("/Users/ching-wenchang/Documents/TCGA_diversity/pan_cancer_HR_All_27.txt")
diversity_pct_Per_MRP_HR<-merge(diversity_pct_Per_MRP, Pan_cancer_HR_Per, by = "Type")
write.table(diversity_pct_Per_MRP_HR, "diversity_pct_Per_MRP_HR.txt", quote = F, sep = "\t")
diversity_pct_Per_MRP_HR <- read.delim("/Users/ching-wenchang/Documents/TCGA_diversity/diversity_pct_Per_MRP_HR1.txt")
library(ggplot2)
diversity_pct_Per_MRP_HR_Plot<- ggplot(diversity_pct_Per_MRP_HR, aes(x=diversity_Per, y=HR1)) + 
  geom_point(aes(col=Type1), size=3) +  # Set color to vary based on state categories.
  geom_smooth(method="lm", col="firebrick", size=2) +  
  labs(title="MRP_deletion_Per_diversity",  y="HR1", x="Diversity_Score", caption="R=0.43, p=0.024")
plot(diversity_pct_Per_MRP_HR_Plot)
ggscatter(diversity_pct_Per_MRP_HR, x = "diversity_Per", y = "HR1", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "diversity_Per", ylab = "Per")


diversity_pct_Per_MRP_HR$diversity_Type<-as.factor(diversity_pct_Per_MRP_HR$diversity_Type)
diversity_pct_Per_MRP_HR$Type<-as.factor(diversity_pct_Per_MRP_HR$Type)


diversity_pct_Per_MRP_HR$diversity_Type1 <-factor(diversity_pct_Per_MRP_HR$diversity_Type, levels=c("Low", "High"))


diversity_pct_Per_MRP_HR_Plot1<-ggplot(diversity_pct_Per_MRP_HR, aes(x= diversity_Type1, y=HR1, fill=diversity_Type1) ) + 
  geom_violin(trim=FALSE, size=0.5)+
  geom_boxplot(width=0.1, fill="white", size=0.5)+
  labs(title=paste("RNA_data_diversity_HR_Plot"),x="Diversity_Group", y = "HR")+stat_compare_means(method = "anova")+
  theme_classic()+ stat_compare_means(label="p.signif", method = "t.test")+ scale_fill_manual(values=c("#F067A6", "#2b8cbe"))


pdf(file="Pan_cancer_HR_Per_diversity.pdf",width=3, height=2,compress=TRUE,useDingbats=F)
diversity_pct_Per_MRP_HR_Plot1
dev.off()
