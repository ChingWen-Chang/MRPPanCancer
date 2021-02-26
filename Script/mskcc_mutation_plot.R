#Install from GitHub repository
BiocManager::install("PoisonAlien/maftools")



setwd("/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation")

data_mutations_mskcc <- read.delim("~/Documents/MSK-IM/msk_impact_2017/data_mutations_mskcc.txt", comment.char="#")
data_mutations_mskcc[1:10,1:30]
#MSKCC_MRP_sur_Group<- read.delim("/Users/ching-wenchang/Documents/MSK-IM/MSKCC_clinical_inf.txt")

MSKCC_MRP_sur_Group <- read.delim("/Users/ching-wenchang/Documents/MSK-IM/MSKCC_MRP_sur_Group2.txt")

#MSKCC_MRP_sur_Group <- read.delim("/Users/ching-wenchang/Documents/MSK-IM/MSKCC_MRP_sur_Group.txt")
#MSKCC_MRP_sur_Group$Group1<-gsub("1,5", "0,1",MSKCC_MRP_sur_Group$Group)

#Pan_cancer_TP53_mutation_group_all <- read.delim("/Users/ching-wenchang/Documents/Pan_cancer_mutation2/Mutation_TP53/Pan_cancer_TP53_mutation_group_all.txt")
#TCGA_MRP_sur_Group<-read.delim("/Users/ching-wenchang/Documents/Sur-All/patientdata_status_20200524.txt",header = T)
#Pan_cancer_TP53_mutation_group_all_per<-merge(Pan_cancer_TP53_mutation_group_all, TCGA_MRP_sur_Group, by="Tumor_Sample_Barcode")
#write.table(Pan_cancer_TP53_mutation_group_all_per,"Pan_cancer_TP53_mutation_group_all_per.txt",sep="\t",quote = F)

MSKCC_mutation_Group<-merge(MSKCC_MRP_sur_Group, data_mutations_mskcc, by="Tumor_Sample_Barcode")
MSKCC_mutation_Group[1:10,1:30]

dataset<-MSKCC_mutation_Group

dataset_WT<-dataset[which(dataset$Group1=="WT"),]
write.table(dataset_WT,"dataset_WT.txt",sep="\t",quote = F)
dataset_WT_clinical<-dataset_WT[1:102]
write.table(dataset_WT_clinical,"dataset_WT_clinical.txt",sep="\t",quote = F)

dataset_MD_HD<- dataset[which(dataset$Group1=="MDH"),]
write.table(dataset_MD_HD,"dataset_MD_HD.txt",sep="\t",quote = F)
dataset_MD_HD_clinical<-dataset_MD_HD[1:102]
write.table(dataset_MD_HD_clinical,"dataset_MD_HD_clinical.txt",sep="\t",quote = F)

dataset_LD <- dataset[which(dataset$Group1=="LD"),]
write.table(dataset_LD,"dataset_LD.txt",sep="\t",quote = F)
dataset_LD_clinical<-dataset_LD[1:102]
write.table(dataset_LD_clinical,"dataset_LD_clinical.txt",sep="\t",quote = F)


dataset_WT1<-dataset[which(dataset$Group=="[0,1)"),]
write.table(dataset_WT1,"dataset_WT1.txt",sep="\t",quote = F)
dataset_WT1_clinical<-dataset_WT1[1:102]
write.table(dataset_WT1_clinical,"dataset_WT1_clinical.txt",sep="\t",quote = F)

library(maftools)
laml.WT <- read.maf(maf = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_WT.maf', clinicalData = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_WT_clinical_MSKCC.tsv')
laml.MD_HD <- read.maf(maf = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_MD_HD.maf', clinicalData = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_MD_HD_clinical_MSKCC.tsv')
laml.LD <- read.maf(maf = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_LD.maf', clinicalData = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_LD_clinical_MSKCC.tsv')

laml.WT1 <- read.maf(maf = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_WT1.maf', clinicalData = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_WT1_clinical_MSKCC.tsv')

plotmafSummary(maf = laml.WT, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = laml.MD_HD, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = laml.LD, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = laml.WT1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

col = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del', 'Translation_Strat_Site')

#Color coding for FAB classification
fabcolors = RColorBrewer::brewer.pal(n = 3,name = 'Spectral')
names(fabcolors) = c("WT", "LD", "MHD")
fabcolors = list(FAB_classification = fabcolors)
#Color coding for Cancer_Type
Cancer_Type_colors = c("#0868ac", "#3182bd", "#2b8cbe", "#43a2ca", "#6baed6",  "#74a9cf", "#bdd7e7", "#b3cde3", "#bdc9e1", "#bdd7e7",
                       "#d7b5d8", "#decbe4", "#decbe4", "#f6eff7",
                       "#e78ac3",
                       "#54278f", "#756bb1", "#9e9ac8", "#cbc9e2",
                       "#c51b8a", "#dd1c77", "#df65b0", "#f1b6da", "#f4cae4",
                       "#31a354", "#edf8e9",
                       "#fd8d3c", "#feedde")
names(Cancer_Type_colors) = c("BRCA", "OV", "UCEC", "LUAD", "PRAD", "COAD", "STAD", "PAAD", "READ", "ESAD",
                              "HNSCC", "LUSC", "CESC", "ESCA",
                              "SARC",
                              "GBM", "LGG", "SKCM", "UVM",
                              "THCA", "BLCA", "LIHC", "ACC", "MESO",
                              "AML", "THYM",
                              "KIRC", "KIRP")

OS_colors = c("#d7b5d8", "#980043")
names(OS_colors) = c("0", "1")
anno_cols = list(Overall_Survival_Status = OS_colors, Cancer_Type = Cancer_Type_colors)
print(anno_cols)

pdf(file="Pan_cancer_mutation_WT.MSKCC1.pdf",width=7, height=8,compress=TRUE,useDingbats=F)
oncoplot(maf = laml.WT,
         clinicalFeatures = c('Cancer_Type', 'Overall_Survival_Status'),
         sortByAnnotation = T,
         #mutsig = '/Users/ching-wenchang/Documents/Pan_cancer_mutation2/WT_OUTPUT.sig_genes.txt',
         #mutsigQval = 0.00000001,
         colors = col,
         annotationColor = anno_cols, top = 20)
dev.off()

pdf(file="Pan_cancer_mutation_MD_HD.MSKCC.pdf",width=7, height=8,compress=TRUE,useDingbats=F)
oncoplot(maf = laml.MD_HD,
         clinicalFeatures = c('Cancer_Type', 'Overall_Survival_Status'),
         sortByAnnotation = T,
         #mutsig = '/Users/ching-wenchang/Documents/Pan_cancer_mutation2/MD_HD_OUTPUT.sig_genes.txt',
         #mutsigQval = 0.0000001,
         colors = col,
         annotationColor = anno_cols, top = 20)
dev.off()

pdf(file="Pan_cancer_mutation_LD.MSKCC.pdf",width=7, height=8,compress=TRUE,useDingbats=F)
oncoplot(maf = laml.LD,
         clinicalFeatures = c('Cancer_Type', 'Overall_Survival_Status'),
         sortByAnnotation = T,
         #mutsig = '/Users/ching-wenchang/Documents/Pan_cancer_mutation2/LD_OUTPUT.sig_genes.txt',
         #mutsigQval = 0.0000001,
         colors = col,
         annotationColor = anno_cols, top = 20)
dev.off()


#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(maf = laml.WT, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE)




pt.vs.rt <- mafCompare(m1 = laml.MD_HD, m2 = laml.WT, m1Name = 'Deletion', m2Name = 'Non_Deletion', minMut = 5)
print(pt.vs.rt)
write.table(pt.vs.rt$results,"Pan_cancer_mutation_compre_MDHD_MSKCC.txt",sep="\t",quote = F)

pt.vs.rt <- mafCompare(m1 = laml.LD, m2 = laml.WT, m1Name = 'Deletion', m2Name = 'Non_Deletion', minMut = 5)
print(pt.vs.rt)
write.table(pt.vs.rt$results,"Pan_cancer_mutation_compre_LD_MSKCC.txt",sep="\t",quote = F)

pdf(file="Pan_cancer_mutation_compre_MDHD_WT_MSKCC1.pdf",width=4, height=3,compress=TRUE,useDingbats=F)
forestPlot1(mafCompareRes = pt.vs.rt, pVal = 0.0000005, color = c('royalblue', 'maroon'), geneFontSize = 0.9)
dev.off()

pdf(file="Pan_cancer_mutation_compre_MDHD_WT_MSKCC_per.pdf",width=4, height=3,compress=TRUE,useDingbats=F)
coBarplot(m1 = laml.MD_HD, m2 = laml.WT, m1Name = "Deletion", m2Name = "Non_Deletion")
dev.off()

pdf(file="Pan_cancer_mutation_compre_MDHD_WT_MSKCC_TP53.pdf",width=4, height=6,compress=TRUE,useDingbats=F)
lollipopPlot2(m1 = laml.MD_HD, m2 = laml.WT, gene = "TP53", AACol1 = "HGVSp", AACol2 = "HGVSp", m1_name = "Deletion", m2_name = "Non_Deletion")
dev.off()


pdf(file="Pan_cancer_mutation_WT_MSKCC_TP53.pdf",width=4, height=1.6,compress=TRUE,useDingbats=F)
#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
WT_TP53<-lollipopPlot(maf = laml.WT, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE, labelPos = 248)
dev.off()

mutation_WT_TP53<-lollipopPlot_inf(maf = laml.WT, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE, labelPos = 248)
write.table(mutation_WT_TP53,"Pan_cancer_mutation_WT_TP53_MSKCC.txt",sep="\t",quote = F)

mutation_WT_TP53_sum<-mutation_WT_TP53 %>%
  group_by(pos) %>%
  summarise(sum=sum(count))
write.table(mutation_WT_TP53_sum,"Pan_cancer_mutation_WT_TP53_sum_MSKCC.txt",sep="\t",quote = F)
TP53_per<-mutation_WT_TP53_sum$sum/1409*100
mutation_WT_TP53_sum_per<-cbind(mutation_WT_TP53_sum,TP53_per)

pdf(file="Pan_cancer_mutation_MDHD_MSKCC_TP53.pdf",width=4, height=3.8,compress=TRUE,useDingbats=F)
#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
lollipopPlot(maf = laml.MD_HD, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE, labelPos = 248)
dev.off()

mutation_MDLD_TP53<-lollipopPlot_inf(maf = laml.MD_HD, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE, labelPos = 248)
write.table(mutation_MDLD_TP53,"Pan_cancer_mutation_MDLD_TP53_MSKCC.txt",sep="\t",quote = F)
mutation_MDLD_TP53_sum<-mutation_MDLD_TP53 %>%
  group_by(pos) %>%
  summarise(sum=sum(count))
write.table(mutation_MDLD_TP53_sum,"Pan_cancer_mutation_MDLD_TP53_sum_MSKCC.txt",sep="\t",quote = F)
MDLD_TP53_per<-mutation_MDLD_TP53_sum$sum/3990*100
mutation_MDLD_TP53_sum_per<-cbind(mutation_MDLD_TP53_sum,MDLD_TP53_per)

MSKCC_mutation_group_Tp53_sum<-merge(mutation_WT_TP53_sum_per, mutation_MDLD_TP53_sum_per, by="pos", all=TRUE)

MSKCC_mutation_group_Tp53_sum[is.na(MSKCC_mutation_group_Tp53_sum)] = 0
MSKCC_mutation_group_Tp53_foldchange<-MSKCC_mutation_group_Tp53_sum$MDLD_TP53_per/MSKCC_mutation_group_Tp53_sum$TP53_per
MSKCC_mutation_group_Tp53_sum_foldchange<-cbind(MSKCC_mutation_group_Tp53_sum, MSKCC_mutation_group_Tp53_foldchange)
write.table(MSKCC_mutation_group_Tp53_sum_foldchange,"MSKCC_mutation_group_Tp53_sum_foldchange.txt",sep="\t",quote = F)

#Requires BSgenome object
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm.WT = trinucleotideMatrix(maf = laml.WT, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = laml.tnm.WT, maf = laml.WT, pVal = 0.2)
library('NMF')

laml.sig.WT = extractSignatures(mat = laml.tnm.WT, n = 3)
laml.og30.cosm.WT = compareSignatures(nmfRes = laml.sig.WT, sig_db = "legacy")
laml.v3.cosm.WT = compareSignatures(nmfRes = laml.sig.WT, sig_db = "SBS")


library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm.WT$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

plotSignatures(nmfRes = laml.sig.WT, title_size = 0.8, sig_db = "SBS")
plotSignatures(nmfRes = laml.sig.WT, title_size = 0.8, sig_db = "legacy")




laml.tnm.MD_HD = trinucleotideMatrix(maf = laml.MD_HD, prefix = 'chr', add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
plotApobecDiff(tnm = laml.tnm.MD_HD, maf = laml.MD_HD, pVal = 0.2)
library('NMF')

laml.sig.MD_HD = extractSignatures(mat = laml.tnm.MD_HD, n = 4)
laml.og30.cosm.MD_HD = compareSignatures(nmfRes = laml.sig.MD_HD, sig_db = "legacy")
laml.v3.cosm.MD_HD = compareSignatures(nmfRes = laml.sig.MD_HD, sig_db = "SBS")


library('pheatmap')
pheatmap::pheatmap(mat = laml.og30.cosm.MD_HD$cosine_similarities, cluster_rows = FALSE, main = "cosine similarity against validated signatures")

plotSignatures(nmfRes = laml.sig.MD_HD, title_size = 0.8, sig_db = "SBS")
plotSignatures(nmfRes = laml.sig.MD_HD, title_size = 0.8, sig_db = "legacy")





pdf(file="Pan_cancer_TP53_MD.LD_MSKCC.pdf",width=3, height=4,compress=TRUE,useDingbats=F)
MD_HD<-mafSurvival(maf = laml.MD_HD, genes = 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = F)
dev.off()

pdf(file="Pan_cancer_TP53_WT_MSKCC.pdf",width=3, height=4,compress=TRUE,useDingbats=F)
WT<-mafSurvival(maf = laml.WT, genes = 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = F)
dev.off()


pdf(file="Pan_cancer_TP53_SD_MSKCC.pdf",width=3, height=4,compress=TRUE,useDingbats=F)
LD<-mafSurvival(maf = laml.LD, genes = 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = F)
dev.off()


library("survival")
library("survminer")

group_new<-MD_HD$Group
group_new<-gsub("WT","MRP_del_LD",group_new)
group_new<-gsub("Mutant","P53_MRP_del",group_new)
data <-cbind(MD_HD,group_new)

group2<-LD$Group
group_new<-gsub("WT","MRP_del_SD",group2)
group_new<-gsub("Mutant","P53_MRP_del_SD",group_new)
data2 <-cbind(LD,group_new)

group1<-WT$Group
group_new<-gsub("WT","non_deletion",group1)
group_new<-gsub("Mutant","P53",group_new)
data1 <-cbind(WT,group_new)

patientdata <- rbind(data1, data)

write.table(patientdata2,"Pan_cancer_TP53_mutation_group_all_MSKCC.txt",sep="\t",quote = F)

patientdata2 <- rbind(data1, data, data2)

fit <- survfit(Surv(Time, Status) ~ group_new, data = patientdata)
print(fit)

fit2 <- survfit(Surv(Time, Status) ~ group_new, data = patientdata2)
print(fit2)

#handle_merge1<-order(patientdata[,2], decreasing=F)
#all_patient_merge_reorder <- patientdata[handle_merge1,]
#all_patient_merge_reorder[1:5,1:5]
#rdu_all_patient_merge_reorder <-all_patient_merge_reorder[-(1:28),]
ggsurvplot(fit, patientdata, facet.by = "Cancer_Type", 
           pval = TRUE, xlim = c(0, 120),palette = c("#f03b20", "#54278f", "#08519c", "#dd1c77"))


#fit.list <- list(
#ph.ecog = fit, sex = fit1
#)
#ggsurv <- ggsurvplot(fit.list, data = prad_data_CNV1_hetMRP_Sur4, censor = FALSE,
#combine = TRUE, keep.data = TRUE, 
#palette = "Dark2", legend = "right",pval = TRUE)
ggsurvplot(fit,
           pval = TRUE, conf.int = F,
           risk.table = TRUE, # Add risk table
           #risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#f03b20", "#54278f", "#08519c", "#dd1c77"),
           censor.size=2,font.main = c(10),
           font.x = c(10),
           font.y = c(10),
           font.tickslab = c(10),
           fontsize=10)


#palette = "Dark2", legend = "right",pval = TRUE)

DF<-ggsurvplot(fit2,#facet.by = "Cancer_Type",
               pval = TRUE, conf.int = F,
               risk.table = F, # Add risk table
               #risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c("#f03b20", "#54278f", "#41b6c4", "#08519c", "#dd1c77", "#beaed4"),
               censor.size=2,font.main = c(10),
               font.x = c(10),
               font.y = c(10),
               font.tickslab = c(10),
               fontsize=10,break.time.by = 12,xlim = c(0, 48))
pdf(file="Pan_cancer_p53_del_Sur_MSKCC.pdf",width=4, height=4,compress=TRUE,useDingbats=F)
DF
dev.off()





library(pheatmap)

MSKCC_PanCancer_compare_3<- read.delim("~/Documents/MSK-IM/MSKCC_mutation/MSKCC_PanCancer_compare_3.txt")
row.names(MSKCC_PanCancer_compare_3)<-MSKCC_PanCancer_compare_3[,1]
MSKCC_PanCancer_compare_3 <- MSKCC_PanCancer_compare_3[,-1]
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
require(graphics)
library(RColorBrewer)
breaks = c(0, 0.001, 0.05, 0.06, 1)

color = colorRampPalette(c("white", "#f1b6da","#f1b6da", "#df65b0","#c2a5cf", "#8856a7","#8856a7","#8856a7"))(30)
p2<-pheatmap(MSKCC_PanCancer_compare_3, color = color,
             show_rownames=T,clustering_method = "median",
             show_colnames =T,  cluster_rows=F, cluster_cols=T)


pdf(file="Pan_cancer_muation_per1.pdf",width=3, height=6)
p2
dev.off()




##################
TCGA_MSKCC_mutation_group_Tp53<- read.delim("/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/TCGA_MSKCC_mutation_group_Tp53.txt")
TCGA_MSKCC_mutation_group_Tp53$pos<- as.factor(TCGA_MSKCC_mutation_group_Tp53$pos)

stat.test1 <- TCGA_MSKCC_mutation_group_Tp53 %>%
  group_by(pos) %>%
  pairwise_t_test(
    MDLD_TP53_per_1~Group, paired = TRUE, 
    p.adjust.method = "bonferroni"
  ) %>%
  select(-df, -statistic, -p) # Remove details
stat.test


stat.test <- TCGA_MSKCC_mutation_group_Tp53 %>%
  group_by(pos) %>%
  t_test(MDLD_TP53_per_1~Group, alternative = "two.sided", var.equal = T)

library(tidyverse)
library(rstatix)
library(ggpubr)
p <- ggplot(TCGA_MSKCC_mutation_group_Tp53, aes(x=pos, y=MDLD_TP53_per_1, fill=pos)) + 
  geom_boxplot()
p+facet_wrap(~Group)+test()
