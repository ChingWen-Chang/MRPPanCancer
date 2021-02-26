
library(maftools)
laml.WT <- read.maf(maf = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_WT.maf', clinicalData = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_WT_clinical_MSKCC.tsv')
laml.MD_HD <- read.maf(maf = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_MD_HD.maf', clinicalData = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_MD_HD_clinical_MSKCC.tsv')
laml.LD <- read.maf(maf = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_LD.maf', clinicalData = '/Users/ching-wenchang/Documents/MSK-IM/MSKCC_mutation/dataset_LD_clinical_MSKCC.tsv')

plotmafSummary(maf = laml.WT, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = laml.MD_HD, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
plotmafSummary(maf = laml.LD, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#####################oncoplot plot
col = RColorBrewer::brewer.pal(n = 9, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del', 'Translation_Strat_Site')

#Color coding for MRP classification
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

#Color coding for Overall_Survival_Status
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


#########################Comparison plot between different deletion groups
pt.vs.rt <- mafCompare(m1 = laml.MD_HD, m2 = laml.WT, m1Name = 'Deletion', m2Name = 'Non_Deletion', minMut = 5)
print(pt.vs.rt)
write.table(pt.vs.rt$results,"Pan_cancer_mutation_compre_MDHD_MSKCC.txt",sep="\t",quote = F)
pdf(file="Pan_cancer_mutation_compre_MDHD_WT_MSKCC1.pdf",width=4, height=3,compress=TRUE,useDingbats=F)
forestPlot1(mafCompareRes = pt.vs.rt, pVal = 0.0000005, color = c('royalblue', 'maroon'), geneFontSize = 0.9)
dev.off()


#######################lollipop plot for TP53 in different MRP groups
#lollipop plot for TP53, which is one of the most frequent mutated gene in MRP deletion group.
pdf(file="Pan_cancer_mutation_compre_MDHD_WT_MSKCC_TP53.pdf",width=4, height=6,compress=TRUE,useDingbats=F)
lollipopPlot2(m1 = laml.MD_HD, m2 = laml.WT, gene = "TP53", AACol1 = "HGVSp", AACol2 = "HGVSp", m1_name = "Deletion", m2_name = "Non_Deletion")
dev.off()

#lollipop plot for TP53 in WT group
pdf(file="Pan_cancer_mutation_WT_MSKCC_TP53.pdf",width=4, height=1.6,compress=TRUE,useDingbats=F)
lollipopPlot(maf = laml.WT, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE, labelPos = 248)
dev.off()

#lollipop plot for TP53 in MDHD group
pdf(file="Pan_cancer_mutation_MDHD_MSKCC_TP53.pdf",width=4, height=3.8,compress=TRUE,useDingbats=F)
lollipopPlot(maf = laml.MD_HD, gene = 'TP53', AACol = 'HGVSp', showMutationRate = TRUE, labelPos = 248)
dev.off()


#############################Survival Analysis

MD_HD<-mafSurvival1(maf = laml.MD_HD, genes = 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = F)
WT<-mafSurvival1(maf = laml.WT, genes = 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = F)
LD<-mafSurvival1(maf = laml.LD, genes = 'TP53', time = 'days_to_last_followup', Status = 'Overall_Survival_Status', isTCGA = F)


library("survival")
library("survminer")

group_new<-MD_HD$Group
group_new<-gsub("WT","HMD",group_new)
group_new<-gsub("Mutant","P53-mut/HMD",group_new)
data <-cbind(MD_HD,group_new)

group2<-LD$Group
group_new<-gsub("WT","LD",group2)
group_new<-gsub("Mutant","P53-mut/LD",group_new)
data2 <-cbind(LD,group_new)

group1<-WT$Group
group_new<-gsub("WT","WT",group1)
group_new<-gsub("Mutant","P53-mut",group_new)
data1 <-cbind(WT,group_new)

patientdata2 <- rbind(data1, data, data2)
fit2 <- survfit(Surv(Time, Status) ~ group_new, data = patientdata2)
print(fit2)
OS<-ggsurvplot(fit2,#facet.by = "Cancer_Type",
               pval = TRUE, conf.int = F,
               risk.table = F, # Add risk table
               #risk.table.col = "strata", # Change risk table color by groups
               linetype = "strata", # Change line type by groups
               surv.median.line = "hv", # Specify median survival
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c("#f03b20", "#54278f", "#08519c" , "#dd1c77", "#beaed4" , "#41b6c4"),
               censor.size=2,font.main = c(10),
               font.x = c(10),
               font.y = c(10),
               font.tickslab = c(10),
               fontsize=10,break.time.by = 12,xlim = c(0, 48))
pdf(file="Pan_cancer_p53_del_Sur_MSKCC.pdf",width=4, height=4,compress=TRUE,useDingbats=F)
OS
dev.off()




