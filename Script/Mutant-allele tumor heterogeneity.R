library(maftools)
library(mclust)
laml.MD.HD <- read.maf(maf = '/Users/changc11/Documents/Pan_cancer_mutation2/mutation/dataset_MD_LD.maf', clinicalData = '/Users/changc11/Documents/Pan_cancer_mutation2/mutation/dataset_MD_LD_clinical.tsv')
laml.WT <- read.maf(maf = '/Users/changc11/Documents/Pan_cancer_mutation2/mutation/dataset_WT.maf', clinicalData = '/Users/changc11/Documents/Pan_cancer_mutation2/mutation/dataset_WT_clinical.tsv')
laml.LD <- read.maf(maf = '/Users/changc11/Documents/Pan_cancer_mutation2/dataset_LD.maf', clinicalData = '/Users/changc11/Documents/Pan_cancer_mutation2/dataset_LD_clinical.tsv')


getFields(laml.WT)
laml.WT<-laml.WT
laml.WT@data$VAF <- laml.WT@data$t_alt_count/laml.WT@data$t_depth
head(laml.WT@data$Tumor_Sample_Barcode)
tcga.ab.2972.het = inferHeterogeneity(maf = laml.WT, tsb = 'TCGA-OR-A5J1', vafCol = 'VAF')
unique(tcga.ab.2972.het$clusterData$MATH)
plotClusters(clusters=tcga.ab.2972.het)

mutation_Heterogeneity <-c()
for (sample_name in na.omit(unique(laml.WT@data[["Tumor_Sample_Barcode"]])))
{
  #sample_name<-"TCGA-ZF-A9R4"
  tem <- try(inferHeterogeneity(maf = laml.WT, tsb = sample_name, vafCol = 'VAF'))
  if(!class(tem) %in% "try-error")
  {
    sample_num<-unique(tem$clusterData$MATH)
    names(sample_num) <-sample_name
    mutation_Heterogeneity<-c(mutation_Heterogeneity,sample_num)
  }
}
mutation_Heterogeneity_WT<-data.frame(mutation_Heterogeneity)
write.table(mutation_Heterogeneity_WT,file = "mutation_Heterogeneity_WT.txt",sep="\t",quote = F)

###########################
getFields(laml.LD)
laml.LD<-laml.LD
laml.LD@data$VAF <- laml.LD@data$t_alt_count/laml.LD@data$t_depth
head(laml.LD@data$Tumor_Sample_Barcode)
tcga.ab.2972.het = inferHeterogeneity(maf = laml.LD, tsb = 'TCGA-OR-A5J4', vafCol = 'VAF')
unique(tcga.ab.2972.het$clusterData$MATH)
plotClusters(clusters=tcga.ab.2972.het)
mutation_Heterogeneity <-c()
for (sample_name in na.omit(unique(laml.LD@data[["Tumor_Sample_Barcode"]])))
{
  #sample_name<-"TCGA-ZF-A9R4"
  tem <- try(inferHeterogeneity(maf = laml.LD, tsb = sample_name, vafCol = 'VAF'))
  if(!class(tem) %in% "try-error")
  {
    sample_num<-unique(tem$clusterData$MATH)
    names(sample_num) <-sample_name
    mutation_Heterogeneity<-c(mutation_Heterogeneity,sample_num)
    
  }
}

mutation_Heterogeneity_LD<-data.frame(mutation_Heterogeneity)
write.table(mutation_Heterogeneity_LD,file = "mutation_Heterogeneity_LD.txt",sep="\t",quote = F)

##########################
getFields(laml.MD.HD)
data<-laml.MD.HD
data@data$VAF <- data@data$t_alt_count/data@data$t_depth
head(data@data$Tumor_Sample_Barcode)
tcga.ab.2972.het = inferHeterogeneity(maf = data, tsb = 'TCGA-OR-A5J7', vafCol = 'VAF')
unique(tcga.ab.2972.het$clusterData$MATH)
plotClusters(clusters=tcga.ab.2972.het)
mutation_Heterogeneity <-c()
for (sample_name in na.omit(unique(data@data[["Tumor_Sample_Barcode"]])))
  
{
  #sample_name<-"TCGA-ZF-A9R4"
  tem <- try(inferHeterogeneity(maf = data, tsb = sample_name, vafCol = 'VAF'))
  if(!class(tem) %in% "try-error")
  {
    sample_num<-unique(tem$clusterData$MATH)
    names(sample_num) <-sample_name
    mutation_Heterogeneity<-c(mutation_Heterogeneity,sample_num)
  }
}

mutation_Heterogeneity_MD.HD<-data.frame(mutation_Heterogeneity)
write.table(mutation_Heterogeneity,file = "mutation_Heterogeneity_MD.HD.txt",sep="\t",quote = F)


mutation_Heterogeneity<-rbind(mutation_Heterogeneity_WT, mutation_Heterogeneity_LD, mutation_Heterogeneity_MD.HD)
