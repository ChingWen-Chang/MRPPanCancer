sample.mut.ref.WT<-read.delim("/Users/changc11/Documents/Pan_cancer_mutation2/mutation/dataset_WT.txt")
unique_genes_bysample <- sapply(split(sample.mut.ref.WT$Hugo_Symbol, sample.mut.ref.WT$Tumor_Sample_Barcode), unique)
# for each gene, calculate the number of samples with a mutation in the gene 
nrsamples_bygene <- table(unlist(unique_genes_bysample))
top10genes <- names(tail(sort(nrsamples_bygene), n=10))
top10genes
unique(sample.mut.ref.WT[,"Chromosome"])
sample.mut.ref.WT[,"Chromosome"] <- paste0("chr", sample.mut.ref.WT[,"Chromosome"])
sample.mut.ref.WT <- sample.mut.ref.WT[!grepl("chrM", sample.mut.ref.WT[,"Chromosome"]), ]
sample.mut.ref.WT <- sample.mut.ref.WT[!grepl("chrGL000212.1", sample.mut.ref.WT[,"Chromosome"]), ]
sample.mut.ref.WT <- sample.mut.ref.WT[!grepl("chrGL000219.1", sample.mut.ref.WT[,"Chromosome"]), ]
sample.mut.ref.WT <- sample.mut.ref.WT[!grepl("chrGL000205.1", sample.mut.ref.WT[,"Chromosome"]), ]
sample.mut.ref.WT <- sample.mut.ref.WT[!grepl("chrGL000209.1", sample.mut.ref.WT[,"Chromosome"]), ]

# filter mutations on some chromosome (e.g. chrM or "random"-chromosomes)
library(deconstructSigs)
trinucs_WT <- mut.to.sigs.input( sample.mut.ref.WT, 
                              sample.id="Tumor_Sample_Barcode", # names of the corresponding columns in data
                              chr="Chromosome", 
                              pos="Start_Position", 
                              ref="Reference_Allele", 
                              alt="Tumor_Seq_Allele2")


trinucs_WT_selected <- trinucs_WT[rowSums(trinucs_WT)>100,]
# initialize a list of the length of samples 

trinucs_WT_selected_t<-t(trinucs_WT_selected)
trinucs_WT_selected_t<-data.frame(trinucs_WT_selected_t)
results_WT<-lapply(split(trinucs_WT_selected_t, names(trinucs_WT_selected_t)), unname)
names(results_WT) <- row.names(trinucs_WT_selected)


# run the estimation of exposures for each sample and save the results_WT in the list
for( sID in row.names(trinucs_WT_selected) ){
  results_WT[[sID]] <- whichSignatures(trinucs_WT_selected, # the matrix generated with mut.to.sigs.input 
                                    sample.id=sID, # the current sample ID
                                    signatures.ref=signatures.cosmic, # the data.frame with the signatures that comes with deconstructSigs
                                    tri.counts.method="exome2genome", # which normalization method to use
                                    contexts.needed=TRUE) # set to TRUE if your input matrix contains counts instead of frequencies
}


makePie(results_WT[[4]])
plotSignatures(results_WT[[1]])

expo <- do.call("rbind", sapply(results_WT, "[", 1))
# add the unknown value to the matrix such that the contributions add up to 1 per sample
Signature.unknown <- unlist(sapply(results_WT, "[", 5))
expo <- cbind(expo, Signature.unknown)
write.table(expo,paste("Pan_cancer_WT_Sing_Patient_100",patient_name,sep = "",".txt"),quote = F,sep = "\t")
