sample.mut.ref.MD.LD<-read.delim("/Users/changc11/Documents/Pan_cancer_mutation2/mutation/dataset_MD_LD.txt")
unique_genes_bysample <- sapply(split(sample.mut.ref.MD.LD$Hugo_Symbol, sample.mut.ref.MD.LD$Tumor_Sample_Barcode), unique)
# for each gene, calculate the number of samples with a mutation in the gene 
nrsamples_bygene <- table(unlist(unique_genes_bysample))
top10genes <- names(tail(sort(nrsamples_bygene), n=10))
top10genes
unique(sample.mut.ref.MD.LD[,"Chromosome"])
sample.mut.ref.MD.LD[,"Chromosome"] <- paste0("chr", sample.mut.ref.MD.LD[,"Chromosome"])
sample.mut.ref.MD.LD <- sample.mut.ref.MD.LD[!grepl("chrM", sample.mut.ref.MD.LD[,"Chromosome"]), ]
sample.mut.ref.MD.LD <- sample.mut.ref.MD.LD[!grepl("chrGL000212.1", sample.mut.ref.MD.LD[,"Chromosome"]), ]
sample.mut.ref.MD.LD <- sample.mut.ref.MD.LD[!grepl("chrGL000209.1", sample.mut.ref.MD.LD[,"Chromosome"]), ]
sample.mut.ref.MD.LD <- sample.mut.ref.MD.LD[!grepl("chrGL000205.1", sample.mut.ref.MD.LD[,"Chromosome"]), ]
sample.mut.ref.MD.LD <- sample.mut.ref.MD.LD[!grepl("chrGL000192.1", sample.mut.ref.MD.LD[,"Chromosome"]), ]
sample.mut.ref.MD.LD <- sample.mut.ref.MD.LD[!grepl("chrGL000213.1", sample.mut.ref.MD.LD[,"Chromosome"]), ]

# filter mutations on some chromosome (e.g. chrM or "random"-chromosomes)
library(deconstructSigs)
trinucs <- mut.to.sigs.input( sample.mut.ref.MD.LD, 
                              sample.id="Tumor_Sample_Barcode", # names of the corresponding columns in data
                              chr="Chromosome", 
                              pos="Start_Position", 
                              ref="Reference_Allele", 
                              alt="Tumor_Seq_Allele2")


trinucs_selected <- trinucs[rowSums(trinucs)>100,]
# initialize a list of the length of samples 
results <- vector("list", nrow(trinucs_selected))

trinucs_selected_t<-t(trinucs_selected)
trinucs_selected_t<-data.frame(trinucs_selected_t)
trinucs_selected_t[1:10,1:10]
results<-lapply(split(trinucs_selected_t, names(trinucs_selected_t)), unname)
names(results) <- row.names(trinucs_selected)


# run the estimation of exposures for each sample and save the results in the list
for( sID in row.names(trinucs_selected) ){
  results[[sID]] <- whichSignatures(trinucs_selected, # the matrix generated with mut.to.sigs.input 
                                          sample.id=sID, # the current sample ID
                                          signatures.ref=signatures.cosmic, # the data.frame with the signatures that comes with deconstructSigs
                                          tri.counts.method="exome2genome", # which normalization method to use
                                          contexts.needed=TRUE) # set to TRUE if your input matrix contains counts instead of frequencies
}


makePie(results[[1]])


expo <- do.call("rbind", sapply(results, "[", 1))
# add the unknown value to the matrix such that the contributions add up to 1 per sample
Signature.unknown <- unlist(sapply(results, "[", 5))
expo <- cbind(expo, Signature.unknown)
write.table(expo,paste("Pan_cancer_MD.LD_Sing_Patient_100",patient_name,sep = "",".txt"),quote = F,sep = "\t")

