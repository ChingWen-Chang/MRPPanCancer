# MRPPanCancer

Analysis Pipeline
Code was written by Ching-Wen Chang to ensure full reproducibility of the results reported in "Pan-cancer analysis uncovers a defective mitoribosomal biogenesis universally affecting tumor malignancy" by Ching-Wen Chang et al. 
For questions and/or comments, please contact changc11@nih.gov

Runtime environment: R 3.6.2

1. Identification and analysis of somatic MRPs copy number alterations 
Clinical and genomic data from TCGA, and MSK-IMPACT projects were downloaded from cBioPortal database (https://www.cbioportal.org/). The TCGA (phs000178.v10. p8) copy number data from GISTIC2, and log2 copy-number values are consist of 10,201 cancer and sample metadata from the NIH Genomic Data Commons (GDC). The MSK-IMPACT copy number data from GISTIC2 and log2 copy-number values consist of 7,976 cancer and sample metadata from Memorial Sloan Kettering Cancer Center. The CCLE copy number data from GISTIC2 and log2 copy-number values consist of 1,657 cancer and sample metadata. The CCLE copy number data was downloaded from CCLE database (https://portals.broadinstitute.org/ccle).
We determine copy number alteration thresholds according to the set of discrete copy number calls provided by GISTIC analysis: -2, deep loss/homozygous deletion; -1, shallow loss/heterozygous deletion; 0, diploid; 1, one copy gain; 2, high-level amplification. We focused on deletion altered samples with either deep loss (–2) or shallow loss (-1) of genes located in regions.
	Circle_plot_MRP.R: generates a relative circle plot of affected MRP genes, and chromosomes are shown clockwise in the outermost circle. The two innermost circles represent the correlation between gene expression and CNV and the deletion frequencies of affected MRP genes in each cancer type.
	Pan_cancer_clinical_group.R: performs comparison of conventional clinical parameters for subjects grouped according to the four MRP deletion status (HD, MD, LD, and WT) was evaluated using the chi-square test

2. Whole Exome Sequencing-based analysis
The mutation data (MAF files) based on whole exome sequencing is obtained from the cBioPortal database. The individual MAF files in different cancer types were concatenated into one combined MAF file for different MRPs defect groups for downstream analysis. Driver genes with statistically significant levels of recurrent mutation were determined by Mutation Significance (MutSig). Visualization and summarization of the mutations were performed by custom scripts in R version 3.6.2, primarily utilizing the packages maftools (version 2.4.05).
2.1 Oncoplot and driver gene analysis: 
	TCGA_mutation_plot.R: 
	Applies filters for variant selection to TCGA-Pan-Cancer, generates oncoplot and MRP deletion subcohorts (HD-MD, LD, and WT).
	Generates list of candidate driver genes using MutSigCV q-values.
	Performs pan-cancer analysis of mutation frequencies of driver genes for subcohorts based on MRP deletion status.
	Performs pan-cancer analysis of mutation frequencies of p53 hotspot loci for subcohorts based on MRP deletion status.
	Performs survival analysis (Kaplan-Meier, log-rank test) for subcohorts based on Wild type or MRPs deletion group (LD and HMD) alone or co-occurring with P53 mutations.
	mskcc_mutation_plot.R:
	Applies filters for variant selection to MSK-IMPACT-Pan-Cancer, generates oncoplot and MRP deletion subcohorts (HD-MD, LD, and WT).
	Generates list of candidate driver genes using MutSigCV q-values.
	Performs pan-cancer analysis of mutation frequencies of driver genes for subcohorts based on MRP deletion status.
	Performs pan-cancer analysis of mutation frequencies of p53 hotspot loci for subcohorts based on MRP deletion status.
	Performs survival analysis (Kaplan-Meier, log-rank test) for subcohorts based on Wild type or MRPs deletion group (LD and HMD) alone or co-occurring with P53 mutations.
	mafsurvival1_function.R: generates survival outcome table.

2.2 Mutational signature analysis:
	WT_Pancancer_mutsing_TCGA.R: generates matrix of subject vs signature weights using COSMIC (v2) for Wild type group.
	MD.LD_Pancancer_mutsing_TCGA.R: generates matrix of subject vs signature weights using COSMIC (v2) for MRPs deletion group (HMD).
	TCGA_Pancancer_mutsing_All_plot.R: generates bar chart of subject vs signature weights in TCGA-Pan-Cancer. And calculates signature weight distributions for
    the difference between MRPs deletion group (HMD) and Wild type or MRPs deletion group (HMD) co-occurring with P53 mutations and P53 mutations alone.
	mskcc_Pancancer_mutsing_All_plot.R: generates bar chart of subject vs signature weights in MSK-IMPACT-Pan-Cancer. And calculates signature weight
    distributions for the difference between MRPs deletion group (HMD) and Wild type or MRPs deletion group (HMD) co-occurring with P53 mutations and P53
    mutations alone.
2.3 Mutant-allele tumor heterogeneity analysis:
	Mutant-allele tumor heterogeneity.R: generates matrix of Mutant-Allele Tumor Heterogeneity (MATH) score.

3. Transcriptomics-based analysis
The gene expression from TCGA projects were downloaded from cBioPortal database (https://www.cbioportal.org/).
	TCGA_diversity_Score.R: The gene expression diversity of each TCGA patient is calculated with Shannon entropy (H(X)), using Equation 2 modified.

