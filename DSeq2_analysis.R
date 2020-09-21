
library(DESeq2)
library(dplyr)
library(limma)
library(readxl)
## Load data
cts <- read.table("data/all.tsv", header=T, sep='\t')
annotation <- read_excel("data/elstar_colon_RNAseq_sample_sheet.xlsx")
cData <- data.frame(annotation)
rownames(cts)<-cts$gene
cts$gene <- NULL

#give name to rows based on their facility_id and patient_id
new_id_patient<-sub("N","C",cData$patient_id)
rownames(cData) <- paste(cData$facility_id,new_id_patient,sep="_")

#For the moment only interested in normal and primary samples 
cData_filtered<-filter(cData, tissue == 'normal' | tissue == 'primary')

#missing samples from count matrix 
noMissing<-rownames(cData_filtered)[rownames(cData_filtered) %in% colnames(cts)]

#excluding those samples 
cts_clean <- cts[,noMissing]
cData_filtered<-cData_filtered[colnames(cts_clean),]

#check the order
all(rownames(cData_filtered) %in% colnames(cts))
all(rownames(cData_filtered)==colnames(cts_clean))

# strip transcripts with expression lower than avg 10 counts per sample
expressed_genes <- cts_clean[rowSums(cts_clean)>80,]

cData_filtered$tissue<-factor(cData_filtered$tissue)

# create a DESeq data structure from raw counts for statistical testing
dds = DESeqDataSetFromMatrix(countData=expressed_genes, colData=cData_filtered, design= ~tissue)
# choose and assign reference samples
dds$tissue = relevel(dds$tissue, ref="normal")

# run DESeq to test data
dds = DESeq(dds)

# report results to res variable
res = results(dds,alpha=0.05)

# write full DESeq2 results to a table
write.table(res, "output/DESeq2_results.txt",sep="\t",quote=F)

#LFC shrinkage 
resLFC <- lfcShrink(dds, coef="tissue_primary_vs_normal", type="apeglm")

# plot log fold change as a function of transcript abundance, your data should be centered around 0 on the y-axis if it is well normalized
# alpha is set to 0.05 because of convention
DESeq2::plotMA(res,0.05,main='alpha = 0.05',ylim=c(-2,2))

DESeq2::plotMA(resLFC,0.05,main='alpha = 0.05 log fold change shrinkage')

# another normalization method with a variance stabilizing transform, a type of log transform
# good for comparing between datasets
vst_norm = vst(dds)


vst_corr_for_table <- assay(vst_norm)

# write a new vst normalized count table
write.table(vst_corr_for_table, "output/DESeq2_vst_normalized_ct.txt", sep="\t",quote=F,row.names=T,col.names=T)

### Count Table Batch Correction

# checking for batch effects in our vst normalized data

plotPCA(vst_norm,"batch")

plotPCA(vst_norm,"tissue")

