
library(DESeq2)
library(dplyr)
library(limma)
library(readxl)
library(ggplot2)
library(forcats)
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
cData_filtered$seq_batch<-factor(cData_filtered$seq_batch)

# create a DESeq data structure from raw counts for statistical testing
dds = DESeqDataSetFromMatrix(countData=expressed_genes, colData=cData_filtered, design= ~tissue+seq_batch)
# choose and assign reference samples
dds$tissue = relevel(dds$tissue, ref="normal")

# run DESeq to test data
dds = DESeq(dds)

# report results to res variable
res = results(dds,alpha=0.05, contrast=c("tissue","normal","primary"))
summary(res)
sum(res$padj < 0.05, na.rm = TRUE)
# write full DESeq2 results to a table
write.table(res, "output/DESeq2_results_batch_correct.txt",sep="\t",quote=F)

#LFC shrinkage 
resLFC <- lfcShrink(dds, coef="tissue_primary_vs_normal", type="apeglm")

# plot log fold change as a function of transcript abundance, your data should be centered around 0 on the y-axis if it is well normalized
# alpha is set to 0.05 because of convention
DESeq2::plotMA(res,0.05,main='padj < 0.05 tissue+batch',ylim=c(-4,4))

DESeq2::plotMA(resLFC,0.05,main='padj < 0.05 log fold change shrinkage tissue+batch', ylim=c(-4,4))

# another normalization method with a variance stabilizing transform, a type of log transform
# good for comparing between datasets
vst_norm = vst(dds)


vst_corr_for_table <- assay(vst_norm)

# write a new vst normalized count table
write.table(vst_corr_for_table, "output/DESeq2_vst_normalized_ct_batch_correct.txt", sep="\t",quote=F,row.names=T,col.names=T)

### Count Table Batch Correction

# checking for batch effects in our vst normalized data

plotPCA(vst_norm,"tissue")
plotPCA(vst_norm,"seq_batch")
plotPCA(vst_norm,"passage")
plotPCA(vst_norm,"day")
plotPCA(vst_norm,"patient_id")

### Format res report
#order data based on logfold change 
res_no_format = res
res_no_format_log_fold_ordered = res_no_format[order(-res_no_format$log2FoldChange),]
#filter for significant and abundant genes
#get a boolean vector
res_of_interest_bool = !is.na(res_no_format_log_fold_ordered$padj) & res_no_format_log_fold_ordered$padj < 0.05 & res_no_format_log_fold_ordered$baseMean >= 100
#filter on boolean vector
res_of_interest = res_no_format_log_fold_ordered[res_of_interest_bool,]
#filter on boolean vector
DESeq2::plotMA(res_of_interest,0.05,main='padj < 0.05 tissue+batch',ylim=c(-4,4))

###### Realize the batch correction 
#actual batch correction
assay(vst_norm) = limma::removeBatchEffect(assay(vst_norm), vst_norm$seq_batch)

# check to look at the corrected structure of the data

plotPCA(vst_norm,"tissue")

plotPCA(vst_norm,"seq_batch")


pcaData <- plotPCA(vst_norm, intgroup=c("passage"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=factor(passage,levels=c('P1','P2','P3','P4','P7','P8','P9','P12')))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(color = "Passage") +scale_color_brewer(palette="Dark2")

library(RColorBrewer)
mycolors = c(brewer.pal(name="Set1", n = 9), brewer.pal(name="Pastel1", n = 5))
pcaData <- plotPCA(vst_norm, intgroup=c("patient_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=factor(patient_id,levels=c('R1N','R1T','R5N','R5T','R6N','R6T','R7T','R8T','R11N', 'R11T','R12N','R12T','R13T','R14T')))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + labs(color = "Patient") +scale_color_manual(values=mycolors)


vst_batch_corr_for_table <- assay(vst_norm)
write.table(vst_batch_corr_for_table, "output/DESeq2_vst_normalized_batch_corrected_ct.txt", sep="\t",quote=F,row.names=T,col.names=T)

