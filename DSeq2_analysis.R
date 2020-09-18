
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
expressed_genes <- cts_clean[rowSums(cts_clean)>10,]

cData_filtered$tissue<-factor(cData_filtered$tissue)

# create a DESeq data structure from raw counts for statistical testing
dds = DESeqDataSetFromMatrix(countData=expressed_genes, colData=cData_filtered, design= ~tissue)
# choose and assign reference samples
dds$tissue = relevel(dds$tissue, ref="normal")
