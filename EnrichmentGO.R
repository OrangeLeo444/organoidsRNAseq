library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(geneplotter)

###Enrichment analysis

#res_no_format_log_fold_ordered <- this is the variable with my result from previous analysis


DESeq2Res <- res_no_format_log_fold_ordered

remove_decimal <- regexpr('ENSMUSG[0-9]*',rownames(DESeq2Res))
new_names <- regmatches(rownames(DESeq2Res), remove_decimal)

sigGenes <- rownames(subset(DESeq2Res, padj < 0.05))
anno <- AnnotationDbi::select(org.Mm.eg.db, 
                              keys=rownames(DESeq2Res), 
                              columns=c("SYMBOL", "GENENAME"),
                              keytype="ENSEMBL")

anSig <- as.data.frame(subset(anno, ENSEMBL %in% sigGenes))

sample_n(anSig, 5)

#identify background genes
overallBaseMean <- as.matrix(DESeq2Res[, "baseMean", drop = F])

sig_idx <- match(anSig$ENSEMBL, rownames(overallBaseMean))

backG <- c()

for(i in sig_idx){
  ind <- genefinder(overallBaseMean, i, 10, method = "manhattan")[[1]]$indices
  backG <- c(backG, ind)
  
}

backG <- unique(backG)
backG <- rownames(overallBaseMean)[backG]

backG <- setdiff(backG,  anSig$ENSEMBL)
length(backG)

multidensity( list( 
  all= log2(DESeq2Res[,"baseMean"]) ,
  foreground =log2(DESeq2Res[anSig$ENSEMBL, "baseMean"]), 
  background =log2(DESeq2Res[backG, "baseMean"])), 
  xlab="log2 mean normalized counts", main = "Matching for enrichment analysis",xlim = c(0,max(log2(DESeq2Res[,"baseMean"]))) )
