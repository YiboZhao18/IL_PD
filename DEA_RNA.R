# transcriptomics
setwd("/Users/zhaoyibo/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCL PostDoc 2024/ILproject_23/IL2023/focus/IL_Project")
condition <- read.delim("updated_pheno.txt")
#train_set <- condition[condition$split=="train",]

BL <- read.delim("rc_BL.txt", check.names = F)
BL_protein <- BL$protein
#BL_training <- BL[,intersect(train_set$patno, colnames(BL))]
#table(train_set$pheno[match(intersect(train_set$patno, colnames(BL)), train_set$patno)])
BL_patno_in <- colnames(BL)
BL_pheno_in <- condition$pheno[match(colnames(BL), condition$patno)]
# plot pheno count
ggplot(data.frame(pheno = BL_pheno_in), aes(x = pheno, fill = pheno))+
        geom_bar(stat = "count")+
        geom_text(stat = "count", aes(label = after_stat(count)),vjust = -0.5, size = 5)+
        theme_minimal()+
        labs(title = "Whole Blood Transcriptomics", x = " ", y = "Count")+
        theme(text = element_text(size = 14))

row.names(BL) <- BL_protein

library(DESeq2)
# load rc matrix-
count_matrix <- BL[,-1]

# read pheno data
condition_matrix <- data.frame(patno = BL_patno_in, pheno = BL_pheno_in)
condition_matrix$age <- scale(condition$age[match(BL_patno_in, condition$patno)])[,1]
condition_matrix$sex <- factor(condition$sex[match(BL_patno_in, condition$patno)])
#condition_matrix$HY <- factor(condition$HY[match(BL_patno_in, condition$patno)])
#condition_matrix$APOE <- factor(condition$APOE[match(BL_patno_in, condition$patno)])

# DESeq2
sample75 <- 0.75*dim(count_matrix)[2]
lowcount <- row.names(count_matrix)[!rowSums(count_matrix >=15) >= sample75]
gene_in <- row.names(count_matrix)[rowSums(count_matrix >=15) >= sample75]
count_matrix_in  <- count_matrix[rowSums(count_matrix >=15) >= sample75,]
condition_matrix$pheno <- as.factor(condition_matrix$pheno)
rownames(count_matrix_in) <- gene_in
# generate a normalised count matrix without design for WGCNA
#summary(is.na(condition_matrix$APOE))
#count_matrix_in <- count_matrix_in[,-match(condition_matrix$patno[is.na(condition_matrix$APOE)], colnames(count_matrix_in))]
condition_matrix <- as.data.frame(na.omit(condition_matrix))
library(DESeq2)
dds <- DESeqDataSetFromMatrix(count_matrix_in, condition_matrix, design = ~pheno+age+sex)
dds <- DESeq(dds)
# dds <- estimateSizeFactors(dds)
norm.counts <- counts(dds, normalized = T)
# sPD vs. Control
res1 <- results(dds, alpha = 0.05, contrast = c("pheno", "sPD", "Control")) # 1 is the reference level so the res is case vs. control
res1$padj <- p.adjust(res1$pvalue, "fdr")
summary(res1$padj<0.05)
write.table(res1, "updated_RNA_sPD_Control", sep = "\t", quote = F, row.names = T)

# LPD vs. Control
res2 <- results(dds, alpha = 0.05, contrast = c("pheno", "LPD", "Control")) # 1 is the reference level so the res is case vs. control
res2$padj <- p.adjust(res2$pvalue, "fdr")
summary(res2$padj<0.05)
write.table(res2, "updated_RNA_LPD_Control", sep = "\t", quote = F, row.names = T)

# LProdromal vs. Control
res3 <- results(dds, alpha = 0.05, contrast = c("pheno", "LPD", "sPD")) # 1 is the reference level so the res is case vs. control
summary(res3$padj<0.05)
write.table(res3, "updated_RNA_LPD_sPD", sep = "\t", quote = F, row.names = T)

library(ggplot2)
# sPD vs. Control
res1$sig <- rep("grey", nrow(res1))
res1$sig[res1$padj<0.05&res1$log2FoldChange>0] <- "indianred"
res1$sig[res1$padj<0.05&res1$log2FoldChange<0] <- "steelblue"
ggplot(data=res1, aes(x=log2FoldChange, y=-log10(pvalue), col=sig)) + 
        geom_point() + 
        theme_minimal() +
        # Add lines as before...
        #geom_vline(xintercept=c(-0.05, 0.1), col="black", linetype = "dashed") +
        geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed") +
        geom_vline(xintercept = 0, col = "grey")+
        geom_hline(yintercept=0, col="grey")+
        ## Change point color 
        scale_color_manual(values=c( "grey","indianred", "steelblue"))+
        xlim(-1.2, 1.2)+
        ylim(0,8)+
        theme(legend.position = "none")+
        labs(y = "-Log10(q)")


# LPD vs. Control
res2 <- get_sig_protein(comb2)
ggplot(data=res2, aes(x=log2FC, y=-log10(p), col=sig)) + 
        geom_point() + 
        theme_minimal() +
        # Add lines as before...
        #geom_vline(xintercept=c(0), col="black", linetype = "dashed") +
        geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
        ## Change point color 
        scale_color_manual(values=c("grey", "indianred", "steelblue"))+
        xlim(-0.1, 0.1)+
        ylim(0,8)+
        theme(legend.position = "none")+
        labs(y = "-Log10(q)")

# # LProdromal vs. Control
# res3 <- get_sig_protein(comb3)
# ggplot(data=res3, aes(x=log2FC, y=-log10(p), col=sig)) + 
#         geom_point() + 
#         theme_minimal() +
#         #geom_vline(xintercept=c(0), col="black", linetype = "dashed") +
#         geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
#         ## Change point color 
#         scale_color_manual(values=c("grey", "indianred", "steelblue"))+
#         xlim(-0.1, 0.1)+
#         ylim(0,8)+
#         theme(legend.position = "none")
