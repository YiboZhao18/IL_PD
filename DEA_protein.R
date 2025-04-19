# Proteomics
setwd("/Users/zhaoyibo/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCL PostDoc 2024/ILproject_23/IL2023/focus/IL_Project")
condition <- read.delim("updated_pheno.txt")
#train_set <- condition[condition$split=="train",]

# IL_list <- read.delim("IL_list.txt", stringsAsFactors = F)
# IL_protein <- IL_list$gene
# library(readxl)
# protein_all <- read_xlsx("ppmi_project_151.xlsx", sheet = "Results")
# probe_in <- colnames(protein_all)[-c(1,2)]
# protein_in <- integer()
# for (i in 1:length(probe_in)){
#          protein_in[i] <- unlist(strsplit(probe_in[i], split = "_"))[3]
#  }
# length(intersect(unique(IL_protein), protein_in))
# protein_IL <- data.frame(protein_all[,-c(1,2)])
# row.names(protein_IL) <- protein_all$PATNO
# write.table(protein_IL, "protein_IL.txt", sep = "\t", quote =  F, row.names = T)

protein_IL <- read.delim("protein_IL.txt")
probe_in <- colnames(protein_IL)
protein_in <- integer()
for (i in 1:length(probe_in)){
        protein_in[i] <- unlist(strsplit(probe_in[i], split = "_"))[3]
}

protein_patno <- rownames(protein_IL)
protein_pheno <- condition$pheno[match(protein_patno, condition$patno)]
table(protein_pheno)

comb1 <- protein_patno[protein_pheno=="Control"|protein_pheno=="sPD"]
comb2 <- protein_patno[protein_pheno=="Control"|protein_pheno=="LPD"]
#comb3 <- protein_patno_in[protein_pheno_in=="Control"|protein_pheno_in=="LProdromal"]
comb3 <- protein_patno[protein_pheno=="sPD"|protein_pheno=="LPD"]
#comb4 <- protein_patno_in[protein_pheno_in=="LPD"|protein_pheno_in=="LProdromal"]


get_sig_protein <- function(comb){
        protein_temp <- data.frame(protein_IL[comb,])
        group <- ifelse(protein_pheno[match(comb, protein_patno)]=="sPD", 0, 1)
        age <- condition$age[match(comb, condition$patno)]
        sex <- as.factor(condition$sex[match(comb, condition$patno)])
        #APOE <- as.factor(condition$APOE[match(comb, condition$patno)])
        p <- integer()
        coef <- integer()
        for (i in 1:ncol(protein_temp)){
                model <- data.frame(as.numeric(protein_temp[,i]), group, age, sex)
                colnames(model)[1] <- protein_in[i]
                res <- glm(group~., data = model, family = binomial)
                p[i] <- summary(res)$coefficients[2,4]
                coef[i] <- summary(res)$coefficients[2,1]
        }
        p.adj <- p.adjust(p, method = "fdr")
        sig <- ifelse(p.adj<0.05&coef>0, "indianred", ifelse(p.adj<0.05&coef<0, "steelblue", "grey"))
        outtable <- data.frame(coef, p, p.adj, sig)
        rownames(outtable) <- probe_in
        return(outtable)
}

# Volcano plots
res3 <- get_sig_protein(comb3)
library(ggplot2)
ggplot(data=res2, aes(x=coef, y=-log10(p), col=sig)) + 
        geom_point() + 
        theme_minimal() +
        # Add lines as before...
        #geom_vline(xintercept=c(-0.05, 0.1), col="black", linetype = "dashed") +
        geom_hline(yintercept=-log10(0.05), col="grey", linetype = "dashed") +
        geom_vline(xintercept = 0, col = "grey")+
        geom_hline(yintercept=0, col="grey")+
        ## Change point color 
        scale_color_manual(values=c("grey", "indianred", "steelblue"))+
        xlim(-10, 10)+
        ylim(0,15)+
        theme(legend.position = "none")+
        labs(x = "Coefficient",
             y = "-Log10(q)")



write.table(res1, "Updated_Protein_sPDvsHC.txt", sep = "\t", row.names = T, quote = F)
write.table(res2, "Updated_Protein_LPDvsHC.txt", sep = "\t", row.names = T, quote = F)
#write.table(res3, "Protein_LProdromalvsHC.txt", sep = "\t", row.names = T, quote = F)
write.table(res3, "Updated_Protein_LPDvssPD.txt", sep = "\t", row.names = T, quote = F)
