# Cohort summary
setwd("/Users/zhaoyibo/Library/CloudStorage/OneDrive-UniversityCollegeLondon/UCL PostDoc 2024/ILproject_23/IL2023/focus/IL_Project")
cohort <- read.delim("pheno.txt")
#head(cohort)
# keep healthy control, sporadic PD and LRRK2 PD patients (G2019S Only)
cohort_control <- cohort[cohort$COHORT_DEFINITION=="Healthy Control",]
cohort_sPD <- cohort[cohort$COHORT_DEFINITION=="PD"&
                             cohort$S.G=="Sporadic",]
cohort_LPD <- cohort[cohort$COHORT_DEFINITION=="PD"&
                             cohort$LRRK2==1&
                             cohort$GBA==0&
                             cohort$SNCA==0&
                             cohort$PRKN==0&
                             cohort$PINK1==0&
                             cohort$LRRK2.mutation=="G2019S",]
cohort_LProdromal <- cohort[cohort$COHORT_DEFINITION=="Prodromal"&
                                    cohort$LRRK2==1&
                                    cohort$GBA==0&
                                    cohort$SNCA==0&
                                    cohort$PRKN==0&
                                    cohort$PINK1==0&
                                    cohort$LRRK2.mutation=="G2019S",]

patno_all <- c(cohort_control$PATNO, cohort_sPD$PATNO, cohort_LPD$PATNO, cohort_LProdromal$PATNO)
condition <- data.frame(patno = patno_all, 
                        pheno = c(rep("Control", nrow(cohort_control)), 
                                  rep("sPD", nrow(cohort_sPD)), 
                                  rep("LPD", nrow(cohort_LPD)), 
                                  rep("LProdromal", nrow(cohort_LProdromal))), 
                        age = c(cohort_control$ENROLL_AGE,
                                cohort_sPD$ENROLL_AGE,
                                cohort_LPD$ENROLL_AGE,
                                cohort_LProdromal$ENROLL_AGE))


# compare baseline features: age, sex, cognitive function, H&Y stage, ethinicity
MoCA <- read.csv("MoCA.csv")
MoCA_BL <- as.data.frame(MoCA[MoCA$EVENT_ID=="SC",])
UPDRS3 <- read.csv("UPDRSIII.csv")
UPDRS3_BL <- UPDRS3[UPDRS3$EVENT_ID=="BL",]

condition$MoCA = rep(-9, nrow(condition))
condition$HY = rep(-9, nrow(condition))
for (i in 1:nrow(condition)){
        posi1 <- match(condition$patno[i], UPDRS3_BL$PATNO)
        if (!is.na(posi1)){
                condition$HY[i] = UPDRS3_BL$NHY[posi1]
        }
        posi2 <- match(condition$patno[i], MoCA_BL$PATNO)
        if(!is.na(posi2)){
                condition$MoCA[i] = MoCA_BL$MCATOT[posi2]
        }
} 
condition$MoCA[condition$MoCA<0] <- NA
condition$HY[condition$HY<0] <- NA
condition$HY[condition$HY>100] <- NA # found one with HY = 101
demographic <- read.csv("Demographics.csv")       
condition$sex <- demographic$SEX[match(condition$patno, demographic$PATNO)] # 0-Female; 1-Male      
condition$white <- demographic$RAWHITE[match(condition$patno, demographic$PATNO)]      
table(condition$white)      
condition <- condition[condition$white==1,]

# disease duration
duration <- read.csv("PD_Diagnosis_History.csv")
duration_mon <- integer()
convert_mon <- function(date){
        year <- as.integer(unlist(strsplit(date, split = "/", fixed = T))[2])
        mon <- as.integer(unlist(strsplit(date, split = "/", fixed = T))[1])
        cal_mon <- (year-2000)*12+mon
        return(cal_mon)
}
for (i in 1:nrow(duration)){
        symptom_start <- convert_mon(duration$SXDT[i])
        sc <- convert_mon(duration$INFODT[i])
        duration_mon[i] <- sc-symptom_start
        
}

condition$duration <- duration_mon[match(condition$patno, duration$PATNO)]
