library(TwoSampleMR)
library(data.table)
library(mr.raps)
library(MRPRESSO)
library(dplyr)
library(MRlap)
library(MendelianRandomization)

#maybe just start with a snplist to filter all data quickly?
analyze <- scan(what=character())
135_SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz.txt


57_RA_GWASmeta_European_v2.txt


59_C_reactive_protein.imp.gz.txt
90_AITD2020.fixed
133_ipscsg2016.result.combined.full.with_header.txt2
91_supar2021.fixed
94_IL6_summary_results_discovery_GWAS_MA_CHARGE.txt
108_pernicious_anemia_Laisketal2021_sumstats.txt
87_dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.txt

54_GCST90018907_buildGRCh37.tsv
52_UKBB.asthma.assoc
18_BCX2_WBC_EA_GWAMA.out.gz.fuma.txt


 pthresh=5e-8
#Load PTSD data
  ptsdU_dat <- fread("sumstats_reformatted/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt",data.table=F)
 #remove HLA, MAPT inversion, CD40 gene |/- 3MB, CD40 is listed but not removed
  ptsdU_dat <- subset(ptsdU_dat, !(CHROM == 6 & POS >= 28477797 - 3000000 & POS <= 33448354 + 3000000) & 
                                 !(CHROM == 17 & POS >= 40928986 - 3000000 & POS <= 42139672 + 3000000) ) # &
                               # !(CHROM == 20 & POS >= 44746906 - 3000000 & POS <= 44758384  + 3000000)) 
 #Only strand unambiguous                            
  ptsdU_dat <- subset(ptsdU_dat, !((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A")  | (A1 == "C" & A2 == "G")  | (A1 == "G" & A2 == "C"))) #Remove ambiguous
 #PTSD filename
  ptsdset=ptsdU_dat$phenotype_col[1]

 #read in exp2
  exp2_data <- fread("sumstats_reformatted/23_BCX2_EOS_EA_GWAMA.out.gz.fuma.txt",data.table=F)
 #Filter bad regions - Make sure you are on hg19!
  exp2_data <- subset(exp2_data, !(CHROM == 6 & POS >= 28477797 - 3000000 & POS <= 33448354 + 3000000) & 
                                       !(CHROM == 17 & POS >= 40928986 - 3000000 & POS <= 42139672 + 3000000) )#&                                      #  !(CHROM == 20 & POS >= 44746906 - 3000000 & POS <= 44758384  + 3000000) ) #remove HLA, MAPT inversion, +/- 3MB
 #Remove ambiguous                                      
  exp2_data <- subset(exp2_data, !((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A")  | (A1 == "C" & A2 == "G")  | (A1 == "G" & A2 == "C"))) 
 
 #PTSD only
  exp2_data <- subset(exp2_data, SNP %in% ptsdU_dat$SNP)
 #Disease filename
  mediator_name=exp2_data$phenotype_col[1]


 #load disease files. Disease name is in phenotype file itself.
 for (diseasefile in analyze)
 {
  #I can 
 # diseasefile="54_GCST90018907_buildGRCh37.tsv"
  print (diseasefile)
  #Load dataset 
   disease_data <- fread(paste('sumstats_reformatted/',diseasefile,sep=""),data.table=F)
  #Filter bad regions - Make sure you are on hg19!
   disease_data <- subset(disease_data, !(CHROM == 6 & POS >= 28477797 - 3000000 & POS <= 33448354 + 3000000) & 
                                       !(CHROM == 17 & POS >= 40928986 - 3000000 & POS <= 42139672 + 3000000) )#&                                      #  !(CHROM == 20 & POS >= 44746906 - 3000000 & POS <= 44758384  + 3000000) ) #remove HLA, MAPT inversion, +/- 3MB
  #Remove ambiguous                               
  #A1 != A2 is an additional step because some snps just harmonized poorly when I used 1000g refs for a1 and a2  
   disease_data <- subset(disease_data, A1 != A2 & !((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A")  | (A1 == "C" & A2 == "G")  | (A1 == "G" & A2 == "C"))) 
  
  #Disease filename
   disease=disease_data$phenotype_col[1]
 

  exposurea <- ptsdU_dat
  mediator <- exp2_data
  outcome <- disease_data
  exposurefname=ptsdset
  mediatorfname=mediator_name
  outcomefname=disease
  

  #Do main analysis based on thresholds

   #Filter to just the good SNPs!
    threshsnps <- which(exposurea[,"P"] <= pthresh)
    exposure <- exposurea[threshsnps,]

   #Format exposure - done here as outocme a programming trick
    exposure_rf <- format_data(exposure,type="outcome",snp_col="SNP",beta_col="BETA",eaf_col="MAF",se_col="SE",effect_allele_col="A1",other_allele_col="A2",pval_col="P",phenotype_col="phenotype_col",samplesize_col="samplesize")

   #Filter exposure to only SNPs in outcome data and in mediator
    exposure_rf <- subset(exposure_rf,!is.na(other_allele.outcome)  & SNP %in% disease_data$SNP  & SNP %in% mediator$SNP ) #prior to clumping, only take overlapping markers
     # outcome <- disease_data #disease data works but not outcome, what the fuck is going on?
  
    
   #Clump exposure
    exposure_clump <- clump_data(exposure_rf, clump_r2 = 0.001)

   #Filter outcome to only clumped SNPs
    outcome2 <-  subset(outcome,SNP %in% exposure_clump$SNP)
  
   #Filter mediator to only clumped SNPS
    mediator2 <-  subset(mediator,SNP %in% exposure_clump$SNP )
 
   #Format mediator
    mediator_rf <- format_data(mediator2,type="exposure",snp_col="SNP",beta_col="BETA",se_col="SE",eaf_col="MAF",effect_allele_col="A1",other_allele_col="A2",pval_col="P",phenotype_col="phenotype_col",samplesize_col="samplesize")
  
   #Format outcome
    outcome_rf <- format_data(outcome2,type="outcome",snp_col="SNP",beta_col="BETA",se_col="SE",eaf_col="MAF",effect_allele_col="A1",other_allele_col="A2",pval_col="P",phenotype_col="phenotype_col",samplesize_col="samplesize")
    

   #harmonize data Seems to harmonize on whatever the exposure is
    harm_exposure_mediator <- harmonise_data(mediator_rf,exposure_rf,action=3) #it is fine to use hte e
    harm_exposure_outcome <- harmonise_data(mediator_rf,outcome_rf,action=3)

      #Remember that in harmonizing, it says exposure and then outcome

    #mediator_beta_mat2 <- subset(harm_exposure_outcome,select=c(SNP,beta.exposure))
    
   # proof - alleles harnmonize the same way when the same variable is set as the exposure across datasets
   # table(mediator_beta_mat == mediator_beta_mat2)
    
    exposure_beta_mat <- subset(harm_exposure_mediator,select=c(SNP,beta.outcome,se.outcome,effect_allele.exposure,other_allele.exposure)) #exposure is listed as the outcome here, due to data formatting
    mediator_beta_mat <- subset(harm_exposure_mediator,select=c(SNP,beta.exposure,se.exposure)) # mediator is listed as the exposure, again due to data formatting
    outcome_beta_mat <- subset(harm_exposure_outcome,select=c(SNP,beta.outcome,se.outcome))
    
    names(exposure_beta_mat) <- c("SNP","beta.exposure","se.exposure","effect_allele.exposure","other_allele.exposure")
    names(mediator_beta_mat) <- c("SNP","beta.mediator","se.mediator")
    names(outcome_beta_mat) <- c("SNP","beta.outcome","se.outcome")
    
    aligned1 <- merge(exposure_beta_mat,mediator_beta_mat,by="SNP") #must merge, some snps will not align across both datasets due to incompatible alleles.
    aligned <- merge(aligned1,outcome_beta_mat,by="SNP") #must merge, some snps will not align across both datasets due to incompatible alleles.
    
    bx_mat <- as.matrix(subset(aligned,select=c(beta.exposure,beta.mediator)))
    bxse_mat <- as.matrix(subset(aligned,select=c(se.exposure,se.mediator)))
    by_mat <- subset(aligned,select=c(beta.outcome))$beta.outcome
    byse_mat <- subset(aligned,select=c(se.outcome))$se.outcome
    
    

    effect_allele=aligned$effect_allele.exposure
    other_allele=aligned$other_allele.exposure
   
   mvmr_obj <- mr_mvinput(bx=bx_mat,bxse=bxse_mat,by=by_mat,byse=byse_mat,exposure=c("VarPut","AdjustmentVar"),outcome=outcomefname,snps=aligned$SNP,effect_allele=aligned$effect_allele.exposure,other_allele=aligned$other_allele.exposure)
  
    test1 <- mr_mvivw(mvmr_obj,model="random")
    
    #write.table(summary(test1),file=paste("results5e-8/mr_mvivwPTSDF3_",diseasefile,".txt",sep=""),row.names=F)
  print(diseasefile)
   print(test1) 
}
  
