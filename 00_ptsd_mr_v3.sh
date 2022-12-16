library(TwoSampleMR)
library(data.table)
library(cause)
library(mr.raps)
library(MRPRESSO)
library(dplyr)
library(MRlap)
library(MendelianRandomization)


 
# #Which analysis still needs to be done?
 # for file in $(cat trait_analysis_mrnames.txt)
 # do
   # found=$(ls results_5e-08 | grep -c complete_"$file"_PTSD.F3)
   # if [ $found -ge 1 ]
   # then
    # echo complete_"$file"_"PTSD.F3" found
   # else 
    # echo complete_"$file"_"PTSD.F3" not found!
    # echo complete_"$file"_"PTSD.F3" >> notfound.txt
   # fi
   
   # found=$(ls results_5e-08 | grep -c complete_PTSD.F3_"$file")
   # if [ $found -ge 1 ]
   # then
    # echo complete_"PTSD.F3"_"$file" found
   # else 
    # echo complete_"PTSD.F3"_"$file" not found!
    # echo complete_"PTSD.F3"_"$file" >> notfound.txt
   # fi
 # done
  


# #Which traits to analyze?
# for analysis in $(cat trait_numbers_analyze.txt) 
# do
 # file=$(ls sumstats_reformatted | grep ^"$analysis"_)
 # echo $file >> trait_analysis_filenames.txt
 # done
 
analyze <-scan(what=character())
135_SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz.txt
91_supar2021.fixed
87_dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.txt


141_GCST90027161-EFO_0000274-Build38.f.tsv.gz.txt


48_26394269-GCST003129-EFO_1001486.h.tsv.gz.txt
90_AITD2020.fixed
135_SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz.txt
91_supar2021.fixed
87_dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.txt


85_GCST90014023_buildGRCh38.tsv
42_EUR.UC.gwas_info03_filtered.assoc

133_ipscsg2016.result.combined.full.with_header.txt2
110_GCST90011871_buildGRCh37.tsv.gz.txt
134_EAGLE_AD_no23andme_results_29072015.txt

125_GCST90010715_buildGRCh37.tsv
112_discovery_metav3.0.meta.gz.txt
108_pernicious_anemia_Laisketal2021_sumstats.txt
129_Lopez-Isac_prePMID_META_GWAS_SSc.meta.txt
111_vitiligo.txt
94_IL6_summary_results_discovery_GWAS_MA_CHARGE.txt
52_UKBB.asthma.assoc
43_EUR.CD.gwas_info03_filtered.assoc
44_EUR.IBD.gwas_info03_filtered.assoc
57_RA_GWASmeta_European_v2.txt
54_GCST90018907_buildGRCh37.tsv
55_lupus.ea.imputed.allchr.out
59_C_reactive_protein.imp.gz.txt
18_BCX2_WBC_EA_GWAMA.out.gz.fuma.txt



#something wrong with asthma 3

#fileslist1=system('ls sumstats_reformatted/* ',intern=TRUE)
#fileslist=fileslist1[which(fileslist1 %in% analyze)]  #Change this if you add more
 #analyze <- fread('trait_analysis_filenames.txt',data.table=F,header=F)$V1
 
 #Set these to 1 if you want analysis performed
 do_main=0 #main MR analysis
 do_harmsave=0#save harmonized data
 do_harmload=0 #placeholder to just load harmonized data..
 do_mrpresso=0 #MR PRESSO
 do_mrlap=0 #MR LAP
 do_cause=1 #CAUSE analysis
 loops=0 # c(1,0) #c(1,0) #1 for PTSD -> trait, 0 for trait -> ptsd, c(0,1) for both
 
 
 # pthresh=5e-8
 pthresh=5e-8
 
#29:44 need a lower threshold (Sliz)
#49 does too (IL6)

#Do main MR of IBD with 

for (ptsdfile in c( "sumstats_reformatted/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt" )) #"sumstats_reformatted/77_eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.txt" ))#  ,"sumstats_reformatted/77_eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.txt"   "sumstats_reformatted/95_eur_ptsd_pcs_v4_aug3_2021_nomvpukbb.fuma.txt","sumstats_reformatted/78_eur_ptsd_pcs_v4_aug3_2021_nomvp.fuma.txt" ))#,)) #
{
 #Load PTSD data
  ptsdU_dat <- fread(ptsdfile,data.table=F)
 #remove HLA, MAPT inversion, CD40 gene |/- 3MB, CD40 is listed but not removed
  ptsdU_dat <- subset(ptsdU_dat, !(CHROM == 6 & POS >= 28477797 - 3000000 & POS <= 33448354 + 3000000) & 
                                 !(CHROM == 17 & POS >= 40928986 - 3000000 & POS <= 42139672 + 3000000)) #&
                              #  !(CHROM == 20 & POS >= 44746906 - 3000000 & POS <= 44758384  + 3000000)) 
 #Only strand unambiguous                            
  ptsdU_dat <- subset(ptsdU_dat, !((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A")  | (A1 == "C" & A2 == "G")  | (A1 == "G" & A2 == "C"))) #Remove ambiguous
 #PTSD filename
  ptsdset=ptsdU_dat$phenotype_col[1]

 #load disease files. Disease name is in phenotype file itself.
 for (diseasefile in analyze)
 {
  try({
  #I can 
  print (diseasefile)
  #Load dataset 
   disease_data <- fread(paste('sumstats_reformatted/',diseasefile,sep=""),data.table=F)
  #Filter bad regions - Make sure you are on hg19!
   disease_data <- subset(disease_data, !(CHROM == 6 & POS >= 28477797 - 3000000 & POS <= 33448354 + 3000000) & 
                                       !(CHROM == 17 & POS >= 40928986 - 3000000 & POS <= 42139672 + 3000000) )#&                                      #  !(CHROM == 20 & POS >= 44746906 - 3000000 & POS <= 44758384  + 3000000) ) #remove HLA, MAPT inversion, +/- 3MB
  #Remove ambiguous                                      
   disease_data <- subset(disease_data, !((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A")  | (A1 == "C" & A2 == "G")  | (A1 == "G" & A2 == "C"))) 
  
  #Disease filename
   disease=disease_data$phenotype_col[1]
 
  #Bidirectional testing
   for (loop in loops)
   {
    if (loop == 0)
    {
     exposurea <-   disease_data
     outcome <-   ptsdU_dat
     exposurefname=disease
     outcomefname=ptsdset
    }
     if (loop == 1)
    {
     exposurea <-   ptsdU_dat
     outcome <-      disease_data
     exposurefname=ptsdset
     outcomefname=disease
    }

  #Do main analysis based on thresholds
   
   if(do_main == 1 || do_harmsave==1 || do_mrpresso==1 || do_mrlap == 1 )
   {
   #Filter to just the good SNPs!
    threshsnps <- which(exposurea[,"P"] <= pthresh)
    exposure <- exposurea[threshsnps,]

   #Format exposure
    exposure_rf <- format_data(exposure,type="exposure",snp_col="SNP",beta_col="BETA",eaf_col="MAF",se_col="SE",effect_allele_col="A1",other_allele_col="A2",pval_col="P",phenotype_col="phenotype_col",samplesize_col="samplesize")

   #Filter exposure to only SNPs in outcome data
    exposure_rf <- subset(exposure_rf,!is.na(other_allele.exposure) & SNP %in% outcome$SNP) #prior to clumping, only take overlapping markers 
    
   #Clump exposure
    exposure_clump <- clump_data(exposure_rf, clump_r2 = 0.001)

   #Filter outcome to only clumped SNPs
    outcome2 <-  subset(outcome,SNP %in% exposure_clump$SNP )
    
   #Format outcome
    outcome_rf <- format_data(outcome2,type="outcome",snp_col="SNP",beta_col="BETA",se_col="SE",eaf_col="MAF",effect_allele_col="A1",other_allele_col="A2",pval_col="P",phenotype_col="phenotype_col",samplesize_col="samplesize")

   #harmonize data
    harmdat <- harmonise_data(exposure_clump ,outcome_rf,action=3)
   }
   #save data table of harmonized data, for later recreation without having to use that website. Need to build in being able to load this data
    if(do_harmsave == 1)
    {
     write.table(harmdat, file = paste('results_',pthresh,"/",'harmonizeddata_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""),quote=F)
    } 
    
    
    
    #main steps
    if(do_main ==1)
    {
     mr_report(harmdat,output_path=paste('results_',pthresh,sep=""))
     resultsmr <- mr(harmdat)
     hetmr <- mr_heterogeneity(harmdat)
     intmr <- mr_pleiotropy_test(harmdat)
     write.table(resultsmr, file = paste('results_',pthresh,"/",'mr_',  exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""),quote=F)
     write.table(hetmr, file = paste('results_',pthresh,   "/",'mrhet_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""),quote=F)
     write.table(intmr, file = paste('results_',pthresh,   "/",'mrint_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""),quote=F)
     
     resultsmrmain <- unlist(c(resultsmr[which(resultsmr$method == "MR Egger"),],
                                 resultsmr[which(resultsmr$method == "Weighted median"),],
                                 resultsmr[which(resultsmr$method == "Inverse variance weighted"),],
                                 hetmr[which(hetmr$method == "Inverse variance weighted"),],
                                 intmr))
     write.table(t(resultsmrmain), file=paste('results_',pthresh,"/",'basicmr_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""),quote=T,row.names=F)
   
    }
   #MR Presso
    if(do_mrpresso == 1)
    {
     mp <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = harmdat, NbDistribution = 10000,  SignifThreshold = 0.05, seed = 17)
     save(mp,file=paste('results_',pthresh,"/",'mrpresso_',exposurefname,"_",outcomefname,"_",pthresh,".r",sep=""))
     capture.output(mp, file = paste('results_',pthresh,"/",'mrpresso_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""))
     write.table( t(unlist(c(mp$`Main MR results`[2,],  
       mp$`MR-PRESSO results`$`Global Test`$Pvalue,
       mp$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`,
       mp$`MR-PRESSO results`$`Distortion Test`$Pvalue))),
       file=paste('results_',pthresh,"/",'mrpressoquick_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""),quote=T,row.names=F)
    } 
    
 #MRlap analysis
    if(do_mrlap == 1)
    {
     snplist <-fread('/mnt/ukbb/adam/tinnitus_gwas/w_hm3.noMHC.snplist',data.table=F)
     
     exposurea$N <- exposurea$samplesize
     exposurea$chr <- as.numeric(exposurea$CHROM)

     exp_pos <- subset(exposurea,select=c(SNP,POS))
     
     outcome$N <- outcome$samplesize
     outcome$chr <- as.numeric(outcome$CHROM)
     outcome$A1 <- toupper(outcome$A1)
     outcome$A2 <- toupper(outcome$A2)
     #This thing hates going across genome builds because the positions change... god damn it.
     #So just use the positions from the exposure.
     outcomea2 <- merge(outcome,exp_pos,by="SNP",suffixes=c("_dontuse",""))
     
     B = MRlap(exposure = exposurea,
               exposure_name = exposurefname,
               outcome =  outcomea2,
               outcome_name = outcomefname,
               ld = "/mnt/ukbb/royce/eur_w_ld_chr/",
               hm3 = "/mnt/ukbb/adam/tinnitus_gwas/w_hm3.noMHC.snplist",
               MR_threshold = pthresh, MR_pruning_dist=10000,
               MR_pruning_LD = 0.001)

     save(B,file=paste('results_',pthresh,"/",'mrlap_',exposurefname,"_",outcomefname,"_",pthresh,".r",sep=""))
     capture.output(B, file=paste('results_',pthresh,"/",'mrlap_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""))
     write.table(t(unlist(c(B$MRcorrection,B$LDSC )))    , file=paste('results_',pthresh,"/",'mrlap2_',exposurefname,"_",outcomefname,"_",pthresh,".txt",sep=""),quote=T,row.names=F)                                    
    }



   # CAUSE analysis (no p threshold needed)
   if(do_cause ==1)
   {
     X <- gwas_merge(  exposurea, outcome, snp_name_cols = c("SNP", "SNP"), 
                          beta_hat_cols = c("BETA", "BETA"), 
                          se_cols = c("SE", "SE"), 
                          A1_cols = c("A1", "A1"), 
                          A2_cols = c("A2", "A2"))
      #perform cause with 1 mil snps
      causesize=1000000
      if(dim(X)[1] < 1000000)
       {
            #use lesss if data too small
            causesize=500000
       }
     set.seed(100)
     varlist <- with(X, sample(snp, size=causesize, replace=FALSE))
     params <- est_cause_params(X, varlist)
     variants <- X %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

     prunelist <- NA
     for (chrom in c(1:22))
     {
      ld <- readRDS(paste("/mnt/ukbb/adam/tinnitus_gwas/cause_mr/chr",chrom,"_AF0.05_0.1.RDS",sep=""))
      snp_info <- readRDS(paste("/mnt/ukbb/adam/tinnitus_gwas/cause_mr/chr",chrom,"_AF0.05_snpdata.RDS",sep=""))
     
      prunelist <- c(prunelist,unlist(ld_prune(variants = variants, 
                                ld = ld, total_ld_variants = snp_info$SNP, 
                                pval_cols = c("pval1"), 
                                pval_thresh = c(1e-3))))
      }
      prunelist <- prunelist[-1] #get rid of the NA at the start
      
      res <- cause(X=X, variants=prunelist,param_ests=params)

      save(res,prunelist,params,file=paste('results/CAUSE_',exposurefname,"_",outcomefname,".r",sep=""))
      save(res,prunelist,params,file=paste('results/CAUSEres_',exposurefname,"_",outcomefname,".r",sep=""))
      
      res$elpd
      capture.output(summary(res, ci_size=0.95),file =paste('results/CAUSE_',exposurefname,"_",outcomefname,".txt",sep=""))
      
      pdf(paste('results/CAUSE_',exposurefname,"_",outcomefname,".pdf",sep=""),10,7)
      plot(res)
      #plot(res,type="data")
      dev.off()

    }
   }
   },silent=TRUE)
  }
 }


 #gamma = causal effect
 #eta = correlated pleiotropy (the bad thing that other models cant account for)
 #sharing model has gamma fixed at 0, allowing for horizontal pleitoropy but no cause
 # causal model allows gamma to vary
 #I think the null model is just of uncorrelated pleiotropy (direct effects of variants on outcome uncorrelated with exposure)

 #q i proportion of variants that exhibit correlated pleiotropy.
 
 
