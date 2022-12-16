library(data.table)

#flip maf based on what allele is coded
  library(plyr)
  flipcheck <- function(x)
  {
   betaval=as.numeric(x["MAF"])
  # print(betaval)
   snp=c(x["A1"],x["A2"])
   pivot=c(x["A1_ldref"],x["A2_ldref"])
  # print(snp)
  # print(pivot)
   flipped=mapvalues(snp, c("A", "C", "G", "T"),  c("T", "G", "C", "A"))

   if (pivot[1] == snp[1] & pivot[2] == snp[2])
   {
    #No changes need to be made, alleles already aligned
     #  print(paste("No changes for", x[,]$Study,sep=" "))
   } else if (pivot[1] == snp[2] & pivot[2] == snp[1]) 
   {
    #A1n and A2n need to be flipped, change sign of beta value
     betaval=(1-betaval)
     #print(paste("Direction flip for", x[,]$Study,sep=" "))
   } else if ( !(flipped[1] %in% pivot) |  !(flipped[2] %in% pivot)) 
   {
    #print(paste("Study is unflippable", x[,]$Study,sep=" "))
    #Check if allele might be flippable, if not, give NA value 
    betaval=NA
   } else if (flipped[1] == snp[1] &  flipped[2] == snp[2]) 
   {
      #print(paste("Wrong strand but correct flip for ", x[,]$Study,sep=" "))
     #It's on the wrong strand but nothing actually needs to be done
   } else if (flipped[1] == snp[2] &  flipped[2] == snp[1]) 
   {
    #It's on the wrong strand and needs to be flipped
    # print(paste("Direction and allele flip for", x[,]$Study,sep=" "))
     betaval=(1-betaval)
   } else {
   betaval=NA 
   }
   return(betaval)
   }
   
#1 CigarettesPerDay.WithoutUKB.txt.gz
 cigperday_dat <- fread('zcat sumstats/CigarettesPerDay.WithoutUKB.txt.gz',data.table=F)
 cigperday_dat$SNP <- cigperday_dat$RSID
 cigperday_dat$A2 <- cigperday_dat$REF
 cigperday_dat$A1 <- cigperday_dat$ALT
 cigperday_dat$MAF <- cigperday_dat$AF
 cigperday_dat$samplesize <- cigperday_dat$N #Need a column for sample size. If there is no column, just write in the overall N
 cigperday_dat$phenotype_col <- "Cigarettes.Per.Day.NoUKBB"  #Make a dummy column just giving the trait name
 cigperday_dat$P <- cigperday_dat$PVALUE
 
 cigperday_exp <- subset(cigperday_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(cigperday_exp,'sumstats_reformatted/1_CigarettesPerDay.WithoutUKB.txt',quote=F,row.names=F)
  
#2 DrinksPerWeek.WithoutUKB.txt.gz
 drinksperweek_dat <- fread('zcat sumstats/DrinksPerWeek.WithoutUKB.txt.gz',data.table=F)
 drinksperweek_dat$SNP <- drinksperweek_dat$RSID
 drinksperweek_dat$A2 <- drinksperweek_dat$REF
 drinksperweek_dat$A1 <- drinksperweek_dat$ALT
 drinksperweek_dat$MAF <- drinksperweek_dat$AF
 drinksperweek_dat$samplesize <- drinksperweek_dat$N 
 drinksperweek_dat$phenotype_col <- "Drinks.Per.Week.NoUKBB"  
 drinksperweek_dat$P <- drinksperweek_dat$PVALUE
 
 drinksperweek_exp <- subset(drinksperweek_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(drinksperweek_exp,'sumstats_reformatted/2_DrinksPerWeek.WithoutUKB.txt',quote=F,row.names=F)
 
#3 CigarettesPerDay.txt.gz
 cigperdayfull_dat <- fread('zcat sumstats/CigarettesPerDay.txt.gz',data.table=F)
 cigperdayfull_dat$SNP <- cigperdayfull_dat$RSID
 cigperdayfull_dat$A2 <- cigperdayfull_dat$REF
 cigperdayfull_dat$A1 <- cigperdayfull_dat$ALT
 cigperdayfull_dat$MAF <- cigperdayfull_dat$AF
 cigperdayfull_dat$samplesize <- cigperdayfull_dat$N 
 cigperdayfull_dat$phenotype_col <- "Cigarettes.Per.Day"  
 cigperdayfull_dat$P <- cigperdayfull_dat$PVALUE
 
 cigperdayfull_exp <- subset(cigperdayfull_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(cigperdayfull_exp,'sumstats_reformatted/3_CigarettesPerDay.txt',quote=F,row.names=F)
 

#4 DrinksPerWeek.WithoutUKB.txt.gz
 drinksperweekfull_dat <- fread('zcat sumstats/DrinksPerWeek.txt.gz',data.table=F)
 drinksperweekfull_dat$SNP <- drinksperweekfull_dat$RSID
 drinksperweekfull_dat$A2 <- drinksperweekfull_dat$REF
 drinksperweekfull_dat$A1 <- drinksperweekfull_dat$ALT
 drinksperweekfull_dat$MAF <- drinksperweekfull_dat$AF
 drinksperweekfull_dat$samplesize <- drinksperweekfull_dat$N 
 drinksperweekfull_dat$phenotype_col <- "Drinks.Per.Week"  
 drinksperweekfull_dat$P <- drinksperweekfull_dat$PVALUE
 
 drinksperweekfull_exp <- subset(drinksperweekfull_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(drinksperweekfull_exp,'sumstats_reformatted/4_DrinksPerWeek.txt',quote=F,row.names=F)
 
#5 AUTOMOBILE_SPEEDING_PROPENSITY_GWAS.txt
 speeding_dat <- fread('sumstats/AUTOMOBILE_SPEEDING_PROPENSITY_GWAS.txt',data.table=F)
 speeding_dat$SNP <- speeding_dat$MarkerName
 speeding_dat$CHROM <- speeding_dat$CHR
 speeding_dat$MAF <- speeding_dat$EAF_A1
 speeding_dat$samplesize <- 404291 
 speeding_dat$phenotype_col <- "Automobile.Speeding.Propensity"  #Make a dummy column just giving the trait name
 speeding_dat$P <- speeding_dat$Pval
 speeding_dat$BETA <- -speeding_dat$Beta # it is reverse coded..
 
 speeding_exp <- subset(speeding_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(speeding_exp,'sumstats_reformatted/5_AUTOMOBILE_SPEEDING_PROPENSITY_GWAS.txt',quote=F,row.names=F)


#6 RISK_GWAS_MA_UKB+replication.txt
 risktol_dat <- fread('sumstats/RISK_GWAS_MA_UKB+replication.txt',data.table=F)
 risktol_dat$SNP <- risktol_dat$MarkerName
 risktol_dat$CHROM <- risktol_dat$CHR
 risktol_dat$MAF <- risktol_dat$EAF_A1
 risktol_dat$samplesize <-  466571 
 risktol_dat$phenotype_col <- "General.Risk.Tolerance"  #Make a dummy column just giving the trait name
 risktol_dat$P <- risktol_dat$Pval
 risktol_dat$BETA <- risktol_dat$Beta
 
 risktol_exp <- subset(risktol_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(risktol_exp,'sumstats_reformatted/6_RISK_GWAS_MA_UKB+replication.txt',quote=F,row.names=F)

#7 DRINKS_PER_WEEK_GWAS.txt
 drinksweek_dat <- fread('sumstats/DRINKS_PER_WEEK_GWAS.txt',data.table=F)
 drinksweek_dat$SNP <- drinksweek_dat$MarkerName
 drinksweek_dat$CHROM <- drinksweek_dat$CHR
 drinksweek_dat$MAF <- drinksweek_dat$EAF_A1
 drinksweek_dat$samplesize <-  414343
 drinksweek_dat$phenotype_col <- "Drinks.Per.Week.Karlson2019"  #Make a dummy column just giving the trait name
 drinksweek_dat$P <- drinksweek_dat$Pval
 drinksweek_dat$BETA <- drinksweek_dat$Beta
 
 drinksweek_exp <- subset(drinksweek_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(drinksweek_exp,'sumstats_reformatted/7_DRINKS_PER_WEEK_GWAS.txt',quote=F,row.names=F)


#8 EVER_SMOKER_GWAS_MA_UKB+TAG.txt
 smokinit_dat <- fread('sumstats/EVER_SMOKER_GWAS_MA_UKB+TAG.txt',data.table=F)
 smokinit_dat$SNP <- smokinit_dat$MarkerName
 smokinit_dat$CHROM <- smokinit_dat$CHR
 smokinit_dat$MAF <- smokinit_dat$EAF_A1
 smokinit_dat$samplesize <-  518633
 smokinit_dat$phenotype_col <- "Ever.Smoker.Karlson2019"  #Make a dummy column just giving the trait name
 smokinit_dat$P <- smokinit_dat$Pval
 smokinit_dat$BETA <- smokinit_dat$Beta
 
 smokinit_exp <- subset(smokinit_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(smokinit_exp,'sumstats_reformatted/8_EVER_SMOKER_GWAS_MA_UKB+TAG.txt',quote=F,row.names=F)

#9 NUMBER_SEXUAL_PARTNERS_GWAS.txt
 numsexpartners_dat <- fread('sumstats/NUMBER_SEXUAL_PARTNERS_GWAS.txt',data.table=F)
 numsexpartners_dat$SNP <- numsexpartners_dat$MarkerName
 numsexpartners_dat$CHROM <- numsexpartners_dat$CHR
 numsexpartners_dat$MAF <- numsexpartners_dat$EAF_A1
 numsexpartners_dat$samplesize <-  370711
 numsexpartners_dat$phenotype_col <- "Number.Sexual.Partners"  #Make a dummy column just giving the trait name
 numsexpartners_dat$P <- numsexpartners_dat$Pval
 numsexpartners_dat$BETA <- numsexpartners_dat$Beta
 
 numsexpartners_exp <- subset(numsexpartners_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(numsexpartners_exp,'sumstats_reformatted/9_NUMBER_SEXUAL_PARTNERS_GWAS.txt',quote=F,row.names=F)

#10 RISK_PC1_GWAS.txt
 riskpc1_dat <- fread('sumstats/RISK_PC1_GWAS.txt',data.table=F)
 riskpc1_dat$SNP <- riskpc1_dat$MarkerName
 riskpc1_dat$CHROM <- riskpc1_dat$CHR
 riskpc1_dat$MAF <- riskpc1_dat$EAF_A1
 riskpc1_dat$samplesize <-  315894
 riskpc1_dat$phenotype_col <- "Risk.PC1"  #Make a dummy column just giving the trait name
 riskpc1_dat$P <- riskpc1_dat$Pval
 riskpc1_dat$BETA <- riskpc1_dat$Beta
 
 riskpc1_exp <- subset(riskpc1_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(riskpc1_exp,'sumstats_reformatted/10_RISK_PC1_GWAS.txt',quote=F,row.names=F)

#11-25 Blood cells (European ancestry)
#Note: these were reformatted already to annotate positions, hence the FUMA appended. That is separate code.
cell_list <- scan(what=character())
RBC
HGB
HCT
MCH
MCV
MCHC
RDW
WBC
NEU
LYM
MON
BAS
EOS
PLT
MPV


 diseasenum=11 #cell types start with this number, I do naming the lazy way - loop at this point on the ordered list
 for (disease in cell_list)
 {
  print(c(disease,diseasenum))
   cell_dat <- fread(paste('sumstats/bloodcellnew/BCX2_',disease,'_EA_GWAMA.out.gz.fuma.txt',sep=''),data.table=F)
   outfname <- paste("sumstats_reformatted/",diseasenum,"_",'BCX2_',disease,'_EA_GWAMA.out.gz.fuma.txt',sep='')
   cell_dat$samplesize <- cell_dat$n_samples 
   cell_dat$phenotype_col <- paste(disease,".","European",sep="")  
   cell_dat$BETA = cell_dat$beta
   cell_dat$SE <- cell_dat$se

   cell_dat$A1 <- toupper(cell_dat$reference_allele)
   cell_dat$A2 <- toupper(cell_dat$other_allele)
   cell_dat$P <- cell_dat[,"p-value"]
   cell_dat$MAF <- cell_dat$eaf
   cell_dat$POS <- cell_dat$BP
   cell_exp <- subset(cell_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
   write.table(cell_exp,file=outfname,quote=F,row.names=F)
   diseasenum=diseasenum+1
  }
  
#26-40 Blood cells (Transethnic)


 diseasenum=26 #cell types start with this number, I do naming the lazy way - loop at this point on the ordered list
 for (disease in cell_list[-c(1:7)])
 {
  print(c(disease,diseasenum))
   cell_dat <- fread(paste('sumstats/bloodcellnew/BCX2_',disease,'_Trans_GWAMA.out.gz.fuma.txt',sep=''),data.table=F)
   outfname <- paste("sumstats_reformatted/",diseasenum,"_",'BCX2_',disease,'_Trans_GWAMA.out.gz.fuma.txt',sep='')
   cell_dat$samplesize <- cell_dat$n_samples 
   cell_dat$phenotype_col <- paste(disease,".","Transethnic",sep="")  
   cell_dat$BETA = cell_dat$beta
   cell_dat$SE <- cell_dat$se

   cell_dat$A1 <- toupper(cell_dat$reference_allele)
   cell_dat$A2 <- toupper(cell_dat$other_allele)
   cell_dat$P <- cell_dat[,"p-value"]
   cell_dat$MAF <- cell_dat$eaf
   cell_dat$POS <- cell_dat$BP
   cell_exp <- subset(cell_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
   write.table(cell_exp,file=outfname,quote=F,row.names=F)
   diseasenum=diseasenum+1
  }
  

#41 nallsEtAl2019_excluding23andMe_allVariants.tab

#Note: There was no SNP info in this file. I use chr and BP to give SNP names. Ideally, I would have used VCFtools and also used the allele codes. But this should be fine.
#Data is also filtered to a minimum N (N>25000).
#Only variants with rs ids are kept.

 #LC_ALL=C join <(awk '{if(NR== 1) $1= "CHR:BP"; if (NR == 1 || $8 >= 0.5*33674) print}' nallsEtAl2019_excluding23andMe_allVariants.tab | LC_ALL=C sort -k1b,1 ) <(awk '{print "chr"$2":"$3,$1}' /mnt/sdb/genetics/parkinsonsmr/phase3.locations2 | LC_ALL=C sort -k1b,1) > nallsEtAl2019_excluding23andMe_allVariants.tab.snp
 # echo CHROM BP A1 A2 MAF BETA SE P N_cases N_controls SNP > parkheader.txt
 # cat parkheader.txt <(cat nallsEtAl2019_excluding23andMe_allVariants.tab.snp | grep rs | dos2unix | sed 's/:/ /g' | sed 's/chr//g')  > nallsEtAl2019_excluding23andMe_allVariants.tab.snp2
  
 park_dat <- fread('sumstats/nallsEtAl2019_excluding23andMe_allVariants.tab.snp2',data.table=F)

 park_dat$samplesize <- park_dat$N_cases + park_dat$N_controls #Need a column for sample size. If there is no column, just write in the overall N
 park_dat$phenotype_col <- "Parkinsons"  #Make a dummy column just giving the trait name
 park_dat$POS <- park_dat$BP 
 parkinsons_exp <- subset(park_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(parkinsons_exp,file="sumstats_reformatted/41_nallsEtAl2019_excluding23andMe_allVariants.tab.snp2",quote=F,row.names=F)


#42 EUR.UC.gwas_info03_filtered.assoc
 uc_dat <- fread('EUR.UC.gwas_info03_filtered.assoc',data.table=F)
 uc_dat <- subset(uc_dat, INFO > 0.7)
 uc_dat$samplesize <- 6968 + 20464  #Need a column for sample size. If there is no column, just write in the overall N
 uc_dat$phenotype_col <- "Ulcerative.Colitis" #Make a dummy column just giving the trait name
 
 uc_dat$BETA = log(uc_dat[,"OR"])
 uc_dat$MAF <- uc_dat$FRQ_U_20464
 uc_dat$POS <- uc_dat$BP
 uc_dat$CHROM <- uc_dat$CHR
 
 uc_exp <- subset(uc_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(uc_exp,'sumstats_reformatted/42_EUR.UC.gwas_info03_filtered.assoc',quote=F,row.names=F)


#43 EUR.CD.gwas_info03_filtered.assoc
 cd_dat <- fread('sumstats/EUR.CD.gwas_info03_filtered.assoc',data.table=F)
 cd_dat <- subset(cd_dat, INFO > 0.7)
 cd_dat$samplesize <- 5956 + 14927  #Need a column for sample size. If there is no column, just write in the overall N
 cd_dat$phenotype_col <- "Crohns" #Make a dummy column just giving the trait name
 
 cd_dat$BETA = log(cd_dat[,"OR"])
 cd_dat$MAF <- cd_dat$FRQ_U_14927
 cd_dat$POS <- cd_dat$BP
 cd_dat$CHROM <- cd_dat$CHR
 
 cd_exp <- subset(cd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(cd_exp,'sumstats_reformatted/43_EUR.CD.gwas_info03_filtered.assoc',quote=F,row.names=F)



#44 EUR.IBD.gwas_info03_filtered.assoc
 ibd_dat <- fread('sumstats/EUR.IBD.gwas_info03_filtered.assoc',data.table=F)
 ibd_dat <- subset(ibd_dat, INFO > 0.7)
 ibd_dat$samplesize <- 12882 + 21770  #Need a column for sample size. If there is no column, just write in the overall N
 ibd_dat$phenotype_col <- "IBD" #Make a dummy column just giving the trait name
 
 ibd_dat$BETA = log(ibd_dat[,"OR"])
 ibd_dat$MAF <- ibd_dat$FRQ_U_21770
 ibd_dat$POS <- ibd_dat$BP
 ibd_dat$CHROM <- ibd_dat$CHR
 
 ibd_exp <- subset(ibd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ibd_exp,'sumstats_reformatted/44_EUR.IBD.gwas_info03_filtered.assoc',quote=F,row.names=F)

#45 UC_trans_ethnic_association_summ_stats_b37.txt.gz

 #I would have to apply metal or plink to this to get an IVW estimate
 #but first do sapply on the RMA package to see if thats fast enough.
 uctrans_dat <- fread('zcat sumstats/UC_trans_ethnic_association_summ_stats_b37.txt.gz',data.table=F)
 uctrans_dat$samplesize <- uctrans_dat$Mantra_n_samples  #Need a column for sample size. If there is no column, just write in the overall N
 uctrans_dat$phenotype_col <- "Ulcerative.Colitis.Transethnic" #Make a dummy column just giving the trait name

# SNP	1
# Chr	2
# Pos	3
# A1_effect	4
# A2_other	5
# Mantra_n_studies	6
# Mantra_log10BF	7
# Mantra_n_samples	8
# Mantra_dir	9
# beta_EUR	10
# se_EUR	11
# P_EUR	12
# beta_EAS	13
# se_EAS	14
# pval_EAS	15
# beta_IND	16
# se_IND	17
# pval_IND	18
# beta_IRA	19
# se_IRA	20
# pval_IRA	21
# EAF_EUR	22
# EAF_EAS	23
# EAF_IND	24
# EAF_IRA	25



#betas thaen ses:
#10,13,16,19,11,14,17,20


 library(metafor)
 rmaf <- function(x)
 {
 #We will do this strict, must be in ALL datasets to run!

  yis=x[1:4]
  
  seis=1/x[5:8]
  
  bm <- sum(yis*seis^2)/sum(seis^2)
  sm <- sqrt(1/sum(seis^2))
  
  if(length(yis)==0  | length(seis) == 0)
  {
  return(c(NA,NA))
  } else {
  return(c(bm,sm))
  }
 }
  
 #to insure this is errorproof, reorder the columns and apply only to the relevant ones..
 #rma(yi=unlist(uctrans_dat[1,c(10,13,16,19)]),sei=unlist(uctrans_dat[1,c(11,14,17,20)]),method="FE") #test single variable
 #apply(uctrans_dat[1:2,c(10,13,16,19,11,14,17,20)],1,rmaf)
 tesft <- as.data.frame(t(apply(uctrans_dat[,c(10,13,16,19,11,14,17,20)],1,rmaf)))
  names(tesft) <- c("beta","se")
 tesft$pval <- 2*pnorm(abs(tesft$beta/tesft$se),lower.tail=F)

 uctrans_dat$BETA = tesft$beta
 uctrans_dat$SE <- tesft$se
 uctrans_dat$P <- tesft$pval
 
 uctrans_dat$MAF <- uctrans_dat$EAF_EUR
 uctrans_dat$POS <- uctrans_dat$Pos
 uctrans_dat$CHROM <- uctrans_dat$Chr
 uctrans_dat$A1 <- uctrans_dat$A1_effect
 uctrans_dat$A2 <- uctrans_dat$A2_other
 
 uctrans_exp <- subset(uctrans_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(uctrans_exp,'sumstats_reformatted/45_UC_trans_ethnic_association_summ_stats_b37.txt',quote=F,row.names=F)

#46 Crohn's  (Transethnic)
#47 IBD (Crohns + UC)  (Transethnic)

#48 26394269-GCST003129-EFO_1001486.h.tsv.gz 
#result dataset is very small? 
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 pbc1_dat <- fread('zcat sumstats/26394269-GCST003129-EFO_1001486.h.tsv.gz',data.table=F)
 pbc1_dat$samplesize <- 2764 +  10475  #Need a column for sample size. If there is no column, just write in the overall N
 pbc1_dat$phenotype_col <- "Primary.Biliary.Cirrhosis" #Make a dummy column just giving the trait name
 
 pbc1_dat$SNP <- pbc1_dat$hm_rsid
 pbc1_dat$A1 <- pbc1_dat$effect_allele
 pbc1_dat$A2 <- pbc1_dat$other_allele
 pbc1_dat$BETA <- pbc1_dat$beta
 pbc1_dat$SE<- pbc1_dat$standard_error
 pbc1_dat$P <- pbc1_dat$p_value 
 pbc1_dat$CHROM <- pbc1_dat$hm_chrom
 pbc1_dat$POS <- pbc1_dat$hm_pos
 
 pbc_dat <- merge(pbc1_dat,AFs,by="SNP",suffixes=c("","_ldref")) #suffix must be 'ldref' for merging on this
 
 
 betaval <- try(apply(pbc_dat,1,flipcheck),silent=TRUE)
   
 pbc_dat$MAF <- betaval
 
 
 pbc_exp <- subset(pbc_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(pbc_exp,'sumstats_reformatted/48_26394269-GCST003129-EFO_1001486.h.tsv.gz.txt',quote=F,row.names=F)

#49 29273806-GCST006862-EFO_0000270-build37.f.tsv.gz
#UNSURE OF ALLELE CODING!

#There are not allele frequencies in this file. I will merge them in based on a reference file
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 asthma1_dat <- fread('sumstats/TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv',data.table=F)
 asthma1_dat$SNP <- asthma1_dat$rsid
 asthma1_dat$BETA = asthma1_dat$Multiancestry_beta_fix
 asthma1_dat$SE <- asthma1_dat$Multiancestry_se_fix
 asthma1_dat$phenotype_col <- "Asthma.Demenais2017.Transethnic"
 asthma1_dat$samplesize <- 142486
 
 asthma1_dat$A1 <- toupper(asthma1_dat$alternate_allele) #the corresopnding excel file confirms this
 asthma1_dat$A2 <- toupper(asthma1_dat$reference_allele)

 asthma1_dat$P <- asthma1_dat$Multiancestry_pval_fix
 
 asthma_dat <- merge(asthma1_dat,AFs,by="SNP",suffixes=c("","_ldref"))


   betaval <- try(apply(asthma_dat,1,flipcheck),silent=TRUE)
   
 asthma_dat$MAF <- betaval
 
 asthma_exp <- subset(asthma_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(asthma_exp,'sumstats_reformatted/49_TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv.transethnic.txt',quote=F,row.names=F)


#50 asthma european

#Just a change of beta, se and p from the last code:
 asthmaeur_dat <- asthma_dat
 asthmaeur_dat$BETA = asthmaeur_dat$European_ancestry_beta_fix
 asthmaeur_dat$SE <- asthmaeur_dat$European_ancestry_se_fix
 asthmaeur_dat$P <- asthmaeur_dat$European_ancestry_pval_fix
 asthmaeur_dat$phenotype_col <- "Asthma.Demenais2017.European"
 asthmaeur_dat$samplesize <- 127669

 asthmaeur_exp <- subset(asthmaeur_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(asthmaeur_exp,'sumstats_reformatted/50_TAGC_Multiancestry_and_European-Ancestry_Meta-analyses_Results.tsv.european.txt',quote=F,row.names=F)



#51 UKBB.allergy.assoc 
 allergy1_dat <- fread('sumstats/UKBB.allergy.assoc',data.table=F)

 allergy1_dat$SE <- allergy1_dat$standard_error
 allergy1_dat$phenotype_col <- "Cross.Allergy"
 allergy1_dat$samplesize <- 110361 

 allergy1_dat$CHROM <- allergy1_dat$CHR
 allergy1_dat$POS <- allergy1_dat$BP

 allergy1_dat$BETA <- log(allergy1_dat$OR)
 allergy1_dat$SE <- allergy1_dat$BETA/allergy1_dat$STAT
 allergy_dat <- merge(allergy1_dat,AFs,by="SNP",suffixes=c("","_ldref"))

 betaval <- try(apply(allergy_dat,1,flipcheck),silent=TRUE)
   
 allergy_dat$MAF <- betaval
 
 allergy_exp <- subset(allergy_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(allergy_exp,'sumstats_reformatted/51_UKBB.allergy.assoc',quote=F,row.names=F)

#52 UKBB.asthma.assoc.gz

 #NEEDS MAF INFORMATION (SAME FORMAT AS 51)
 #Asthma.Zhu2018
 
 asthma1_dat <- fread('zcat sumstats/UKBB.asthma.assoc.gz',data.table=F)

 asthma1_dat$SE <- asthma1_dat$standard_error
 asthma1_dat$phenotype_col <- "Asthma.Zhu2018"
 asthma1_dat$samplesize <- 90853

 asthma1_dat$CHROM <- asthma1_dat$CHR
 asthma1_dat$POS <- asthma1_dat$BP

 asthma1_dat$BETA <- log(asthma1_dat$OR)
 asthma1_dat$SE <- asthma1_dat$BETA/asthma1_dat$STAT
 asthma_dat <- merge(asthma1_dat,AFs,by="SNP",suffixes=c("","_ldref"))

 betaval <- try(apply(asthma_dat,1,flipcheck),silent=TRUE)
   
 asthma_dat$MAF <- betaval
 
 asthma_exp <- subset(asthma_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(asthma_exp,'sumstats_reformatted/52_UKBB.asthma.assoc',quote=F,row.names=F)

#54 psoriasis GCST90018907_buildGRCh37.tsv.gz
psoriasis_dat <- fread('zcat sumstats/GCST90018907_buildGRCh37.tsv.gz',data.table=F)
psoriasis_dat$phenotype_col <- "Psoriasis"
psoriasis_dat$CHROM <- psoriasis_dat$chromosome
psoriasis_dat$POS <- psoriasis_dat$base_pair_location
psoriasis_dat$A1 <- psoriasis_dat$effect_allele
psoriasis_dat$A2 <- psoriasis_dat$other_allele
psoriasis_dat$MAF <- psoriasis_dat$effect_allele_frequency
psoriasis_dat$BETA <- psoriasis_dat$beta
psoriasis_dat$SE <- psoriasis_dat$standard_error
psoriasis_dat$P <- psoriasis_dat$p_value

psoriasis_dat$samplesize <- 2472 +	346318 

psoriasis_dat$SNP <- #they forgot the snp...fuckers
psoriasis_dat$CHRBP <- paste(psoriasis_dat$CHROM,":",psoriasis_dat$POS,sep="")

 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)

 AFs$CHRBP <-  paste(AFs$CHROM,":",AFs$POS,sep="")
psoriasis_dat2 <- merge(psoriasis_dat,AFs,by="CHRBP",suffixes=c("","_AFS"))

 psoriasis_exp <- subset(psoriasis_dat2, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(psoriasis_exp,'sumstats_reformatted/54_GCST90018907_buildGRCh37.tsv',quote=F,row.names=F)

#55 lupus sle_immunochip9-22.zip
for chr in {1..22}
do
 cut -d " " -f 2,4,5,6,14,15,16,40,41,42 ea.imputed.chr"$chr".out  | tail -n+14    > lupus.ea.imputed.chr"$chr".out
done

cat lupus.ea.imputed.chr*.out | sort -g -k 6 > lupus.ea.imputed.allchr.out

#rs has : after it.

sed 's/:/ /g' lupus.ea.imputed.allchr.out | awk '{print $1}'  | paste -d " " - lupus.ea.imputed.allchr.out | grep rs    | grep -v completed    >  lupus.ea.imputed.allchr.out.rs
 

lupus_dat <- fread('sumstats/lupus.ea.imputed.allchr.out.rs',data.table=F,skip=21)

lupus_dat$phenotype_col <- "Lupus"
lupus_dat$A1 <- lupus_dat$alleleB #checked by looking up rs1132200 , https://www.nature.com/articles/ncomms16021/tables/2
lupus_dat$A2 <- lupus_dat$alleleA
lupus_dat$MAF <- lupus_dat$all_maf
lupus_dat$BETA <- lupus_dat[,"frequentist_add_beta_1:add/sle=1"]
lupus_dat$SE <- lupus_dat[,"frequentist_add_se_1"]
lupus_dat$P <- lupus_dat[,"frequentist_add_wald_pvalue_1"]
lupus_dat$SNP <- lupus_dat$rsid

lupus_dat$samplesize <- lupus_dat$all_total

 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)

 lupus_dat2 <- merge(lupus_dat,AFs,by="SNP",suffixes=c("","_ldref"))

 lupus_exp <- subset(lupus_dat2, MAF >= 0.01 & MAF <=0.99, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(lupus_exp,'sumstats_reformatted/lupus.ea.imputed.allchr.out',quote=F,row.names=F)



#skip 43 lines of lupus

#56 Trans-ethnic RA GWAS meta-analysis (19,234 RA cases and 61,565 conrols) [Link]


 rheum_dat1 <- fread('zcat sumstats/RA_GWASmeta_TransEthnic_v2.txt.gz',data.table=F)
 names(rheum_dat1)[1] <- "SNP"

 rheum_dat1$BETA = log(rheum_dat1[,"OR(A1)"])
 rheum_dat1$SE <- rheum_dat1$BETA / sqrt(qchisq(rheum_dat1[,"P-val"],1,lower.tail=F))

 rheum_dat1$P <- rheum_dat1[,"P-val"]
 rheum_dat1$phenotype_col <- "Rheumatoid.Arthritis.Transethnic"
 rheum_dat1$samplesize <- 19234 + 61565

 rheum_dat <- merge(rheum_dat1,AFs,by="SNP",suffixes=c("","_ldref"))
 betaval <- try(apply(rheum_dat,1,flipcheck),silent=TRUE)
 rheum_dat$MAF <- betaval
 
 
 rheum_exp <- subset(rheum_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(rheum_exp,'sumstats_reformatted/56_RA_GWASmeta_TransEthnic_v2.txt',quote=F,row.names=F)


#57 Eurpean RA GWAS meta-analysis (14,361 RA cases and 43,923 conrols) [Link]




for chr in {1..22}
do
 zcat sumstats/RA_GWASmeta_European_v2.txt.gz |  awk -v chr=$chr 'BEGIN{OFS="\t"}{beta=log($6); if (NR == 1) {$2= "chr";$3="bp";$1="MarkerName";$4="Allele1";$5="Allele2";beta="beta";$9="p"}; if(NR==1 || $2 == chr) print $1,$2,$3, $4,$5,beta,$9}'  > sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_"$chr"
done

for chr in {9..22} #remaining
do
/mnt/sdb/genetics/ardissvenv/ssimp_software-master/ssimp --gwas sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_"$chr"_imputed_chr"$chr" 
done

cat sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_*_imputed_chr* | awk '{if (NR==1 || $1 !="chr") print}' > sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_allchr_ssimp

R
f.z2b <- function(z, af, n)
{
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
    se.b <- 1/sqrt(2* af * (1-af) * n)
    b <- z * se.b
    return(c(b, se.b))
}
library(data.table)
output <- fread("sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_allchr_ssimp",data.table=F) ## reading ssimp output
N.max <- 52654 ## maximum sample size in the study
output$N_imp <- output$r2.pred * N.max ## calculating the effective sample size


## applying z-to-b function from above
output.b.se <- lapply(1:nrow(output), function(x) 
  f.z2b(output$z_imp[x], output$maf[x], output$N_imp[x])
)

## assigning b_imp and se_imp to dataframe
output$b_imp <- sapply(output.b.se, function(x) x[1])
output$se_imp <- sapply(output.b.se, function(x) x[2])

## store out extended output file
write.table(output, "sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_allchr_ssimp.betas",quote=F,row.names=F)




 rheumea_dat1 <- fread('zcat sumstats/RA_GWASmeta_European_v2.txt.gz',data.table=F)
 names(rheumea_dat1)[1] <- "SNP"


 rheumea_dat1$BETA = log(rheumea_dat1[,"OR(A1)"])
 rheumea_dat1$SE <- rheumea_dat1$BETA / sqrt(qchisq(rheumea_dat1[,"P-val"],1,lower.tail=F))

 rheumea_dat1$P <- rheumea_dat1[,"P-val"]
 rheumea_dat1$phenotype_col <- "rheumeaatoid.Arthritis.Europeanb"
 rheumea_dat1$samplesize <-  14361 + 43923

 rheumea_dat <- merge(rheumea_dat1,AFs,by="SNP",suffixes=c("","_ldref"))
 betaval <- try(apply(rheumea_dat,1,flipcheck),silent=TRUE)
 rheumea_dat$MAF <- betaval
 
 
 rheumea_exp <- subset(rheumea_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(rheumea_exp,'sumstats_reformatted/57_RA_GWASmeta_European_v2.txt',quote=F,row.names=F)



#58 1KG_CRP_GWAS_AJHG_2018.txt  #pivot errors in export

#needs maf
#awk '{FS=":"; OFS=" "; print $1,$2}' 1KG_CRP_GWAS_AJHG_2018.txt | awk '{CHRPOS=$1":"$2; if (NR==1) {$1="CHROM";$2="POS";CHRPOS="CHR:BP"}; print CHRPOS,$1,$2}' > 1KG_CRP_GWAS_AJHG_2018.txt_chrbp

 #LC_ALL=C join <(paste 1KG_CRP_GWAS_AJHG_2018.txt_chrbp 1KG_CRP_GWAS_AJHG_2018.txt | LC_ALL=C sort -k1b,1  ) <(awk '{print $2":"$3,$1}' /mnt/sdb/genetics/parkinsonsmr/phase3.locations2 | LC_ALL=C sort -k1b,1) | sort -g -k 7 > 1KG_CRP_GWAS_AJHG_2018.txt_rs
  AFs <- fread('zcat eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 
 crp2_dat <- fread('sumstats/1KG_CRP_GWAS_AJHG_2018.txt_rs',data.table=F)

 crp2_dat$samplesize <- 204402  #Need a column for sample size. If there is no column, just write in the overall N
 crp2_dat$phenotype_col <- "C.Reactive.Protein.Ligthart2018" #Make a dummy column just giving the trait name
 
 #crp2_dat$SNP <- crp2_dat$SNP
 crp2_dat$A1 <- toupper(crp2_dat$Allele1)
 crp2_dat$A2 <-  toupper(crp2_dat$Allele2)
 crp2_dat$BETA <- crp2_dat$Effect
 crp2_dat$SE<- crp2_dat$StdErr
 crp2_dat$P <- 2*pnorm(abs(crp2_dat$BETA/crp2_dat$SE),lower.tail=F)


 crp2a_dat <- merge(crp2_dat,AFs,by="SNP",suffixes=c("","_ldref"))

  
 betaval <- try(apply(crp2a_dat,1,flipcheck),silent=TRUE)
 crp2a_dat$MAF <- betaval
 
 
 crp2_exp <- subset(crp2a_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(crp2_exp,'sumstats_reformatted/58_1KG_CRP_GWAS_AJHG_2018.txt_rs',quote=F,row.names=F)


 
#59 C_reactive_protein.imp.gz

#LC_ALL=C join <(zcat C_reactive_protein.imp.gz | awk '{if(NR== 1) CHRBP = "CHR:BP"; else CHRBP =$1":"$2; print CHRBP,$0}' | LC_ALL=C sort -k1b,1 | sed 's/#//g' ) <(awk '{print $2":"$3,$1}' /mnt/sdb/genetics/parkinsonsmr/phase3.locations2 | LC_ALL=C sort -k1b,1) | sort -g -k 8 > C_reactive_protein.imp.gz.txt

 crp_dat <- fread('sumstats/C_reactive_protein.imp.gz.txt',data.table=F)

 crp_dat$samplesize <- 363000  #Need a column for sample size. If there is no column, just write in the overall N
 crp_dat$phenotype_col <- "C.Reactive.Protein.Armstrong2019" #Make a dummy column just giving the trait name
 
 
 crp_dat$A2 <- crp_dat$REF
 crp_dat$A1 <- crp_dat$ALT #THis file is coded in terms of alt, checked correspondce with other CRP GWAS that was in standard format to know this
 crp_dat$BETA <- crp_dat$Effect
 crp_dat$SE<- crp_dat$StdErr
 crp_dat$P <- as.numeric(crp_dat[,"P-value"])

 crp_exp <- subset(crp_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(crp_exp,'sumstats_reformatted/59_C_reactive_protein.imp.gz.txt',quote=F,row.names=F)


#sliz traits
# 60	VEGF_META
# 61	TNFa_META
# 62	MCP1_META
# 63	IP10_META
# 64	IL8_META
# 65	IL6_META
# 66	IL4_META
# 67	IL1RA_META
# 68	IL1B_META
# 69	IL17_META

cell_list <-scan(what=character())
VEGF
TNFa
MCP1
IP10
IL8
IL6
IL4
IL1RA
IL1B
IL17


 diseasenum=60 #biomarker types start with this number, I do naming the lazy way - loop at this point on the ordered list
 for (disease in cell_list)
 {
  print(c(disease,diseasenum))
   cell_dat <- fread(paste('sumstats/',disease,'_META_adj.age.sex.bmi.pc1-10.txt.gz',sep=''),data.table=F)
   outfname <- paste("sumstats_reformatted/",diseasenum,"_",disease,'_META_adj.age.sex.bmi.pc1-10.txt',sep='')
   cell_dat$SNP <- cell_dat$RSID
   cell_dat$samplesize <- cell_dat$N
   cell_dat$phenotype_col <- paste(disease,".Sliz",sep="")  
   cell_dat$MAF <- cell_dat$Freq1
   
   cell_dat2 <- merge(cell_dat,AFs,by="SNP",suffixes=c("","_ldref"))
   
   cell_exp <- subset(cell_dat2, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
   write.table(cell_exp,file=outfname,quote=F,row.names=F)
   diseasenum=diseasenum+1
  }
  
cell_list <-scan(what=character())
sVCAM1
sICAM1
sEselectin
sCD40L
IL1A
PAI1


 diseasenum=70 #biomarker types start with this number, I do naming the lazy way - loop at this point on the ordered list
 for (disease in cell_list)
 {
  print(c(disease,diseasenum))
   cell_dat <- fread(paste('sumstats/',disease,'_NFBC66_adj.age.sex.bmi.pc1-10.txt.gz',sep=''),data.table=F)
   outfname <- paste("sumstats_reformatted/",diseasenum,"_",disease,'_NFBC66_adj.age.sex.bmi.pc1-10.txt',sep='')
   cell_dat$SNP <- cell_dat$MARKER
   cell_dat$samplesize <- cell_dat$N
   cell_dat$phenotype_col <- paste(disease,".Sliz",sep="")  
   cell_dat$MAF <- cell_dat$EAF
   cell_dat$CHROM <- cell_dat$CHR
   cell_dat$A1 <- cell_dat$EA
   cell_dat$A2 <- cell_dat$NEA

   
   cell_exp <- subset(cell_dat, IMPUTATION_INFO > 0.6, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
   write.table(cell_exp,file=outfname,quote=F,row.names=F)
   diseasenum=diseasenum+1
  }
  



#76 PTSD Meta analysis
 ptsd_dat <- fread('zcat ../eur_ptsd_pcs_v4_aug3_2021.fuma.gz',data.table=F)
 #ptsd_dat <- subset(ptsd_dat, !(Chromosome == 6 & Position >= 25000000 & Position <= 35000000))
 ptsd_dat$samplesize <- ptsd_dat$Weight  
 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 ptsd_dat$BETA = ptsd_dat$Zscore * sed / sqrt(2*ptsd_dat$Weight*ptsd_dat$Freq1*(1-ptsd_dat$Freq1))
 ptsd_dat$SE <- ptsd_dat$B / ptsd_dat$Zscore 
 ptsd_dat$A1 <- toupper(ptsd_dat$Allele1)
 ptsd_dat$A2 <- toupper(ptsd_dat$Allele2)
 ptsd_dat$MAF <- ptsd_dat$Freq1
 ptsd_dat$P <- ptsd_dat[,"P-value"]
 ptsd_dat$SNP <- ptsd_dat$MarkerName
 ptsd_dat$POS <- ptsd_dat$Position
 ptsd_dat$CHROM <- ptsd_dat$Chromosome
 ptsd_dat$phenotype_col <- "PTSD.F3"
 
 ptsd_exp <- subset(ptsd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ptsd_exp,'sumstats_reformatted/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt',quote=F,row.names=F)

#77 ptsd meta, no ukbb
 ptsd_dat <- fread('zcat ../eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.gz',data.table=F)
 #ptsd_dat <- subset(ptsd_dat, !(Chromosome == 6 & Position >= 25000000 & Position <= 35000000))
 ptsd_dat$samplesize <- ptsd_dat$Weight  
 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 ptsd_dat$BETA = ptsd_dat$Zscore * sed / sqrt(2*ptsd_dat$Weight*ptsd_dat$Freq1*(1-ptsd_dat$Freq1))
 ptsd_dat$SE <- ptsd_dat$B / ptsd_dat$Zscore 
 ptsd_dat$A1 <- toupper(ptsd_dat$Allele1)
 ptsd_dat$A2 <- toupper(ptsd_dat$Allele2)
 ptsd_dat$MAF <- ptsd_dat$Freq1
 ptsd_dat$P <- ptsd_dat[,"P-value"]
 ptsd_dat$SNP <- ptsd_dat$MarkerName
 ptsd_dat$POS <- ptsd_dat$Position
 ptsd_dat$CHROM <- ptsd_dat$Chromosome
 ptsd_dat$phenotype_col <- "PTSD.F3.NoUKBB"
 
 ptsd_exp <- subset(ptsd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ptsd_exp,'sumstats_reformatted/77_eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.txt',quote=F,row.names=F)

#78 ptsd meta, no mvp
  ptsd_dat <- fread('zcat ../eur_ptsd_pcs_v4_aug3_2021_nomvp.fuma.gz',data.table=F)
# ptsd_dat <- subset(ptsd_dat, !(Chromosome == 6 & Position >= 25000000 & Position <= 35000000))
 ptsd_dat$samplesize <- ptsd_dat$Weight  
 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 ptsd_dat$BETA = ptsd_dat$Zscore * sed / sqrt(2*ptsd_dat$Weight*ptsd_dat$Freq1*(1-ptsd_dat$Freq1))
 ptsd_dat$SE <- ptsd_dat$B / ptsd_dat$Zscore 
 ptsd_dat$A1 <- toupper(ptsd_dat$Allele1)
 ptsd_dat$A2 <- toupper(ptsd_dat$Allele2)
 ptsd_dat$MAF <- ptsd_dat$Freq1
 ptsd_dat$P <- ptsd_dat[,"P-value"]
 ptsd_dat$SNP <- ptsd_dat$MarkerName
 ptsd_dat$POS <- ptsd_dat$Position
 ptsd_dat$CHROM <- ptsd_dat$Chromosome
 ptsd_dat$phenotype_col <- "PTSD.F3.NoMVP"
 
 ptsd_exp <- subset(ptsd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ptsd_exp,'sumstats_reformatted/78_eur_ptsd_pcs_v4_aug3_2021_nomvp.fuma.txt',quote=F,row.names=F)

#79 Educational attainment
 ea_dat <- fread('sumstats/education_GWAS_EA_excl23andMe.txt',data.table=F)

 ea_dat$samplesize <- 766345  #Need a column for sample size. If there is no column, just write in the overall N
 ea_dat$phenotype_col <- "Education" #Make a dummy column just giving the trait name
 
 ea_dat$BETA = ea_dat$Beta
 ea_dat$MAF <- ea_dat$EAF
 ea_dat$CHROM <- ea_dat$CHR
 ea_dat$SNP <- ea_dat$MarkerName
 ea_dat$P <- ea_dat$Pval
 
 ea_exp <- subset(ea_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ea_exp,'sumstats_reformatted/79_education_GWAS_EA_excl23andMe.txt',quote=F,row.names=F)


#80 cognitive performance
 cp_dat <- fread('sumstats/GWAS_CP_all.txt',data.table=F)

 cp_dat$samplesize <- 257828  #Need a column for sample size. If there is no column, just write in the overall N
 cp_dat$phenotype_col <- "Cognitive.Performance" #Make a dummy column just giving the trait name
 
 cp_dat$BETA = cp_dat$Beta
 cp_dat$MAF <- cp_dat$EAF
 cp_dat$CHROM <- cp_dat$CHR
 cp_dat$SNP <- cp_dat$MarkerName
 cp_dat$P <- cp_dat$Pval
 
 cp_exp <- subset(cp_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(cp_exp,'sumstats_reformatted/80_GWAS_CP_all.txt',quote=F,row.names=F)

#81 Irritable bowel disease 

 ibd2_dat <- fread('sumstats/GCST90016564_buildGRCh37.tsv',data.table=F)
 ibd2_dat$phenotype_col <- "IBD.2021"
 ibd2_dat$SNP <- ibd2_dat$variant_id
 ibd2_dat$CHROM <- ibd2_dat$chromosome
 ibd2_dat$POS <- ibd2_dat$base_pair_location
 ibd2_dat$A1 <- ibd2_dat$effect_allele
 ibd2_dat$A2 <- ibd2_dat$other_allele
 ibd2_dat$MAF <- ibd2_dat$effect_allele_frequency
 ibd2_dat$BETA <- ibd2_dat$beta
 ibd2_dat$SE <- ibd2_dat$standard_error
 ibd2_dat$P <- ibd2_dat$p_value
 ibd2_dat$samplesize <- ibd2_dat$N_CASE + ibd2_dat$N_CONTROL

 ibd2_exp <- subset(ibd2_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ibd2_exp,'sumstats_reformatted/81_GCST90016564_buildGRCh37.tsv',quote=F,row.names=F)

 #93 finnish biomarker GWAS - only has pvalues, cant use this.
 #please note, is in hg38
 inflam_dat <- fread('sumstats/GCST90000584_GRCh38.tsv',data.table=F)

#82 MS 24076602-GCST005531-EFO_0003885.h.tsv.gz . Only top associations! bullshit~!
 ms_dat <- fread('zcat sumstats/24076602-GCST005531-EFO_0003885.h.tsv.gz',data.table=F)

#84 graves disease
#rsid annotated
  graves_dat <- fread('zcat sumstats/GCST90018847_buildGRCh37.tsv.gz.2',data.table=F)
 graves_dat$phenotype_col <- "Graves"
 #graves_dat$SNP <- graves_dat$variant_id
 graves_dat$CHROM <- graves_dat$chromosome
 graves_dat$POS <- graves_dat$base_pair_location
 graves_dat$A1 <- graves_dat$effect_allele
 graves_dat$A2 <- graves_dat$other_allele
 graves_dat$MAF <- graves_dat$effect_allele_frequency
 graves_dat$BETA <- graves_dat$beta
 graves_dat$SE <- graves_dat$standard_error
 graves_dat$P <- 2*pnorm(abs(graves_dat$beta/graves_dat$standard_error),lower.tail=F) #just correct the p-value notation by calculating manually, original outputs had 1 set to 1.00E-00 which was annoying.
 
 graves_dat$CHRPOS <- paste(graves_dat$CHROM,graves_dat$POS,sep=":")
   AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 AFs$CHRPOS <- paste(AFs$CHROM,AFs$POS,sep=":")

 graves_dat$samplesize <- 340073 


 dim(graves_dat)
 

 graves_dat2 <- merge(graves_dat,AFs,by="CHRPOS",suffixes=c("","_ldref"))

 graves_exp <- subset(graves_dat2, MAF >=0.01 & MAF <= 0.99, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(graves_exp,'sumstats_reformatted/GCST90018847_buildGRCh37.tsv',quote=F,row.names=F)

 

#85 type 1 diabetes GCST90014023_buildGRCh38.tsv
#Positions need to be downgraded to hg37
 t1d_dat <- fread('sumstats/GCST90014023_buildGRCh38.tsv',data.table=F)
 t1d_dat$phenotype_col <- "T1D.2021"
 t1d_dat$SNP <- t1d_dat$variant_id
 #t1d_dat$CHROM <- t1d_dat$chromosome
 #t1d_dat$POS <- t1d_dat$base_pair_location
 t1d_dat$A1 <- t1d_dat$effect_allele
 t1d_dat$A2 <- t1d_dat$other_allele
 t1d_dat$MAF <- t1d_dat$effect_allele_frequency
 t1d_dat$BETA <- t1d_dat$beta
 t1d_dat$SE <- t1d_dat$standard_error
 t1d_dat$P <- 2*pnorm(abs(t1d_dat$beta/t1d_dat$standard_error),lower.tail=F) #just correct the p-value notation by calculating manually, original outputs had 1 set to 1.00E-00 which was annoying.
 

 t1d_dat$samplesize <- t1d_dat$sample_size
 t1d_dat <- subset(t1d_dat, sample_size > 400000)
 dim(t1d_dat)
 
  AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 t1d_dat2 <- merge(t1d_dat,AFs,by="SNP",suffixes=c("","_ldref"))

 t1d_exp <- subset(t1d_dat2, MAF >=0.01 & MAF <= 0.99, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(t1d_exp,'sumstats_reformatted/85_GCST90014023_buildGRCh38.tsv',quote=F,row.names=F)

#87 celiac disease 
#seems to be unimputed!!


for chr in {1..22}
do
 zcat sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz |  awk -v chr=$chr 'BEGIN{OFS="\t"}{if (NR == 1) {$1= "chr";$2="bp";$3="MarkerName";$5="Allele1";$4="Allele2";$7="beta";$6="p"}; if(NR==1 || $1 == chr) print $1,$2,$3, $5,$4,$7,$6}'  > sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_"$chr"
done

for chr in {1..22} #remaining
do
/mnt/sdb/genetics/ardissvenv/ssimp_software-master/ssimp --gwas sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_"$chr"_imputed_chr"$chr" 
done

cat sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_*_imputed_chr* | awk '{if (NR==1 || $1 !="chr") print}' > sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_imputed_allchr

R
f.z2b <- function(z, af, n)
{
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
    se.b <- 1/sqrt(2* af * (1-af) * n)
    b <- z * se.b
    return(c(b, se.b))
}
library(data.table)
output <- fread("sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_imputed_allchr",data.table=F) ## reading ssimp output
N.max <- 15283 ## maximum sample size in the study
output$N_imp <- output$r2.pred * N.max ## calculating the effective sample size


## applying z-to-b function from above
output.b.se <- lapply(1:nrow(output), function(x) 
  f.z2b(output$z_imp[x], output$maf[x], output$N_imp[x])
)

## assigning b_imp and se_imp to dataframe
output$b_imp <- sapply(output.b.se, function(x) x[1])
output$se_imp <- sapply(output.b.se, function(x) x[2])

## store out extended output file
write.table(output, "sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_imputed_allchr.betas",quote=F,row.names=F)



celiac_dat <- subset(output,N_imp > 10461.366)
 # fread('sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_imputed_allchr.betas',data.table=F)

celiac_dat$CHROM <- celiac_dat$chr
celiac_dat$POS <- celiac_dat$pos
celiac_dat$A1 <- celiac_dat$Allele1
celiac_dat$A2 <- celiac_dat$Allele2
celiac_dat$P <- celiac_dat$P.imp
celiac_dat$BETA <- celiac_dat$b_imp
celiac_dat$SE <- celiac_dat$se_imp
celiac_dat$samplesize <- 15283
celiac_dat$phenotype_col <- "Celiac"
celiac_dat$MAF <- celiac_dat$maf

 # AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 # AFs$MAF <-as.numeric(AFs$MAF)
 # celiac_dat2 <- merge(celiac_dat,AFs,by="SNP",suffixes=c("","_ldref"))
 # betaval <- try(apply(celiac_dat2,1,flipcheck),silent=TRUE)
 
 # celiac_dat2$MAF <- betaval
 

 celiac_exp <- subset(celiac_dat,  select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(celiac_exp,'sumstats_reformatted/87_dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz',quote=F,row.names=F)


# 90 autoimmune thyroid 
zcat sumstats/AITD2020.gz | awk '{if (NR>1) { gsub("chr","",$1)} print}' | awk '{if (NR>1) { gsub("na","NA")} print}' | awk '{if (NR==1 || $13 != "NA") print}' > sumstats/AITD2020.fixed


autothyroid_dat <- fread('sumstats/AITD2020.fixed',data.table=F)
#Remove any marker with SE = Inf

 autothyroid_dat$phenotype_col <- "Autoimmune.Thyroid"
 autothyroid_dat$SNP <- autothyroid_dat$rsID
 autothyroid_dat$CHROM <- autothyroid_dat$Chr
 autothyroid_dat$POS <- autothyroid_dat$Pos
 autothyroid_dat$A1 <- autothyroid_dat$A1
 autothyroid_dat$A2 <- autothyroid_dat$A0
 autothyroid_dat$MAF <- autothyroid_dat[,"IS-frq"]/100
 autothyroid_dat$BETA <- log(autothyroid_dat[,"OR-A1"])
 autothyroid_dat[which(autothyroid_dat$P == 1),]$P <- 0.999 #this function below cant handle p values of exactly 1 because it will make it infinite, so I modify P a bit..
 autothyroid_dat[which(autothyroid_dat[,"OR-A1"] == 1),]$BETA <- 0.001 #this function below cant handle berta of exactly 0 either.
 
 autothyroid_dat$SE <- abs((autothyroid_dat$BETA)/qnorm(autothyroid_dat$P/2)) 
 autothyroid_dat[which(autothyroid_dat$SE > 1),]$SE <- 1
 

 #autothyroid_dat$P <- 2*pnorm(abs(autothyroid_dat$BETA/autothyroid_dat$SE),lower.tail=F)
 #H
 autothyroid_dat$samplesize <- 755406
 dim(autothyroid_dat)
 
  # AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 # AFs$MAF <-as.numeric(AFs$MAF)
 # autothyroid_dat2 <- merge(autothyroid_dat,AFs,by="SNP",suffixes=c("","_ldref"))
 
  autothyroid_exp <- subset(autothyroid_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(autothyroid_exp,'sumstats_reformatted/90_AITD2020.fixed',quote=F,row.names=F)


#91 suPAR - unknown missing line
zcat sumstats/supar2021.gz | awk 'BEGIN{OFS="\t"}{if (NR>1) { gsub("chr","",$1)} print}' > sumstats/supar2021.fixed

supar_dat <- fread('sumstats/supar2021.fixed',data.table=F,header=T)
 supar_dat$phenotype_col <- "suPAR"
 supar_dat$SNP <- supar_dat$rsName
 #supar_dat$CHROM <- supar_dat$Chrom # I dont want HG38!
 #supar_dat$POS <- supar_dat$Pos_hg38
 supar_dat$A1 <- supar_dat$Effect_Allele
 supar_dat$A2 <- supar_dat$Other_Allele
 #supar_dat$MAF <- supar_dat[,"IS-frq"]
 supar_dat$BETA <- supar_dat$Effect
 supar_dat$SE <- supar_dat$StdErr
 supar_dat$P <- supar_dat$Pvalue

 supar_dat$samplesize <- supar_dat$N

 dim(supar_dat)
 
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 
 supar_dat2 <- merge(supar_dat,AFs,by="SNP",suffixes=c("","_ldref"))

 
 betaval <- try(apply(supar_dat2,1,flipcheck),silent=TRUE)
   
 supar_dat2$MAF <- betaval
 
 
 supar_exp <- subset(supar_dat2, MAF >=0.01 & MAF <= 0.99, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(supar_exp,'sumstats_reformatted/91_supar2021.fixed',quote=F,row.names=F)


# 94 IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz
#NOTE: readme says that these are HG18 POSITIONS! 

for chr in {1..22}
do
 zcat sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz |  awk -v chr=$chr 'BEGIN{OFS="\t"}{if (NR == 1) {$9= "chr";$10="bp";$1="MarkerName";$2="Allele1";$3="Allele2";$5="beta";$7="p"}; if(NR==1 || $9 == chr) print $1,$9,$10, $2,$3,$5,$7}'  > sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_"$chr"
done

for chr in {1..22} #remaining
do
/mnt/sdb/genetics/ardissvenv/ssimp_software-master/ssimp --gwas sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_"$chr"_imputed_chr"$chr" 
done

cat sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_*_imputed_chr* | awk '{if (NR==1 || $1 !="chr") print}' > sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_allchr_ssimp

R
f.z2b <- function(z, af, n)
{
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
    se.b <- 1/sqrt(2* af * (1-af) * n)
    b <- z * se.b
    return(c(b, se.b))
}
library(data.table)
output <- fread("sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_allchr_ssimp",data.table=F) ## reading ssimp output
N.max <- 52654 ## maximum sample size in the study
output$N_imp <- output$r2.pred * N.max ## calculating the effective sample size


## applying z-to-b function from above
output.b.se <- lapply(1:nrow(output), function(x) 
  f.z2b(output$z_imp[x], output$maf[x], output$N_imp[x])
)

## assigning b_imp and se_imp to dataframe
output$b_imp <- sapply(output.b.se, function(x) x[1])
output$se_imp <- sapply(output.b.se, function(x) x[2])

## store out extended output file
write.table(output, "sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_allchr_ssimp.betas",quote=F,row.names=F)


 
 il6_dat <-  fread('sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_allchr_ssimp.betas',data.table=F)
 il6_dat <- subset(il6_dat,N_imp >= 30000)
  # AFs <- fread('zcat eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
  # AFs <- subset(AFs,select=c(SNP,CHROM,POS))
  
 # il6_dat <- merge(il6_dat1,AFs,by="SNP",suffixes=c("_il6",""))
 # il6_dat$A1 <- il6_dat$EA
 il6_dat$CHROM <- il6_dat$chr
 il6_dat$POS <- il6_dat$pos
 il6_dat$A1 <- il6_dat$Allele1
 il6_dat$A2 <- il6_dat$Allele2
 il6_dat$P <- il6_dat$P.imp
 il6_dat$BETA <- il6_dat$b_imp
 il6_dat$SE <- il6_dat$se_imp
 il6_dat$samplesize <- 52654
 il6_dat$phenotype_col <- "il6"
 il6_dat$MAF <- il6_dat$maf

 il6_exp <- subset(il6_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(il6_exp,'sumstats_reformatted/94_IL6_summary_results_discovery_GWAS_MA_CHARGE.txt',quote=F,row.names=F)

#95 PTSD with no MVP or UKBB 
#Check this alternative beta:
* Beta/SE were calculated from METAL Z-scores using the formula from Zhu et al (Nature Genetics, 2016): - only minutely different from mine, z score^2 is small relative to N.
Beta = Zscore / sqrt( 2 * MAF * ( 1 - MAF) * ( N + Zscore^2 ) )
SE = 1 / sqrt( 2 * MAF * ( 1 - MAF ) * ( N + Zscore^2 ) )

 ptsd_dat <- fread('zcat sumstats/eur_ptsd_pcs_v4_aug3_2021_nomvpukbb.fuma.gz',data.table=F)
 ptsd_dat$samplesize <- ptsd_dat$Weight  
 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 ptsd_dat$BETA = ptsd_dat$Zscore * sed / sqrt(2*ptsd_dat$Weight*ptsd_dat$Freq1*(1-ptsd_dat$Freq1))
 ptsd_dat$SE <- ptsd_dat$B / ptsd_dat$Zscore 
 ptsd_dat$A1 <- toupper(ptsd_dat$Allele1)
 ptsd_dat$A2 <- toupper(ptsd_dat$Allele2)
 ptsd_dat$MAF <- ptsd_dat$Freq1
 ptsd_dat$P <- ptsd_dat[,"P-value"]
 ptsd_dat$SNP <- ptsd_dat$MarkerName
 ptsd_dat$POS <- ptsd_dat$Position
 ptsd_dat$CHROM <- ptsd_dat$Chromosome
 ptsd_dat$phenotype_col <- "PTSD.F3.NoMVPUKBB"
 
 ptsd_exp <- subset(ptsd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ptsd_exp,'sumstats_reformatted/95_eur_ptsd_pcs_v4_aug3_2021_nomvpukbb.fuma.txt',quote=F,row.names=F)

#96 MVP
 ptsd_dat <- fread('zcat /mnt/ukbb/adam/ptsd/ehr_reformat/TotalPCL_MVP_eur.gz',data.table=F)
 ptsd_dat$samplesize <- 186689 

 ptsd_dat$A1 <- toupper(ptsd_dat$Allele1)
 ptsd_dat$A2 <- toupper(ptsd_dat$Allele2)
 ptsd_dat$MAF <- ptsd_dat$Freq1

 ptsd_dat$SNP <- ptsd_dat$rsid
 ptsd_dat$POS <- ptsd_dat$BP
 ptsd_dat$CHROM <- ptsd_dat$CHR
 ptsd_dat$phenotype_col <- "PTSD.MVP"
 
 ptsd_exp <- subset(ptsd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ptsd_exp,'sumstats_reformatted/96_TotalPCL_MVP_eur',quote=F,row.names=F)

#98 MDD
 
 mdd_dat1 <- fread('sumstats/PGC_UKB_depression_genome-wide.txt',data.table=F)
 mdd_dat1$SNP <- mdd_dat$MarkerName

 AFs <- fread('zcat eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs <- subset(AFs,select=c(SNP,CHROM,POS))
  
 mdd_dat <- merge(mdd_dat1,AFs,by="SNP",suffixes=c("","_ldref"))


 mdd_dat$phenotype_col <- "MDD.Howard"

 #mdd_dat$A1 <- mdd_dat$A1
 #mdd_dat$A2 <- mdd_dat$A2
 mdd_dat$MAF <- mdd_dat$Freq
 mdd_dat$BETA <- mdd_dat$LogOR
 mdd_dat$SE <- mdd_dat$StdErrLogOR
 #mdd_dat$P <- mdd_dat$P
 mdd_dat$samplesize <- 500199

 mdd_exp <- subset(mdd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(mdd_exp,'sumstats_reformatted/98_PGC_UKB_depression_genome-wide.txt',quote=F,row.names=F)


#99	ADHD
 adhd_dat1 <- fread('zcat sumstats/adhd_eur_jun2017.gz',data.table=F)
 
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 adhd_dat1$samplesize <- 53293 #Need a column for sample size. If there is no column, just write in the overall N
 adhd_dat1$phenotype_col <- "ADHD.Demontis" #Make a dummy column just giving the trait name
 
 #adhd_dat$SNP <- adhd_dat$SNP
 #adhd_dat$A1 <- adhd_dat$A1
 #adhd_dat$A2 <- adhd_dat$A2
 adhd_dat1$BETA <- log(adhd_dat1$OR)
 #adhd_dat$SE<- adhd_dat$SE
 #adhd_dat$P <- adhd_dat$P

 adhd_dat <- merge(adhd_dat1,AFs,by="SNP",suffixes=c("","_ldref")) #suffix must be 'ldref' for merging on this
 
 
 betaval <- try(apply(adhd_dat,1,flipcheck),silent=TRUE)
   
 adhd_dat$MAF <- betaval
 adhd_exp <- subset(adhd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 
 write.table(adhd_exp,'sumstats_reformatted/99_adhd_eur_jun2017.txt',quote=F,row.names=F)

 
#100	Alcohol Dependence

 alcdep_dat1 <- fread('zcat sumstats/pgc_alcdep.eur_discovery.aug2018_release_FIXED.txt.gz',data.table=F)
 alcdep_dat1 <- alcdep_dat1[,1:8]
 
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 

 alcdep_dat1$samplesize <- 53293 #Need a column for sample size. If there is no column, just write in the overall N
 alcdep_dat1$phenotype_col <- "Alcohol.Dependence" #Make a dummy column just giving the trait name
 
 #alcdep_dat$SNP <- alcdep_dat$SNP
 #alcdep_dat$A1 <- alcdep_dat$A1
 #alcdep_dat$A2 <- alcdep_dat$A2
 #alcdep_dat$P <- alcdep_dat$P

 alcdep_dat <- merge(alcdep_dat1,AFs,by="SNP",suffixes=c("","_ldref")) #suffix must be 'ldref' for merging on this
 
 betaval <- try(apply(alcdep_dat,1,flipcheck),silent=TRUE) 
 alcdep_dat$MAF <- betaval
 
 alcdep_dat$BETA = alcdep_dat$Z  / sqrt(2*alcdep_dat$Weight*alcdep_dat$MAF*(1-alcdep_dat$MAF))
 alcdep_dat$SE <- alcdep_dat$BETA / alcdep_dat$Z
 
 

 alcdep_exp <- subset(alcdep_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 
 write.table(alcdep_exp,'sumstats_reformatted/100_pgc_alcdep.eur_discovery.aug2018_release.txt',quote=F,row.names=F)


#101	PTSD 2.5
 ptsd_dat <- fread('zcat sumstats/eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz',data.table=F)
 ptsd_dat$samplesize <- ptsd_dat$Weight  
 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 ptsd_dat$BETA = ptsd_dat$Zscore * sed / sqrt(2*ptsd_dat$Weight*ptsd_dat$Freq1*(1-ptsd_dat$Freq1))
 ptsd_dat$SE <- ptsd_dat$B / ptsd_dat$Zscore 
 ptsd_dat$A1 <- toupper(ptsd_dat$Allele1)
 ptsd_dat$A2 <- toupper(ptsd_dat$Allele2)
 ptsd_dat$MAF <- ptsd_dat$Freq1
 ptsd_dat$P <- ptsd_dat[,"P-value"]
 ptsd_dat$SNP <- ptsd_dat$MarkerName
 ptsd_dat$POS <- ptsd_dat$Position
 ptsd_dat$CHROM <- ptsd_dat$Chromosome
 ptsd_dat$phenotype_col <- "PTSD.F25"
 
 ptsd_exp <- subset(ptsd_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ptsd_exp,'sumstats_reformatted/101_eur_PTSD_Continuous_m0_pcs_alldata_may8_2020.fuma.gz',quote=F,row.names=F)


#102	Alzheimers
 alzheimers_dat1 <- fread('zcat sumstats/PGCALZ2sumstatsExcluding23andMe.txt.gz',data.table=F)


 alzheimers_dat1$samplesize <- alzheimers_dat1$N #Need a column for sample size. If there is no column, just write in the overall N
 alzheimers_dat1$phenotype_col <- "Alzheimers.2021" #Make a dummy column just giving the trait name
 
 alzheimers_dat1$CHROM <- alzheimers_dat1$chr
 alzheimers_dat1$POS <- alzheimers_dat1$PosGRCh37
 
 #alzheimers_dat$SNP <- alzheimers_dat$SNP
 alzheimers_dat1$A1 <- alzheimers_dat1$testedAllele
 alzheimers_dat1$A2 <- alzheimers_dat1$otherAllele 
 alzheimers_dat1$P <- alzheimers_dat1$p
 alzheimers_dat1[which(alzheimers_dat1$p == 0),]$P <- 6.464721e-304
 alzheimers_dat1$Z <- alzheimers_dat1$z
 alzheimers_dat1[which(alzheimers_dat1$z == Inf),]$Z <- 40
 alzheimers_dat1[which(alzheimers_dat1$z == -Inf),]$Z <- -40


 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 alzheimers_dat <- merge(alzheimers_dat1,AFs,by=c("CHROM","POS"),suffixes=c("","_ldref")) #suffix must be 'ldref' for merging on this
 
 betaval <- try(apply(alzheimers_dat,1,flipcheck),silent=TRUE) 
 alzheimers_dat$MAF <- betaval
 
 alzheimers_dat$BETA = alzheimers_dat$Z  / sqrt(2*alzheimers_dat$N*alzheimers_dat$MAF*(1-alzheimers_dat$MAF))
 alzheimers_dat$SE <- alzheimers_dat$BETA / alzheimers_dat$Z
 
 alzheimers_dat2 <- subset(alzheimers_dat,!is.na(MAF))
 which(is.na(alzheimers_dat2$BETA))
 
  alzheimers_exp <- subset(alzheimers_dat2, samplesize >= 0.7*max(alzheimers_dat$samplesize) & !(CHROM==19 & POS > (45337918 - 1000000)& POS  < (45337918 + 1000000)), select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
  
 write.table(alzheimers_exp,'sumstats_reformatted/102_PGCALZ2sumstatsExcluding23andMeR.txt',quote=F,row.names=F)
 
 #EXCLUDE APOE due to its leverage
 write.table(alzheimers_exp,'sumstats_reformatted/102_PGCALZ2sumstatsExcluding23andMeNOAPOE.txt',quote=F,row.names=F)


#103	Neuroticism
#104	Impulsivity

#105	Reaction time
reaction_dat <- fread('sumstats/Davies_NC_2018/Davies2018_UKB_RT_summary_results_29052018.txt',data.table=F)
reaction_dat$SNP <- reaction_dat$MarkerName
reaction_dat$CHROM <- reaction_dat$CHR
reaction_dat$POS <- reaction_dat$BP
reaction_dat$A1 <- reaction_dat$Effect_allele
reaction_dat$A2 <- reaction_dat$Other_allele
reaction_dat$BETA <- reaction_dat$Beta
reaction_dat$SE <- abs(reaction_dat$BETA/qnorm(reaction_dat$P/2,lower.tail=F))
reaction_dat$samplesize <- 330069
reaction_dat$phenotype_col <- "Reaction.Time"
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 reaction_dat2 <- merge(reaction_dat,AFs,by="SNP",suffixes=c("","_ldref"))
 betaval <- try(apply(reaction_dat2,1,flipcheck),silent=TRUE) 
 reaction_dat2$MAF <- betaval
 
   reaction_exp <- subset(reaction_dat2, , select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(reaction_exp,'sumstats_reformatted/105_Davies2018_UKB_RT_summary_results_29052018.txt',quote=F,row.names=F)


 



Davies_NC_2018/Davies2018_UKB_VNR_summary_results_29052018.txt 
#106	Verbal Numerical Reasoning
vnr_dat <- fread('sumstats/Davies_NC_2018/Davies2018_UKB_VNR_summary_results_29052018.txt',data.table=F)
vnr_dat$SNP <- vnr_dat$MarkerName
vnr_dat$CHROM <- vnr_dat$CHR
vnr_dat$POS <- vnr_dat$BP
vnr_dat$A1 <- toupper(vnr_dat$Effect_allele)
vnr_dat$A2 <- toupper(vnr_dat$Other_allele)
vnr_dat$samplesize <- 168033
vnr_dat$phenotype_col <- "VerbNumReason"

 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 vnr_dat2 <- merge(vnr_dat,AFs,by="SNP",suffixes=c("","_ldref"))
 betaval <- try(apply(vnr_dat2,1,flipcheck),silent=TRUE) 
 vnr_dat2$MAF <- betaval
 vnr_dat2$BETA = vnr_dat2$Z  / sqrt(2*vnr_dat2$samplesize*vnr_dat2$MAF*(1-vnr_dat2$MAF))
 vnr_dat2$SE <- vnr_dat2$BETA / vnr_dat2$Z
 
   vnr_exp <- subset(vnr_dat2, , select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(vnr_exp,'sumstats_reformatted/106_Davies2018_UKB_VNR_summary_results_29052018.txt',quote=F,row.names=F)



#107	Neuroticism 2
 neuroticismpost_dat <- fread('sumstats/sumstats_neuroticism_ctg_format.txt.gz',data.table=F)
 neuroticismpost_dat$SNP <- neuroticismpost_dat$RSID
 neuroticismpost_dat$CHROM <- neuroticismpost_dat$CHR
 neuroticismpost_dat$MAF <- neuroticismpost_dat$EAF_UKB
 neuroticismpost_dat$samplesize <- neuroticismpost_dat$N
 neuroticismpost_dat$phenotype_col <- "Neuroticism.2018"
 neuroticismpost_dat$BETA = neuroticismpost_dat$Z  / sqrt(2*neuroticismpost_dat$samplesize*neuroticismpost_dat$MAF*(1-neuroticismpost_dat$MAF))
 neuroticismpost_dat$SE <- neuroticismpost_dat$BETA / neuroticismpost_dat$Z
 
  
 neuroticismpost_exp <- subset(neuroticismpost_dat, , select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(neuroticismpost_exp,'sumstats_reformatted/107_sumstats_neuroticism_ctg_format.txt',quote=F,row.names=F)

 #108 Pernicious Anemia
 LC_ALL=C join <(zcat sumstats/pernicious_anemia_Laisketal2021_sumstats.gz | sed 's/rs_number/SNP/g' | LC_ALL=C sort -k1b,1  )  <(cat /mnt/sdb/genetics/parkinsonsmr/phase3.locations2 | LC_ALL=C sort -k1b,1) | sort -g -k 11 >  sumstats/pernicious_anemia_Laisketal2021_sumstats.gz2

 anemia_dat <- fread('sumstats/pernicious_anemia_Laisketal2021_sumstats.gz2',data.table=F)

 #anemia_dat$CHROM <- anemia_dat$CHR
 anemia_dat$MAF <- anemia_dat$eaf
 anemia_dat$samplesize <- anemia_dat$n_samples
 anemia_dat$phenotype_col <- "Pernicious.Anemia"
 anemia_dat$BETA = anemia_dat$beta
 anemia_dat$SE <- anemia_dat$se
  anemia_dat$A1 = anemia_dat$reference_allele
 anemia_dat$A2<- anemia_dat$other_allele
 anemia_dat$P <- anemia_dat[,"p-value"]
 anemia_dat$CHROM <- anemia_dat$CHR
 anemia_dat$POS <- anemia_dat$BP
  
 anemia_exp <- subset(anemia_dat, , select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col")) #
 write.table(anemia_exp,'sumstats_reformatted/108_pernicious_anemia_Laisketal2021_sumstats.txt',quote=F,row.names=F)

#110 addisons disease

 addisons_dat <- fread('zcat sumstats/GCST90011871_buildGRCh37.tsv.gz',data.table=F)
 addisons_dat$samplesize <-  1223 	+ 4097 
  #Need a column for sample size. If there is no column, just write in the overall N
 addisons_dat$phenotype_col <- "Addisons" #Make a dummy column just giving the trait name
 
 addisons_dat$SNP <- addisons_dat$variant_id
 addisons_dat$A1 <- addisons_dat$effect_allele
 addisons_dat$A2 <- addisons_dat$other_allele
 addisons_dat$BETA <- log(addisons_dat$odds_ratio)
 addisons_dat$SE<- addisons_dat$standard_error
 addisons_dat$P <- addisons_dat$p_value 
 addisons_dat$CHROM <- addisons_dat$chromosome
 addisons_dat$POS <- addisons_dat$base_pair_location
 addisons_dat$MAF <- addisons_dat$effect_allele_frequency
 

 addisons_exp <- subset(addisons_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(addisons_exp,'sumstats_reformatted/110_GCST90011871_buildGRCh37.tsv.gz.txt',quote=F,row.names=F)

#111 vitiligo
zcat GWAS123chr*cmh.txt.gz | sort -g -k 1 | awk '{if (NR==1 ||$1!="CHR")print}' | sed 's/RS/rs/g' > vitiligo.txt

 vitiligo_dat <- fread('sumstats/vitiligo.txt',data.table=F)
 vitiligo_dat$samplesize <-   4680 +	 39586 

  #Need a column for sample size. If there is no column, just write in the overall N
 vitiligo_dat$phenotype_col <- "Vitiligo" #Make a dummy column just giving the trait name
 
 vitiligo_dat$SNP <- vitiligo_dat$SNP
 vitiligo_dat$A1 <- vitiligo_dat$A1
 vitiligo_dat$A2 <- vitiligo_dat$A2
 vitiligo_dat$BETA <- log(vitiligo_dat$ORX)
 vitiligo_dat$SE<- vitiligo_dat$SE
 vitiligo_dat$P <- vitiligo_dat$P
 vitiligo_dat$CHROM <- vitiligo_dat$CHR
 vitiligo_dat$POS <- vitiligo_dat$BP
 vitiligo_dat$MAF <- vitiligo_dat$MAF
 

 vitiligo_exp <- subset(vitiligo_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(vitiligo_exp,'sumstats_reformatted/111_vitiligo.txt',quote=F,row.names=F)


#129 Systemic sclerosis
 systemicscler1_dat <- fread('sumstats/Lopez-Isac_prePMID_META_GWAS_SSc.meta.txt',data.table=F)
 systemicscler1_dat$samplesize <-   9095 +	17584


  #Need a column for sample size. If there is no column, just write in the overall N
 systemicscler1_dat$phenotype_col <- "systemicscler1" #Make a dummy column just giving the trait name
 
 systemicscler1_dat$SNP <- systemicscler1_dat$SNP
 systemicscler1_dat$A1 <- systemicscler1_dat$A1
 systemicscler1_dat$A2 <- systemicscler1_dat$A2
 systemicscler1_dat$BETA <- log(systemicscler1_dat$OR)

 systemicscler1_dat$P <- systemicscler1_dat$P
 systemicscler1_dat$SE <- abs(systemicscler1_dat$BETA/qnorm(systemicscler1_dat$P/2))
  
 systemicscler1_dat$CHROM <- systemicscler1_dat$CHR
 systemicscler1_dat$POS <- systemicscler1_dat$BP
 systemicscler1_dat$MAF <- systemicscler1_dat$MAF
 
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 systemicscler_dat <- merge(systemicscler1_dat,AFs,by="SNP",suffixes=c("","_ldref"))
 systemicscler_dat[which(systemicscler_dat$A1 == systemicscler_dat$A1_ldref),]$A2 <- systemicscler_dat[which(systemicscler_dat$A1 == systemicscler_dat$A1_ldref),]$A2_ldref
 systemicscler_dat[which(systemicscler_dat$A1 == systemicscler_dat$A2_ldref),]$A2 <- systemicscler_dat[which(systemicscler_dat$A1 == systemicscler_dat$A2_ldref),]$A1_ldref
 #most are assignable because the strand matches, so it is fine..
 
 betaval <- try(apply(systemicscler_dat,1,flipcheck),silent=TRUE)
 systemicscler_dat$MAF <- betaval
  
   
 systemicscler_exp <- subset(systemicscler_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(systemicscler_exp,'sumstats_reformatted/129_Lopez-Isac_prePMID_META_GWAS_SSc.meta.txt',quote=F,row.names=F)

#133 Primary scleorsing cholangitis 

awk 'NR>=12{print}' ipscsg2016.result.combined.full.with_header.txt | sed 's/#//g' > ipscsg2016.result.combined.full.with_header.txt2

 psc_dat <- fread('sumstats/ipscsg2016.result.combined.full.with_header.txt2',data.table=F)
 psc_dat$samplesize <-   4796	+ 19955 



  #Need a column for sample size. If there is no column, just write in the overall N
 psc_dat$phenotype_col <- "Primary.sclerosing.cholangitis" #Make a dummy column just giving the trait name
 
 psc_dat$SNP <- psc_dat$SNP
 psc_dat$A1 <- psc_dat$allele_1 #corrected on 4/13
 psc_dat$A2 <- psc_dat$allele_0
 psc_dat$BETA <- log(psc_dat$or)

 psc_dat$P <- psc_dat$p
 psc_dat$SE <- psc_dat$se
  
 psc_dat$CHROM <- psc_dat$chr
 psc_dat$POS <- psc_dat$pos
 psc_dat$MAF <- psc_dat$freq_1
 
    
 psc_exp <- subset(psc_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(psc_exp,'sumstats_reformatted/133_ipscsg2016.result.combined.full.with_header.txt2',quote=F,row.names=F)

#134 atopic dermatitis

 

 eczema_dat <- fread('sumstats/EAGLE_AD_no23andme_results_29072015.txt',data.table=F)
 eczema_dat$samplesize <-   21399 +	95464
 



  #Need a column for sample size. If there is no column, just write in the overall N
 eczema_dat$phenotype_col <- "Atopic.Dermatitis" #Make a dummy column just giving the trait name
 
 eczema_dat$SNP <- eczema_dat$rsID
 eczema_dat$A1 <- eczema_dat$reference_allele
 eczema_dat$A2 <- eczema_dat$other_allele
 eczema_dat$BETA <- eczema_dat$beta

 eczema_dat$P <- eczema_dat$p.value
 eczema_dat$SE <- eczema_dat$se
  
 eczema_dat$CHROM <- eczema_dat$chromosome
 eczema_dat$POS <- eczema_dat$position
 eczema_dat$MAF <- eczema_dat$eaf
 
    
 eczema_exp <- subset(eczema_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(eczema_exp,'sumstats_reformatted/134_EAGLE_AD_no23andme_results_29072015.txt',quote=F,row.names=F)

 
 #125 juvenile arthiritis
 
 
 juvenarthritis_dat <- fread('sumstats/GCST90010715_buildGRCh37.tsv',data.table=F)
 juvenarthritis_dat$samplesize <-    3305 +	 9196 

 juvenarthritis_dat$phenotype_col <- "Juvenile.Arthritis" #Make a dummy column just giving the trait name
 
 juvenarthritis_dat$SNP <- juvenarthritis_dat$variant_id
 juvenarthritis_dat$A1 <- juvenarthritis_dat$alleleB
 juvenarthritis_dat$A2 <- juvenarthritis_dat$alleleA #FIXED
 juvenarthritis_dat$BETA <- log(juvenarthritis_dat$all_OR)

 juvenarthritis_dat$P <- juvenarthritis_dat$p_value
 juvenarthritis_dat$SE <- abs( juvenarthritis_dat$BETA/qnorm( juvenarthritis_dat$P/2))
  
 juvenarthritis_dat$CHROM <- juvenarthritis_dat$chromosome
 juvenarthritis_dat$POS <- juvenarthritis_dat$position
 juvenarthritis_dat$MAF <- juvenarthritis_dat$all_maf
 
    
 juvenarthritis_exp <- subset(juvenarthritis_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(juvenarthritis_exp,'sumstats_reformatted/125_GCST90010715_buildGRCh37.tsv',quote=F,row.names=F)

#112 multiple sclerosis..

 ms_dat <- fread('sumstats/discovery_metav3.0.meta.gz',data.table=F)
 write.table(ms_dat,f
 
 ms_dat$BETA <- log(ms_dat$OR)
 ms_dat$SE <- abs( ms_dat$BETA/qnorm( ms_dat$P/2))
 ms_dat$samplesize <- 14802 +	26703
 ms_dat$phenotype_col <- "Multiple.Sclerosis" #Make a dummy column just giving the trait name
 
 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
  ms2_dat <- merge(ms_dat,AFs,by="SNP",suffixes=c("","_ldref"))
  ms2_dat$BETA <- log(ms2_dat$OR)
 ms2_dat$SE <- abs( ms2_dat$BETA/qnorm( ms2_dat$P/2))
  
   betaval <- try(apply(ms2_dat,1,flipcheck),silent=TRUE)
   ms2_dat$MAF <- betaval
 
     
 ms2_exp <- subset(ms2_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ms2_exp,'sumstats_reformatted/112_discovery_metav3.0.meta.gz',quote=F,row.names=F)

#132 ankylosing spondylitis - must be imputed
 igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv
 igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv
 
 
for chr in {1..22}
do
 cat sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv | grep -v '?' |  awk -v chr=$chr 'BEGIN{OFS="\t"}{if (NR == 1) {$1= "chr";$2="bp";$3="MarkerName";$5="Allele1";$4="Allele2";$7="beta";$6="p"}; if(NR==1 || $1 == chr) print $1,$2,$3, $5,$4,$7,$6}'  > sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_"$chr"
done

for chr in {1..22} #remaining
do
/mnt/sdb/genetics/ardissvenv/ssimp_software-master/ssimp --gwas sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_"$chr"_imputed_chr
done

cat sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_*_imputed_chr | awk '{if (NR==1 || $1 !="chr") print}' > sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_imputed_allchr

R
f.z2b <- function(z, af, n)
{
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
    se.b <- 1/sqrt(2* af * (1-af) * n)
    b <- z * se.b
    return(c(b, se.b))
}
library(data.table)
output <- fread("sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_imputed_allchr",data.table=F) ## reading ssimp output
N.max <- 15283 ## maximum sample size in the study
output$N_imp <- output$r2.pred * N.max ## calculating the effective sample size


## applying z-to-b function from above
output.b.se <- lapply(1:nrow(output), function(x) 
  f.z2b(output$z_imp[x], output$maf[x], output$N_imp[x])
)

## assigning b_imp and se_imp to dataframe
output$b_imp <- sapply(output.b.se, function(x) x[1])
output$se_imp <- sapply(output.b.se, function(x) x[2])

## store out extended output file
write.table(output, "sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_imputed_allchr.betas",quote=F,row.names=F)



il6_dat <- subset(output,N_imp > 10461.366)
 # fread('sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_imputed_allchr.betas',data.table=F)

il6_dat$CHROM <- il6_dat$chr
il6_dat$POS <- il6_dat$pos
il6_dat$A1 <- il6_dat$Allele1
il6_dat$A2 <- il6_dat$Allele2
il6_dat$P <- il6_dat$P.imp
il6_dat$BETA <- il6_dat$b_imp
il6_dat$SE <- il6_dat$se_imp
il6_dat$samplesize <- 15283
il6_dat$phenotype_col <- "il6"
il6_dat$MAF <- il6_dat$maf

 # AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 # AFs$MAF <-as.numeric(AFs$MAF)
 # il6_dat2 <- merge(il6_dat,AFs,by="SNP",suffixes=c("","_ldref"))
 # betaval <- try(apply(il6_dat2,1,flipcheck),silent=TRUE)
 
 # il6_dat2$MAF <- betaval
 

 il6_exp <- subset(il6_dat,  select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(il6_exp,'sumstats_reformatted/87_igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv',quote=F,row.names=F)


 #135 cross allergy
 
 
 allerg_dat <- fread('zcat sumstats/SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz',data.table=F)
 
 allerg_dat$samplesize <-   180129 	+ 180709 


 allerg_dat$phenotype_col <- "Cross.Allergy2" #Make a dummy column just giving the trait name
 
 allerg_dat$SNP <- allerg_dat$RS_ID
 allerg_dat$A1 <- allerg_dat$EFFECT_ALLELE
 allerg_dat$A2 <- allerg_dat$OTHER_ALLELE #check if this is right..
# allerg_dat$BETA <- log(allerg_dat$all_OR)

 allerg_dat$P <- allerg_dat$PVALUE
 #allerg_dat$SE <- abs( allerg_dat$BETA/qnorm( allerg_dat$P/2))
  
 allerg_dat$CHROM <- allerg_dat$CHR
 allerg_dat$POS <- allerg_dat$BP
 allerg_dat$MAF <- allerg_dat[,"1000G_ALLELE_FREQ"]
 
 
  allerg_exp <- subset(allerg_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(allerg_exp,'sumstats_reformatted/135_SHARE-without23andMe.LDSCORE-GC.SE-META.v0.gz.txt',quote=F,row.names=F)

 #139 asthma 3
 #this is a bolt analysis so effect would have to be rescaled based on prevalence to the OR scale.
 asthmav3_dat <- fread('zcat sumstats/ADULT1_ADULT2_ONSET_ASTHMA.20180716.allchr.assoc.GC.gz',data.table=F)
 
 asthmav3_dat$samplesize <-    26582 +	 300671 

 asthmav3_dat$phenotype_col <- "Asthma.v3" #Make a dummy column just giving the trait name
 
 #asthmav3_dat$SNP <- asthmav3_dat$variant_id
 asthmav3_dat$A1 <- asthmav3_dat$ALLELE1
 asthmav3_dat$A2 <- asthmav3_dat$ALLELE0 #check if this is right..

 asthmav3_dat$P <- asthmav3_dat$P_BOLT_LMM_INF
 
 asthmav3_dat$CHROM <- asthmav3_dat$CHR
 asthmav3_dat$POS <- asthmav3_dat$BP
 asthmav3_dat$MAF <- asthmav3_dat$A1FREQ
 
 
 asthmav3_exp <- subset(asthmav3_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(asthmav3_exp,'sumstats_reformatted/139_ADULT1_ADULT2_ONSET_ASTHMA.20180716.allchr.assoc.GC.gz.txt',quote=F,row.names=F)

#140 IBD version 2 -- no for the last time, this is unusable due to only having leading associations outputted for hte transethnic..
 ibdv2_dat <- fread('zcat sumstats/iibdgc-trans-ancestry-filtered-summary-stats.tgz',data.table=F)
 
 ibdv2_dat$samplesize <-    42950 +	53536


 ibdv2_dat$phenotype_col <- "Asthma.v3" #Make a dummy column just giving the trait name
 
 #ibdv2_dat$SNP <- ibdv2_dat$variant_id
 ibdv2_dat$A1 <- ibdv2_dat$ALLELE1
 ibdv2_dat$A2 <- ibdv2_dat$ALLELE0 #check if this is right..

 ibdv2_dat$P <- ibdv2_dat$P_BOLT_LMM_INF
 
 ibdv2_dat$CHROM <- ibdv2_dat$CHR
 ibdv2_dat$POS <- ibdv2_dat$BP
 ibdv2_dat$MAF <- ibdv2_dat$A1FREQ
 
 
 ibdv2_exp <- subset(ibdv2_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(ibdv2_exp,'sumstats_reformatted/140_iibdgc-trans-ancestry-filtered-summary-stats.tgz.txt',quote=F,row.names=F)

#141 eczema 2 https://www.ebi.ac.uk/gwas/studies/GCST90027161
 eczemav2_dat <- fread('zcat sumstats/34454985-GCST90027161-EFO_0000274-Build38.f.tsv.gz',data.table=F)
 
 eczemav2_dat$samplesize <-     22474 + 774187 

4/((1/22474 + (1/774187 ))


 eczemav2_dat$phenotype_col <- "AtopicDermatitis.v2" #Make a dummy column just giving the trait name
 


 eczemav2_dat$SNP <- eczemav2_dat$variant_id
 eczemav2_dat$A1 <- eczemav2_dat$effect_allele
 eczemav2_dat$A2 <- eczemav2_dat$other_allele #check if this is right..

 eczemav2_dat$P <- eczemav2_dat$p_value
 
 eczemav2_dat$CHROM <- eczemav2_dat$chromosome
 eczemav2_dat$POS <- eczemav2_dat$base_pair_location
 eczemav2_dat$MAF <- eczemav2_dat$effect_allele_frequency
  eczemav2_dat$BETA <- eczemav2_dat$beta
 eczemav2_dat$SE <- eczemav2_dat$standard_error
 
 
   AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
  AFs <- AFs[,1:3]
 eczemav2_dat2 <- merge(eczemav2_dat,AFs,by="SNP",suffixes=c("_original",""))
 
 
 eczemav2_exp <- subset(eczemav2_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(eczemav2_exp,'sumstats_reformatted/141_GCST90027161-EFO_0000274-Build38.f.tsv.gz.txt',quote=F,row.names=F)


 

IFS=$'\n'
for files in $(cat /mnt/ukbb/adam/ptsd/sumner_mendelian_randomization/traitlist.txt)
do
 trait=$(echo $files | awk '{print $1}'
 mv $trait 