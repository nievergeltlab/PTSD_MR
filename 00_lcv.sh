#Make LDSC files
# cp -r /mnt/ukbb/adam/ptsd/eur_w_ld_chr .
# for files in $(ls | grep gz) ; do gzip -d $files; done
#cat *.ldscore | sort -g -k 1 | awk '{if (NR==1 ||$1 != "CHR") print}' | LC_ALL=C sort -g -k1,1 -k3,3 > all.scored

library(data.table)
library(plyr)

#The faster way to do this would be al imited set of array operations, where you just compare the four different cases manually then combine over all

#Flip Z score 
 flipcheck <- function(x)
 {
  betaval=as.numeric(x["Z.x"])
 # print(betaval)
  snp=c(x["A1.x"],x["A2.x"])
  pivot=c(x["A1.y"],x["A2.y"]) #this will be the reference
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
    betaval=-betaval
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
    betaval=-betaval
  } else {
  betaval=NA 
  }
  return(betaval)
  }
  
# LCV example script written by Katie Siewert 
source('LCV-master/R/MomentFunctions.R')
source('LCV-master/R/RunLCV.R')

#Load LD scores 
 d3=fread('LCV-master/eur_w_ld_chr/all.scored',data.table=F)

#List of files to check with PTSD
 fileslist=scan(what=character())
81_GCST90016564_buildGRCh37.tsv


 59_C_reactive_protein.imp.gz.txt


107_sumstats_neuroticism_ctg_format.txt


77_eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.txt

81_GCST90016564_buildGRCh37.tsv
52_UKBB.asthma.assoc
54_GCST90018907_buildGRCh37.tsv
94_IL6_summary_results_discovery_GWAS_MA_CHARGE.txt
18_BCX2_WBC_EA_GWAMA.out.gz.fuma.txt
55_lupus.ea.imputed.allchr.out
57_RA_GWASmeta_European_v2.txt
21_BCX2_MON_EA_GWAMA.out.gz.fuma.txt
44_EUR.IBD.gwas_info03_filtered.assoc
43_EUR.CD.gwas_info03_filtered.assoc
19_BCX2_NEU_EA_GWAMA.out.gz.fuma.txt

 #fileslist=fileslist1[which(fileslist1 %in% analyze)]  #Change this if you add more
 
#Loop over PTSD (if needed), then over trait
for (ptsdfile in c( "sumstats_reformatted/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt" ))#  ,"sumstats_reformatted/77_eur_ptsd_pcs_v4_aug3_2021_noukbb.fuma.txt"   "sumstats_reformatted/95_eur_ptsd_pcs_v4_aug3_2021_nomvpukbb.fuma.txt","sumstats_reformatted/78_eur_ptsd_pcs_v4_aug3_2021_nomvp.fuma.txt" ))#,)) #
{
   ptsdU_dat <- fread(ptsdfile,data.table=F)
   #cd40 should be added, as sensitivity analysis.chr20:44746906-44758384 
   ptsdU_dat <- subset(ptsdU_dat, !(CHROM == 6 & POS >= 28477797 - 3000000 & POS <= 33448354 + 3000000) & 
                                  !(CHROM == 17 & POS >= 40928986 - 3000000 & POS <= 42139672 + 3000000)) #&
                               #  !(CHROM == 20 & POS >= 44746906 - 3000000 & POS <= 44758384  + 3000000)) #remove HLA, MAPT inversion, CD40 gene |/- 3MB
   ptsdU_dat <- subset(ptsdU_dat, !((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A")  | (A1 == "C" & A2 == "G")  | (A1 == "G" & A2 == "C"))) #Remove ambiguous variants
   
   ptsdset=ptsdU_dat$phenotype_col[1]

  #Calculate Zs
   d1 = ptsdU_dat
   d1$Z <- d1$BETA/d1$SE
   d1 <- d1[!is.na(d1$Z),]
  #Recommended MAF filter
   d1 <- subset(d1,MAF >=0.05 & MAF <= 0.95)

  #Merge with LD scores
  m = merge(d3,d1,by="SNP")
       
  #load disease files. Disease name is in phenotype file itself
  for (diseasefile in fileslist)
  {
   print (diseasefile)
   #Load dataset - this wont work on datasets not updated to hg19! check which ones are not.
   disease_data <- fread(paste('sumstats_reformatted/',diseasefile,sep=''),data.table=F)
   disease_data <- subset(disease_data, !(CHROM == 6 & POS >= 28477797 - 3000000 & POS <= 33448354 + 3000000) & 
                                        !(CHROM == 17 & POS >= 40928986 - 3000000 & POS <= 42139672 + 3000000) )#&
                                      #  !(CHROM == 20 & POS >= 44746906 - 3000000 & POS <= 44758384  + 3000000) ) #remove HLA, MAPT inversion, +/- 3MB
   disease_data <- subset(disease_data, !((A1 == "A" & A2 == "T") | (A1 == "T" & A2 == "A")  | (A1 == "C" & A2 == "G")  | (A1 == "G" & A2 == "C"))) #Remove ambiguous
   
   disease=disease_data$phenotype_col[1]
   exposurefname=ptsdset
   outcomefname=disease
   
#Load trait 2 data and calculate Zs
 d2 = disease_data
 d2$Z <- d2$BETA/d2$SE
 d2 <- d2[!is.na(d2$Z),]
#Recommended MAF filter
 d2 <- subset(d2,MAF >=0.05 & MAF <= 0.95)


 #Merge data
 data = merge(m,d2,by="SNP")

 #Sort by position 
 data = data[order(data[,"CHR"],data[,"BP"]),]

 #Flip sign of one z-score if opposite alleles-shouldn't occur with UKB data
 #If not using munged data, will have to check that alleles match-not just whether they're opposite A1/A2
  mismatch = which(data$A1.x!=data$A1.y | data$A2.x!=data$A2.y) #mismatch detection for both A1 and A2, insures that tri-allelics are not errenously included

 #Adam: Strand correction. Only apply strand correction over mismatches, to save time
 if(length(mismatch) > 0)
 {
  bvs <-  try(apply(data[mismatch,],1,flipcheck),silent=TRUE)
  data[mismatch,]$Z.x <- bvs
  print(length(bvs)) #how many were attempted
  print(table(!is.na(bvs))) #how many were successful? true is success
 }

 # apply(data[mismatch[1:4],],1,flipcheck) sample code

 #Some will be incompatible, only keep numeric, non NA entries
 print(dim(data)[1])
  data <- subset(data, is.numeric(Z.x) & is.numeric(Z.y) & !is.na(Z.x))
 print(dim(data)[1])
 #Run LCV-need (modified so that the helper script is already called, no path needed)

 LCV = RunLCV(data$L2,data$Z.x,data$Z.y)
 capture.output(sprintf("Estimated posterior gcp=%.2f(%.2f), log10(p)=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, log(LCV$pval.gcpzero.2tailed)/log(10), LCV$rho.est, LCV$rho.err),
 file=paste('results/abridged_LCV_',exposurefname,"_",outcomefname,".txt",sep=""))


 write.table(t(LCV),file =paste('results/LCV_',exposurefname,"_",outcomefname,".txt",sep=""),row.names=F)
 
 #LCV <- NULL
   }
   }
    