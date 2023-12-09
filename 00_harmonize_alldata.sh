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

psoriasis_dat$SNP <- 
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

#57 Eurpean RA GWAS meta-analysis (14,361 RA cases and 43,923 conrols) [Link]
for chr in {1..22}
do
 zcat sumstats/RA_GWASmeta_European_v2.txt.gz | awk '{$1=$1; print}' | grep -v  "A T" | grep -v  "T A"  |  grep -v  "G C" | grep -v  "C G"  |  awk -v chr=$chr 'BEGIN{OFS="\t"}{beta=log($6); if (NR == 1) {$2= "chr";$3="bp";$1="MarkerName";$4="Allele1";$5="Allele2";beta="beta";$9="p"}; if(NR==1 || $2 == chr) print $1,$2,$3, $4,$5,beta,$9}'  > sumstats/RA_GWASmeta_European_v2.txt.gz.ssimp_"$chr"
done

for chr in {1..22} #remaining
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




#76 PTSD Meta analysis
 ptsd_dat <- fread('zcat ../eur_ptsd_pcs_v4_aug3_2021.fuma.gz',data.table=F)
 #ptsd_dat <- subset(ptsd_dat, !(Chromosome == 6 & Position >= 25000000 & Position <= 35000000))
 ptsd_dat$samplesize <- ptsd_dat$Weight  
 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 ptsd_dat$BETA = ptsd_dat$Zscore * sed / sqrt(2*(ptsd_dat$Weight + ptsd_dat$Zscore^2)*ptsd_dat$Freq1*(1-ptsd_dat$Freq1))
 ptsd_dat$SE <- ptsd_dat$BETA / ptsd_dat$Zscore 
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
 ptsd_dat$BETA = ptsd_dat$Zscore * sed / sqrt(2*(ptsd_dat$Weight + + ptsd_dat$Zscore^2)*ptsd_dat$Freq1*(1-ptsd_dat$Freq1))
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
for chr in {1..22}
do
 zcat sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz | awk '{$1=$1; print}' | grep -v  "A T" | grep -v  "T A"  |  grep -v  "G C" | grep -v  "C G"  |  awk -v chr=$chr 'BEGIN{OFS="\t"}{if (NR == 1) {$1= "chr";$2="bp";$3="MarkerName";$5="Allele1";$4="Allele2";$7="beta";$6="p"}; if(NR==1 || $1 == chr) print $1,$2,$3, $5,$4,$7,$6}'  > sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_"$chr"
done

for chr in {1..22} #remaining
do
/mnt/sdb/genetics/ssimp_software-master/ssimp --gwas sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_"$chr"_imputed_chr"$chr" 
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
 # celiac_dat <- fread('/mnt/ukbb/sumstats/dubois_2010_20190752_cel_efo0001060_1_gwas.sumstats.tsv.gz.ssimp_imputed_allchr.betas',data.table=F)

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




# 94 IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz
#NOTE: readme says that these are HG18 POSITIONS! 

#Liftover to hg19

LC_ALL=C join <(zcat eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk '{print $1,$2,$3}' | LC_ALL=C sort -k1b,1 ) <(zcat sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz | cut -f 1-8,11- | LC_ALL=C sort -k1b,1) > sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.lifted

for chr in {1..22}
do
 cat  sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.lifted | awk '{$1=$1; print}' | grep -v  "A T" | grep -v  "T A"  |  grep -v  "G C" | grep -v  "C G"    |  awk -v chr=$chr 'BEGIN{OFS="\t"}{if (NR == 1) {$2= "chr";$3="bp";$1="MarkerName";$4="Allele1";$5="Allele2";$7="beta";$9="p"}; if(NR==1 || $2 == chr) print $1,$2,$3, $4,$5,$7,$9}'  > sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_"$chr"
done

for chr in  {1..22} #remaining
do
/mnt/sdb/genetics/ssimp_software-master/ssimp --gwas sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_"$chr" --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out sumstats/IL6_summary_results_discovery_GWAS_MA_CHARGE.txt.gz.ssimp_"$chr"_imputed_chr"$chr" 
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
/mnt/sdb/genetics/ssimp_software-master/ssimp --gwas sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_"$chr"_imputed_chr
done


cat sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_*_imputed_chr | awk '{if (NR==1 || $1 !="chr") print}' > sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_imputed_allchr

R
f.z2b <- function(z, af, n)
{
    # z = imputed z statistics
    # af = allele frequency
    # n = sample size (effective)
    se.b <- 1/sqrt(2* af * (1-af) * (n+ z^2))
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

il6_dat <- fread("sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_imputed_allchr.betas",data.table=F)


il6_dat2 <- subset(il6_dat ,N_imp > 10461.366)
 # fread('sumstats/igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv.ssimp_imputed_allchr.betas',data.table=F)

il6_dat$CHROM <- il6_dat$chr
il6_dat$POS <- il6_dat$pos
il6_dat$A1 <- il6_dat$Allele1
il6_dat$A2 <- il6_dat$Allele2
il6_dat$P <- il6_dat$P.imp
il6_dat$BETA <- il6_dat$b_imp
il6_dat$SE <- il6_dat$se_imp
il6_dat$samplesize <- 15283
il6_dat$phenotype_col <- "ankylosing"
il6_dat$MAF <- il6_dat$maf

 # AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 # AFs$MAF <-as.numeric(AFs$MAF)
 # il6_dat2 <- merge(il6_dat,AFs,by="SNP",suffixes=c("","_ldref"))
 # betaval <- try(apply(il6_dat2,1,flipcheck),silent=TRUE)
 
 # il6_dat2$MAF <- betaval
 

 il6_exp <- subset(il6_dat,  select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(il6_exp,'sumstats_reformatted/87_igas_2103_23749187_as_efo0003898_1_ichip.sumstats.tsv',quote=F,row.names=F)




#143 myasthenia gravis

 myastenia_dat <- fread('sumstats/GCST90093061_buildGRCh38.tsv',data.table=F)
 
 myastenia_dat$samplesize <-    1873 +	36370


 myastenia_dat$phenotype_col <- "myasthenia.gravis" #Make a dummy column just giving the trait name
 
 myastenia_dat$SNP <- myastenia_dat$variant_id
 myastenia_dat$A1 <- myastenia_dat$effect_allele
 myastenia_dat$A2 <- myastenia_dat$other_allele #check if this is right..

 myastenia_dat$P <- myastenia_dat$p_value
 
 myastenia_dat$CHROM <- myastenia_dat$chromosome
 myastenia_dat$POS <- myastenia_dat$base_pair_location 
 myastenia_dat$MAF <- round(myastenia_dat$effect_allele_frequency,4)
 myastenia_dat$BETA <- myastenia_dat$beta
 myastenia_dat$SE <- myastenia_dat$standard_error
 
 
 myastenia_exp <- subset(myastenia_dat, MAF >= 0.005 & MAF <= 0.995, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(myastenia_exp,'sumstats_reformatted/143_GCST90093061_buildGRCh38.tsv',quote=F,row.names=F)

#Neuromyelitis optic

 nmo_dat <- fread('zcat sumstats/NMO.Combined.gz',data.table=F)
 
 #Need external freqs, 
 #calculate SE myself
 #capitalize A1/A2
 
nmo_dat$A1 <- toupper(nmo_dat$Allele1)
nmo_dat$A2 <- toupper(nmo_dat$Allele2)
nmo_dat$SNP <- nmo_dat$rsID
nmo_dat$CHROM <- nmo_dat$Chr
nmo_dat$POS <- nmo_dat$Pos
nmo_dat$BETA <- nmo_dat$Effect
nmo_dat$P <- nmo_dat$P
nmo_dat$phenotype_col <- "NMO"
nmo_dat$SE <- abs(nmo_dat$BETA)/ abs(qnorm(nmo_dat$P/2))


 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 nmo_dat2 <- merge(nmo_dat,AFs[,c(1,6)],by="SNP",suffixes=c("_original",""))
nmo_dat2$samplesize <- nmo_dat2$N 


 nmo_exp <- subset(nmo_dat2 , select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(nmo_exp,'sumstats_reformatted/144_NMO.Combined.gz',quote=F,row.names=F)

#145 Euosinophilic granualamatosis
 unlist_split <- function(x, ...)
  {
	toret <- unlist(strsplit(x, ...) )[1]
	return(t(toret))
  }
  
esog_dat <- fread('sumstats/all.egpa.vs.controls.txt',data.table=F)
SNPS  <- sapply(esog_dat$Variant,unlist_split,split=":")
esog_dat$SNP <- SNPS

esog_dat$A1 <- toupper(esog_dat$ALLELE1)
esog_dat$A2 <- toupper(esog_dat$ALLELE0)
esog_dat$CHROM <- esog_dat$CHR
esog_dat$POS <- esog_dat$BP
esog_dat$MAF <- esog_dat$A1FREQ
esog_dat$phenotype_col <- "eosinophilic_granulomatosis"
esog_dat$samplesize <- 676 +	6809

 esog_exp <- subset(esog_dat, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(esog_exp,'sumstats_reformatted/145_all.egpa.vs.controls.txt',quote=F,row.names=F)


#146 kawasaki
kawasaki_dat <- fread('sumstats/GCST90013537_buildGRCh37.tsv',data.table=F)

kawasaki_dat$CHRPOS <- paste(kawasaki_dat$chromosome,kawasaki_dat$base_pair_location,sep=":")

 AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 AFs$CHRPOS <- paste(AFs$CHROM,AFs$POS,sep=":")

 kawasaki_dat2 <- merge(AFs,kawasaki_dat,by="CHRPOS")

 kawasaki_dat2$A1 <- kawasaki_dat2$effect_allele
 kawasaki_dat2$A2 <- kawasaki_dat2$other_allele
 kawasaki_dat2$samplesize <- 400	+6101
 kawasaki_dat2$phenotype_col <- "kawasaki"

 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 kawasaki_dat2$BETA = kawasaki_dat2$z.meta * sed / sqrt(2*kawasaki_dat2$samplesize*kawasaki_dat2$MAF*(1-kawasaki_dat2$MAF))
 kawasaki_dat2$SE <- kawasaki_dat2$BETA / kawasaki_dat2$z.meta
 kawasaki_dat2$P <- kawasaki_dat2$p_value
 
 kawasaki_exp <- subset(kawasaki_dat2, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(kawasaki_exp,'sumstats_reformatted/146_GCST90013537_buildGRCh37.tsv',quote=F,row.names=F)

#147 sarcoidosis LS #has to be imputed


for chr in {1..22}
do
 cat sumstats/LS_dosage_HCR_summary_statistics.TAB | awk '{$1=$1; print}' | grep -v  "A T" | grep -v  "T A"  |  grep -v  "G C" | grep -v  "C G"   |  awk -v chr=$chr 'BEGIN{OFS="\t"}{beta=log($8); if (NR == 1) {$1= "chr";$3="bp";$2="MarkerName";$4="Allele1";$5="Allele2";beta="beta";$12="p"}; if(NR==1 || $1 == chr) print $1,$3,$2, $4,$5,beta,$12}'  > sumstats/LS_dosage_HCR_summary_statistics.TAB.ssimp_"$chr"
done

for chr in {1..22} #22 # {9..22} #remaining
do
/mnt/sdb/genetics/ssimp_software-master/ssimp --gwas /mnt/ukbb/sumstats/LS_dosage_HCR_summary_statistics.TAB.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out /mnt/ukbb/sumstats/LS_dosage_HCR_summary_statistics.TAB.ssimp_"$chr"_imputed_chr"$chr" 
done

cat /mnt/ukbb/sumstats/LS_dosage_HCR_summary_statistics.TAB.ssimp_*_imputed_chr* | awk '{if (NR==1 || $1 !="chr") print}' > /mnt/ukbb/sumstats/LS_dosage_HCR_summary_statistics.TAB.ssimp_allchr_ssimp

library(data.table)

 sarcoidLS_dat <- fread('/mnt/ukbb/sumstats/LS_dosage_HCR_summary_statistics.TAB.ssimp_allchr_ssimp',data.table=F)
 
 
 sarcoidLS_dat$samplesize <-   384	+ 2025


 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 sarcoidLS_dat$BETA = sarcoidLS_dat$z_imp * sed / sqrt(2*sarcoidLS_dat$samplesize*sarcoidLS_dat$maf*(1-sarcoidLS_dat$maf))
 sarcoidLS_dat$SE <- sarcoidLS_dat$BETA / sarcoidLS_dat$z_imp 

 sarcoidLS_dat$BETA <- round(sarcoidLS_dat$BETA,4)
 sarcoidLS_dat$SE <- round(sarcoidLS_dat$SE,4)
 sarcoidLS_dat$phenotype_col <- "sarcoidosis.lofgrens" #Make a dummy column just giving the trait name
 

 sarcoidLS_dat$A1 <- sarcoidLS_dat$Allele1
 sarcoidLS_dat$A2 <- sarcoidLS_dat$Allele2 #check if this is right..

 sarcoidLS_dat$P <- sarcoidLS_dat$P.imp
 
 sarcoidLS_dat$CHROM <- sarcoidLS_dat$chr
 sarcoidLS_dat$POS <- sarcoidLS_dat$pos
 sarcoidLS_dat$MAF <- round(sarcoidLS_dat$maf,4)

 
 
 sarcoidLS_exp <- subset(sarcoidLS_dat, r2.pred > 0.6, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(sarcoidLS_exp,'sumstats_reformatted/147_LS_dosage_HCR_summary_statistics.TAB.txt',quote=F,row.names=F)


#148 sarcoidosis no LS
for chr in {1..22}
do
 cat /mnt/ukbb/sumstats/nonLS_dosage_HCR_summary_statistics.TAB | awk '{$1=$1; print}' | grep -v  "A T" | grep -v  "T A"  |  grep -v  "G C" | grep -v  "C G"   |  awk -v chr=$chr 'BEGIN{OFS="\t"}{beta=log($8); if (NR == 1) {$1= "chr";$3="bp";$2="MarkerName";$4="Allele1";$5="Allele2";beta="beta";$12="p"}; if(NR==1 || $1 == chr) print $1,$3,$2, $4,$5,beta,$12}'  > /mnt/ukbb/sumstats/nonLS_dosage_HCR_summary_statistics.TAB.ssimp_"$chr"
done

for chr in {1..22} #22 # {9..22} #remaining
do
/mnt/sdb/genetics/ssimp_software-master/ssimp --gwas /mnt/ukbb/sumstats/nonLS_dosage_HCR_summary_statistics.TAB.ssimp_"$chr"  --impute.range $chr --ref ~/reference_panels/1000genomes/ALL.chr"$chr".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --impute.maf 0.01 --sample.names ~/reference_panels/1KG/integrated_call_samples_v3.20130502.ALL.panel/sample/super_pop=EUR --out /mnt/ukbb/sumstats/nonLS_dosage_HCR_summary_statistics.TAB.ssimp_"$chr"_imputed_chr"$chr" 
done

cat /mnt/ukbb/sumstats/nonLS_dosage_HCR_summary_statistics.TAB.ssimp_*_imputed_chr* | awk '{if (NR==1 || $1 !="chr") print}' > /mnt/ukbb/sumstats/nonLS_dosage_HCR_summary_statistics.TAB.ssimp_allchr_ssimp

library(data.table)

 sarcoidnonLS_dat <- fread('/mnt/ukbb/sumstats/nonLS_dosage_HCR_summary_statistics.TAB.ssimp_allchr_ssimp',data.table=F)
 
 
 sarcoidnonLS_dat$samplesize <-   568 + 2086


 sed=3.09 #Give standard error for PTSD measure- This is for CAPS. Converted to Pseudo betas
 sarcoidnonLS_dat$BETA = sarcoidnonLS_dat$z_imp * sed / sqrt(2*sarcoidnonLS_dat$samplesize*sarcoidnonLS_dat$maf*(1-sarcoidnonLS_dat$maf))
 sarcoidnonLS_dat$SE <- sarcoidnonLS_dat$BETA / sarcoidnonLS_dat$z_imp 

 sarcoidnonLS_dat$BETA <- round(sarcoidnonLS_dat$BETA,4)
 sarcoidnonLS_dat$SE <- round(sarcoidnonLS_dat$SE,4)
 sarcoidnonLS_dat$phenotype_col <- "sarcoidosis.NONlofgrens" #Make a dummy column just giving the trait name
 

 sarcoidnonLS_dat$A1 <- sarcoidnonLS_dat$Allele1
 sarcoidnonLS_dat$A2 <- sarcoidnonLS_dat$Allele2 #check if this is right..

 sarcoidnonLS_dat$P <- sarcoidnonLS_dat$P.imp
 
 sarcoidnonLS_dat$CHROM <- sarcoidnonLS_dat$chr
 sarcoidnonLS_dat$POS <- sarcoidnonLS_dat$pos
 sarcoidnonLS_dat$MAF <- round(sarcoidnonLS_dat$maf,4)

 
 
 sarcoidnonLS_exp <- subset(sarcoidnonLS_dat, r2.pred > 0.6, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(sarcoidnonLS_exp,'sumstats_reformatted/148_nonLS_dosage_HCR_summary_statistics.TAB.txt',quote=F,row.names=F)


#150 rheumatoid arthritis 2023
#Merge on this
 locations <- fread('sumstats/GCST90132222_buildGRCh37.tsv.gz',data.table=F)
 arthritis_dat <- fread('sumstats/EUR_all_auto-10-2021.txt.gz',data.table=F)
 names(arthritis_dat) <- c("variant_id_hg19","BETA","SE","P")
 
 arthritis_dat2 <- merge(arthritis_dat,locations,by="variant_id_hg19")
 arthritis_dat2$samplesize <- 97173
 arthritis_dat2$phenotype_col <- "arthritis_2023"
 arthritis_dat2$A1 <- arthritis_dat2$effect_allele
 arthritis_dat2$A2 <- arthritis_dat2$other_allele
 arthritis_dat2$CHROM <- arthritis_dat2$chromosome
 arthritis_dat2$POS <- arthritis_dat2$base_pair_location
  arthritis_dat2$SNP <- arthritis_dat2$variant_id
   
  AFs <- fread('eurmaf005_ALL.allchr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',data.table=F)
 AFs$MAF <-as.numeric(AFs$MAF)
 
 arthritis_dat3 <- merge(arthritis_dat2,AFs,by="SNP",suffixes=c("","_afs"))
 arthritis_dat3$P <- as.numeric(arthritis_dat3$P)
 
 arthritis_exp <- subset(arthritis_dat3, MAF >= 0.01 & MAF <=0.99, select=c("SNP","CHROM","POS","A1","A2","MAF","BETA","SE","P","samplesize","phenotype_col"))
 write.table(arthritis_exp,'sumstats_reformatted/150_EUR_all_auto-10-2021.txt.gz.txt',quote=F,row.names=F)

