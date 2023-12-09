#Perform genetic correlations between a trait and all existing PGC psychiatric traits (as of july 2021)

#need to make a file where the possible alleles are given - will crash otherwise, because sometimes (rare) the same rs has different coded alleles 
#Basically just make a file consisting of the SNP, A1, A2 from your main dataset where you want to see what traits are genetically correlated to it.


LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/software/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat sumstats_reformatted/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt | awk '{if(NR==1) {$4="A1"; $5="A2"}; print $1,toupper($4),toupper($5)}' | LC_ALL=C sort -k1b,1 | sed 's/MarkerName/SNP/g')  > ldscrg/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt.alleles

#Reformat main summary statistics

#HLA and CHR 17 inversion will be excluded from rg estimates

                                  
LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/software/w_hm3.snplist.sorted | LC_ALL=C sort -k1b,1 ) <(cat sumstats_reformatted/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt |  awk {'if (NR==1 ||  !($2 == 6 && $3 >= 28477797 - 3000000 && $2 <= 33448354 + 3000000) && !($2 == 17 && POS >= 40928986 - 3000000 && $2 <= 42139672 + 3000000)) print}'  | sed 's/MarkerName/SNP/g'  | LC_ALL=C sort -u -k1b,1 ) > ldscrg/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt.premunge
#mhc 

python2  /mnt/ukbb/software/ldsc/munge_sumstats.py --chunksize 500000  --sumstats ldscrg/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt.premunge  --N-col samplesize --out ldscrg/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt.munge.gz
 
IFS=$'\n'
for dataset in  $(cat prevalence_data.txt)
do
 echo $file
 file=$(echo $dataset  | awk '{print $1}')
 #if [ ! -f ldscrg/"$file".munge.gz.sumstats.gz ]
 #then
  LC_ALL=C join -1 1 -2 1  <(awk '{print $1}' /mnt/ukbb/software/w_hm3.noMHC.snplist | LC_ALL=C sort -k1b,1 ) <(cat sumstats_reformatted/$file   |  LC_ALL=C sort -u -k1b,1 ) > ldscrg/"$file".premunge
  python2 /mnt/ukbb/software/ldsc/munge_sumstats.py --chunksize 500000  --sumstats ldscrg/"$file".premunge --N-col samplesize --merge-alleles ldscrg/76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt.alleles --out ldscrg/"$file".munge.gz
# fi
done

#Estimate heritabilities

for dataset in $(cat prevalence_data.txt) 
do
 echo $file
  file=$(echo $dataset  | awk '{print $1}')
  samp_prev=$(echo $dataset  | awk '{print $2}')
  pop_prev=$(echo $dataset  | awk '{print $3}')
  
 if [ $samp_prev != "NA" ]
 then
   echo "Using sample prevaelnce for $dataset"
     python2 /mnt/ukbb/software/ldsc/ldsc.py \
  --h2 ldscrg/"$file".munge.gz.sumstats.gz \
  --ref-ld-chr /mnt/ukbb/software/ldsc/eur_w_ld_chr/ \
  --w-ld-chr  /mnt/ukbb/software/ldsc/eur_w_ld_chr/ \
  --samp-prev $samp_prev --pop-prev $pop_prev \
  --out h2s/h2_$file
  else 
 
 # if [ ! -f h2s/h2_$file.log ]
 #  then
    python2 /mnt/ukbb/software/ldsc/ldsc.py \
   --h2 ldscrg/"$file".munge.gz.sumstats.gz \
   --ref-ld-chr /mnt/ukbb/software/ldsc/eur_w_ld_chr/ \
   --w-ld-chr  /mnt/ukbb/software/ldsc/eur_w_ld_chr/ \
   --out h2s/h2_$file
 #fi
 fi
 done

for file in $(ls h2s)
do
 h2=$(grep "scale h2" h2s/$file) 
 echo $file $h2 - >> h2s.txt
 done
 
  
for ptsd in 76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt 
do
for dataset in $(cat prevalence_data.txt) 
do
 echo $file
  file=$(echo $dataset  | awk '{print $1}')
  samp_prev=$(echo $dataset  | awk '{print $2}')
  pop_prev=$(echo $dataset  | awk '{print $3}')
  
# if [ $samp_prev != "NA" ]
 #then
  #prev_flag= --samp-prev $samp_prev --pop-prev $pop_prev
 # fi
 # if [ ! -f rgs/rg_"$ptsd"_"$file".log ]
#  then
   python2 /mnt/ukbb/software/ldsc/ldsc.py \
  --rg  ldscrg/"$ptsd".munge.gz.sumstats.gz,ldscrg/"$file".munge.gz.sumstats.gz \
  --ref-ld-chr /mnt/ukbb/software/ldsc/eur_w_ld_chr/ \
  --w-ld-chr  /mnt/ukbb/software/ldsc/eur_w_ld_chr/  \
  --out rgs/rg_"$ptsd"_$file
 #fi
 done
done

#Send rg outputs to a single file

for ptsd in 76_eur_ptsd_pcs_v4_aug3_2021.fuma.txt 
do

for analysis in $(ls rgs |grep -v unused | grep -v old) 
do
  echo $analysis
  grep "Summary of Genetic Correlation Results" -A2  rgs/$analysis | tail -n1 >> "$ptsd"_rgs_june7_2021.txt
 done
 done
 

#sample code for h2
 python2 /mnt/ukbb/adam/tinnitus_gwas/ldsc-master/ldsc.py \
 --rg ldscrg/81_GCST90016564_buildGRCh37.tsv.munge.gz.sumstats.gz,ldscrg/44_EUR.IBD.gwas_info03_filtered.assoc.munge.gz.sumstats.gz \
 --ref-ld-chr /mnt/ukbb/royce/eur_w_ld_chr/ \
 --w-ld-chr  /mnt/ukbb/royce/eur_w_ld_chr/ \
 --out rgs/ibd_to_ibd
 

 
 
 
  