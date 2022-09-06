## Covert plink text files to VCF
## Note: the default behaviour of Plink is to set the major allele to A2.
## Note: default in Plink --horse = --chr-set 31 no-xy no-mt
##       if you give --chr-set n, and your bam file

suffix="70K_Equcab2_Haflinger_Bellone"
#plink --file "$suffix" --horse --output-chr chrM --not-chr y --recode vcf --out "$suffix"
plink --file "$suffix" --horse --make-bed --out "$suffix"
#plink --file "$suffix" --chr-set 36 --make-bed --out "$suffix"

# fix the bim file where both alleles are abscent
cp "$suffix".bim "$suffix".bim_temp
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $3;next;}{if($5 == "0" && $6 == "0" && a[$2])print $1,$2,$3,$4,a[$2];else print $0}' $arr70k "$suffix".bim_temp > "$suffix".bim ## #else if($6 == "0")print $2,$5,$6,$8,$9; ==> this line should be useless now

# generate the SNP list to be updated
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{if(a[$2])print $0,a[$2]}' $arr70k  "$suffix".bim > "$suffix".bim.ext
cat "$suffix".bim.ext | awk 'BEGIN{FS="\t"}{if($5 != $8 && $5 != $9){
  if($5 == "0"){
    if($6 == $9)print $2,$5,$6,$8,$9;
    else if($6 == $8)print $2,$5,$6,$9,$8;
    else if($6 == "0")print $2,"5",$6,$8,$9;
    else if($6 == $11) print $2,$5,$6,$8,$9;
    else if($6 == $10) print $2,$5,$6,$9,$8;
  }else{
    if($5 == $10)print $2,$5,$6,$8,$9;
    else if($5 == $11)print $2,$5,$6,$9,$8;
  }
}
else if($6 == "0"){
  if($5 == $8)print $2,$5,$6,$8,$9;
  else if($5 == $9)print $2,$5,$6,$9,$8;
}
}' > SNP_list_update.txt

# update the bim/bed files (complete the missing allele and fix strand issues)
plink --bfile "$suffix" --horse --update-alleles SNP_list_update.txt --make-bed --out "$suffix".update


# Fix the REF/ALT alleles if needed
plink --bfile "$suffix".update --horse --make-bed --a1-allele $arr70k 3 1 --out "$suffix".update.swab
#plink --bfile "$suffix" --horse --make-bed --a2-allele $arr70k 2 1 --out "$suffix".2


# convert you plink files to VCF
#plink2 --bfile "$suffix" --horse --output-chr chrM --not-chr y --ref-from-fa --fa $equCab2_ref --recode vcf --out "$suffix"
plink2 --bfile "$suffix".update.swab --horse --not-chr y --ref-from-fa --fa $equCab2_ref --recode vcf --out "$suffix"


# check for reference matching ## exclude unmatched records
#bcftools norm -c ws -f $equCab2_ref $suffix.vcf 1> $suffix.bcftools.check.vcf 2> $suffix.bcftools.log
bcftools norm -c ws -f $equCab2_ref $suffix.vcf 1> $suffix.bcftools.check.vcf 2> $suffix.bcftools.log ## 65155 out of 65157 variants loaded from 70K_Equcab2_Haflinger_Bellone.bim.

grep -v "^#" $suffix.vcf | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5}' > input.simple
grep -v "^#" $suffix.bcftools.check.vcf | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5}' > output.simple
paste input.simple output.simple > paste.simple
cat paste.simple | awk -F"\t" '{if($3 != $8)print}' | wc -l              ## 0
cat paste.simple | awk -F"\t" '{if($4 == $9)print}' | wc -l              ## 57782 Ref Allele did not change
cat paste.simple | awk -F"\t" '{if($4 == $9 && $5 == $10)print}' | wc -l ## 57782 &Alt allele did not change
cat paste.simple | awk -F"\t" '{if($4 == $9 && $5 == $10 && $10 != ".")print}' | wc -l            ## 50068    Alt != .
cat paste.simple | awk -F"\t" '{if($4 == $9 && $5 == $10 && $10 == ".")print}' | wc -l            ##  7714    Alt  = .  (i.e. sample homo for REF)


cat paste.simple | awk -F"\t" '{if($5 == $10)print}' | wc -l             ## 65155  Alt allele did not change
cat paste.simple | awk -F"\t" '{if($5 == $10 && $10 != ".")print}' | wc -l             ## 50078  Alt != .
cat paste.simple | awk -F"\t" '{if($5 == $10 && $10 == ".")print}' | wc -l             ## 15077  Alt  = .

cat paste.simple | awk -F"\t" '{if($4 == $10)print}' | wc -l             ## 0  Ref became Alt Alleles
cat paste.simple | awk -F"\t" '{if($4 != $10 && $5 == $9)print}' | wc -l ## 0

cat paste.simple | awk -F"\t" '{if($4 != $9)print}' | wc -l              ## 7373  Ref changed
cat paste.simple | awk -F"\t" '{if($4 != $9 && $4 != $10)print}' | wc -l ## 7373  &did not become Alt
cat paste.simple | awk -F"\t" '{if($4 != $9 && $4 != $10 && $5 == ".")print}' | wc -l ## 7363  xxx
cat paste.simple | awk -F"\t" '{if($4 != $9 && $4 != $10 && $5 != ".")print}' | wc -l ## 10    xxx


cat 70K_Equcab2_Haflinger_Bellone.bim_temp | awk 'BEGIN{FS=OFS="\t"}{if($1!="33")print $5,$6}' > bim.simple
paste paste.simple bim.simple > paste_bim.simple
cat paste_bim.simple | awk -F"\t" '{if($4 == $11 && $5 == $12)print}' | wc -l ## 24710  ## variants swiched by Plink2
cat paste_bim.simple | awk -F"\t" '{if($4 == $9 && $5 == $10 && $4 == $11 && $5 == $12)print}' | wc -l ## 24710 ## BCFtools has no effect on them
cat paste_bim.simple | awk -F"\t" '{if(($4 != $9 || $5 != $10) && $4 == $11 && $5 == $12)print}' | wc -l ## 0

cat paste_bim.simple | awk -F"\t" '{if($4 == $12 && ($5 == $11 || ($5 == "." && $11 == "0")))print}' | wc -l ## 39880 variants kept as it is by Plink2
cat paste_bim.simple | awk -F"\t" '{if($4 == $9 && $4 == $12 && ($5 == $11 || ($5 == "." && $11 == "0")))print}' | wc -l ## 33072 ## match Ref (difference=6808)
cat paste_bim.simple | awk -F"\t" '{if($4 != $9 && $4 == $12 && $5 == "." && $11 == "0")print}' | wc -l ## 6798 # one allele that does not match the reference
cat paste_bim.simple | awk -F"\t" '{if($4 != $9 && $4 == $12 && $5 == $11)print}' | wc -l ## 10 #strand error


cat paste_bim.simple | awk -F"\t" '{if($4 != $11 && $4 != $12)print}' | wc -l ## 565 ## ungenotyped variants
cat paste_bim.simple | awk -F"\t" '{if($4=="N" && $4 != $11 && $4 != $12)print}' | wc -l ## 565 ## variants changed by Plink2


awk -F "\t" 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next}{if(a[$3])print $0,a[$3];}' $work_dir/SNParrays/70k/GGP_Equine.simple.csv paste.simple > paste.simple.manifest


## conclusion:
VCF generated by Plink2 choose the proper Ref allele if there are 2 alleles & the strand is correct (+ve).
However,
a) if the strand is -ve, wrong alleles will be used
b) if there is only one allele,it will be used as Ref allele and the Alt allele will be ".". This will be error if the allele does match the ref (either because it is ALT or it REF but on -ve starnd)
c) if it ungenotyped, the Ref allele will be N and ALT will be "." ==> this is the only problem that BCFtools solve correctly
--ref-from-fa: 24710 variants changed, 33072 validated.  (total of 57782 out of 65155. i.e. unvalidated=7373)
 

Fixing this by BCFtools
REF/ALT total/modified/added:   65155/0/7373

##################
## explore
suffix="70K_Equcab2_PeurtoRicanPasoFino_Bellone"
cut -d" " -f1-20 $suffix.ped | head     ## sample_IDs
head -n5 $suffix.map                    ## SNP_IDs

grep "^#CH" $suffix.vcf | cut -f1-15          ## grep headers (compsite sample IDs)
grep -w "UKUL2" $suffix.vcf | cut -f10-20     ## grep by SNP_ID

head -n5 $suffix.hh  ## bad call: Family_ID Sample_ID SNP_ID

# print target sample genotypes in ped
snp_line=$(grep -n -w "BIEC2_1107310" $suffix.map | cut -d":" -f1)
snp_column=$(($((snp_line * 2)) + 5))
snp_column2=$((snp_column + 1))
grep -w "09-ES-20" $suffix.ped |cut -d" " -f"$snp_column"-"$snp_column2"

# print target sample genotypes in vcf
sample_column=$(grep "^#CH" $suffix.vcf | awk -F"\t" '{for (i = 10; i <= NF; i++){if($i == "11_09-ES-20")print i}}')
grep -w "BIEC2_1107310" $suffix.vcf | cut -f "$sample_column"

# from ./. in vcf to genotype in ped
grep -w "BIEC2_1107310" $suffix.vcf | cut -f22
grep "^#CH" $suffix.vcf | awk -F"\t" '{print $22}' ## 09-ES-97
grep -w "09-ES-97" $suffix.ped |cut -d" " -f"$snp_column"-"$snp_column2"

# print 1st 10 lines of target sample genotypes in vcf
grep "^#CH" $suffix.vcf | awk -F"\t" '{for (i = 10; i <= NF; i++){if($i == "3_09-ES-100")print i}}' ## 11
grep -v "^##" $suffix.vcf | head | awk -F"\t" '{print $11}'

#################################
#temp
dist="allStudies.nonAmb.dedup.cand1.pruned.dist.mibs"
cat $dist.id | tr '\t' '|' > $dist.id2
paste $dist.id $dist | grep -v "UN_" > $dist.comp

while read study sample;do
 id="$study|$sample"
 found=$(grep -wn "$id" $dist.id2)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 newBreed=$(cut -f 1,2,$((index+2)) $dist.comp | grep -v "$study$'\t'$sample" | sort -k3,3nr | head -n10 | cut -f2 | sed 's/_.*//' | sort | uniq -c | sort -k1,1nr | head -n1 | sed 's/ *//')
 echo $study $sample "$newBreed" | tr ' ' '\t'
 fi
done < <(cat *.newIds | cut -f3,4) > br_bestMatch

cat br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);gsub("\\..*","",$2);gsub("\\..*","",$4);if($2!=$4)print $2,$4}' | sort | uniq -c > br_bestMatch_summary
cat br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{q=$2;gsub("_.*","",q);gsub("\\..*","",q);m=$4;gsub("\\..*","",m);if(q!=m && q=="AP")print $0}' > br_bestMatch_AP




while read study sample;do
 id="$study.$sample"
 dist="allStudies.nonAmb.dedup.cand1.pruned.dist.mibs"
 found=$(grep -wn "$id" $dist.id)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 #newBreed=$(paste $dist.id  <(cut -f $index $dist) | grep -wv $study | grep -v "UN_" | sort -k3,3nr | head -n1 | cut -f2,3 | sed 's/_.*\t/\t/')
 newBreed=$(paste $dist.id  <(cut -f $index $dist) | grep -wv $study | grep -v "UN_" | sort -k3,3nr | head -n5 | cut -f2 | sed 's/_.*//' | sort | uniq -c | sort -k1,1nr | head -n1 | sed 's/ *//')
 echo $study $sample "$newBreed" | tr ' ' '\t'
 fi
done < <(cat 58TBD_52KRD_24PAR_272724_Unpruned.newIds | cut -f3,4) > 58TBD_52KRD_24PAR_272724_Unpruned.UNbreeds ## 54TB animals. The remaining 52KD,24AR.Per, and 4TB are duplicates in 670k_Pet2

while read study sample;do
 id="$study.$sample"
 dist="allStudies.nonAmb.dedup.cand1.pruned.dist.mibs"
 found=$(grep -wn "$id" $dist.id)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 newBreed=$(paste $dist.id  <(cut -f $index $dist) | grep -wv $study | grep -v "UN_" | sort -k3,3nr | head -n1 | cut -f2,3 | sed 's/_.*\t/\t/')
 echo $sample "$newBreed" | tr ' ' '\t'
 fi
done < <(cat MNEc2M_EquCab2.newIds | cut -f3,4 | grep "UN_") > MNEc2M_EquCab2.UNbreeds ## 8FM, 6HF, 6NOR, 1SN

while read study sample;do
 id="$study.$sample"
 dist="allStudies.nonAmb.dedup.cand1.pruned.dist.mibs"
 found=$(grep -wn "$id" $dist.id)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 newBreed=$(paste $dist.id  <(cut -f $index $dist) | grep -wv $study | grep -v "UN_" | sort -k3,3nr | head -n1 | cut -f2,3 | sed 's/_.*\t/\t/')
 echo $sample "$newBreed" | tr ' ' '\t'
 fi
done < <(cat GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.newIds | cut -f3,4) > GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.UNbreeds ## 5AR, 344FM, 7FR, 79FT, 1HF, 8HN, 72HO, 44ISH, 8KNB, 1LU, 116MRM, 19OL, 14OS, 4STB, 676SWB, 4TB, 7TB.Jpn, 3TK, 28WF

echo 670k_Bro1 AR_379
while read study sample;do
 id="$study.$sample"
 dist="allStudies.nonAmb.dedup.cand1.pruned.dist.mibs"
 found=$(grep -wn "$id" $dist.id)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 newBreed=$(paste $dist.id  <(cut -f $index $dist) | grep -wv $study | grep -v "UN_" | sort -k3,3nr | head -n1 | cut -f2,3 | sed 's/_.*\t/\t/')
 echo $sample "$newBreed" | tr ' ' '\t'
 fi
done


# tranform the distance matrix into one colum
cat allStudies.nonAmb.dedup.cand1.pruned.dist.mibs | awk 'BEGIN{FS="\t";OFS="\n";}{for (i = 1; i <= NF; i++)print $i}' > temp.pruned.1
# Prep the ids
awk 'BEGIN{OFS="\t"}FNR==NR{a[FNR]=$1;next}{for (i in a)print $1,a[i]}'  <(cat allStudies.nonAmb.dedup.cand1.pruned.dist.mibs.id | tr '\t' '.') <(cat allStudies.nonAmb.dedup.cand1.pruned.dist.mibs.id | tr '\t' '.') > temp.pruned.2
paste temp.pruned.2 temp.pruned.1 > ibs.pruned.all
cat ibs.pruned.all | awk -F"\t" '{if($3>0.7)print}' > ibs.pruned.all.reduced
cat ibs.pruned.all.reduced | tr '_' '\t' | awk 'BEGIN{FS=OFS="\t"}{if($1"_"$2!=$4"_"$5)print $1"_"$2"_"$3,$4"_"$5"_"$6,$7}' > ibs.pruned.all.reduced2
cat ibs.pruned.all.reduced2 | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2"."$1]){a[$1"."$2]=1;print}}' > ibs.pruned.all.reduced3

cat ibs.pruned.all.reduced3 | tr '_' '\t' | sed 's/\./\t/' | awk 'BEGIN{FS=OFS="\t"}{print $5,$6,$7,$1,$2,$3,$4,$8}'  | sed 's/\./\t/' | awk 'BEGIN{FS=OFS="\t"}{if($3!=$7){if(!a[$3"."$7])a[$7"."$3]=1;if(a[$7"."$3])print $7,$3;else print $3,$7;}}' | sort | uniq -c > breedSim


     99 AR.Per  AR
      8 AR.Per  KD
     27 AR      TB
      1 AR      TB.Jpn
      2 AR      FT

   2525 BE      CL

     83 ICH     EX

      6 KNB     AP
      5 KNB     SWB

      3 QH      AP
      2 QH      MIY
      2 QH      TB

      3 SWB     FT

    426 TB      TB.Jpn
      4 TB      STB

     82 WB      SWB
     22 WB      FT
     13 WB      FR
      2 WB      HN
      9 WB      HO
      2 WB      ISH
      4 WB      OL
      3 WB      OS
      3 WB      WF



#######################################
##temp
## Assess similarities
## prep data set from selected SNPs
tail -n+2 sel.list | awk -F"\t" '{print $2}' > sel.list.in
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --extract sel.list.in --make-bed --out allStudies.nonAmb.dedup.cand1.sel1000

## remove extra-chr & make indvidual id unique
plink --bfile allStudies.nonAmb.dedup.cand1.sel1000 --chr-set 31 no-xy --allow-extra-chr 0 --allow-no-sex --recode --out allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft
cp allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft.ped allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft.ped_temp
cat allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft.ped_temp | awk '{br=$2;gsub("_.*","",br);$2=$1"."$2;$1=br;$6=1;print $0}' > allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft.ped
rm allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft.ped_temp

## pca again
#plink --file allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft  --pca header tabs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --out allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft.pca ## keep failing

## try eigensoft
conda install -c bioconda eigensoft
cp $scripts/spca_parfile_sel1000 .
smartpca -p spca_parfile_sel1000 > smartpca_sel1000.log

##visualize
cp $scripts/plot_sel1000.R .
cat allStudies.nonAmb.dedup.cand1.sel1000.forEigensoft.ped | \
    awk '{if(!a[$1]){n+=1;a[$1]=n;} \
    if(a[$1]<=18)x=a[$1];else x=(a[$1]-1)%18+1; \
    c=int((a[$1]-1)/18);if(c==0)col="black";else if(c==1)col="red";else if(c==2)col="green";else col="blue";
    print $1,$2,x,col}' > st_sample.map
cat st_sample.map | awk '{print $1,$3,$4}' | sort | uniq > st.dat
cat st_sample.map | cut -d" " -f1 | sort | uniq -c | sort -k1,1nr > br.freq
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/br.freq remote_UCDavis_GoogleDr:Horse_parentage_share

Rscript plot_sel1000.R "PC1" "PC2"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC2.out_sel1000.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share
Rscript plot_sel1000.R "PC1" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC3.out_sel1000.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share
Rscript plot_sel1000.R "PC1" "PC4"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC4.out_sel1000.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share
Rscript plot_sel1000.R "PC1" "PC5"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC5.out_sel1000.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share
Rscript plot_sel1000.R "PC2" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC2.PC3.out_sel1000.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share




#########################
## temp
while read study sample;do
 id="$study.$sample" #"50k_McCue1.BE_31"
 dist="allStudies.nonAmb.dedup.cand1.pruned.dist.mibs"
 found=$(grep -wn "$id" $dist.id)
 if [ ! "$found" ];then
 grep $sample ibs.gen.95
 fi
done < <(cat 58TBD_52KRD_24PAR_272724_Unpruned.newIds | cut -f3,4) > 58TBD_52KRD_24PAR_272724_Unpruned.UNbreeds2


50k_McCue1.BE_31        80k_Pet2.CL_1   0.746893
50k_McCue1.BE_31        80k_Pet2.CL_2   0.745885
50k_McCue1.BE_31        80k_Pet2.CL_3   0.740412
50k_McCue1.BE_31        80k_Pet2.CL_4   0.735814
50k_McCue1.BE_31        80k_Pet2.CL_5   0.721375
50k_McCue1.BE_31        80k_Pet2.CL_6   0.727483
50k_McCue1.BE_31        80k_Pet2.CL_7   0.719545
50k_McCue1.BE_31        80k_Pet2.CL_8   0.728487
50k_McCue1.BE_31        80k_Pet2.CL_9   0.751873
50k_McCue1.BE_31        80k_Pet2.CL_10  0.733906
50k_McCue1.BE_31        80k_Pet2.CL_11  0.737981
50k_McCue1.BE_31        80k_Pet2.CL_12  0.724351
50k_McCue1.BE_31        80k_Pet2.CL_13  0.739711
50k_McCue1.BE_31        80k_Pet2.CL_14  0.735269
50k_McCue1.BE_31        80k_Pet2.CL_15  0.73013


670k_Faw1.TB.Jpn_17     MNEc2M_McCue2.TB_257    0.730423

50k_McCue9.TB_135       670k_Bro1.AR_384        0.706089
50k_McCue9.TB_125       670k_Bro1.AR_391        0.711918


80k_Bel6.AP_162 80k_Bel6.KNB_77 0.753411
80k_Bel6.AP_162 80k_Bel6.KNB_118        0.700881
80k_Bel6.AP_162 80k_Bel6.KNB_163        0.76024
80k_Bel6.AP_162 80k_Bel6.KNB_167        0.769765
80k_Bel6.AP_162 80k_Bel6.KNB_169        0.764093
80k_Bel6.AP_172 80k_Bel6.KNB_163        0.702025
80k_Bel6.KNB_93 670k_Bel3.SP_167        0.700266
80k_Bel6.KNB_128        670k_Mik1.SWB_178       0.717087
80k_Bel6.KNB_128        670k_Mik1.SWB_181       0.712799
80k_Bel6.KNB_128        670k_Mik1.SWB_305       0.703026
80k_Bel6.KNB_128        670k_Mik1.SWB_310       0.707497

80k_Bel7.ICH_2  670k_LV1.EX_92  0.703982
80k_Bel7.ICH_2  670k_LV1.EX_180 0.710291
80k_Bel7.ICH_3  670k_LV1.EX_180 0.707648
80k_Bel7.ICH_6  670k_LV1.EX_181 0.703914
80k_Bel7.ICH_10 670k_LV1.EX_134 0.701518
80k_Bel7.ICH_10 670k_LV1.EX_180 0.708313
80k_Bel7.ICH_17 670k_LV1.EX_92  0.702893
80k_Bel7.ICH_17 670k_LV1.EX_97  0.701397
80k_Bel7.ICH_20 670k_LV1.EX_92  0.708976
670k_LV1.EX_92  670k_LV2.ICH_2  0.703759
670k_LV1.EX_92  670k_LV2.ICH_3  0.724613
670k_LV1.EX_92  670k_LV2.ICH_9  0.701833
670k_LV1.EX_92  670k_LV2.ICH_11 0.708856
670k_LV1.EX_92  670k_LV2.ICH_13 0.700359
670k_LV1.EX_92  670k_LV2.ICH_14 0.705748
670k_LV1.EX_92  670k_LV2.ICH_15 0.711111
670k_LV1.EX_92  670k_LV2.ICH_18 0.710558
670k_LV1.EX_92  670k_LV2.ICH_20 0.708699
670k_LV1.EX_92  670k_LV2.ICH_22 0.706181
670k_LV1.EX_92  670k_LV2.ICH_27 0.707369
670k_LV1.EX_92  670k_LV2.ICH_28 0.701873
670k_LV1.EX_92  670k_LV2.ICH_37 0.703734
670k_LV1.EX_92  670k_LV2.ICH_40 0.716009

80k_KWPN1.WB_14 MNEc2M_McCue2.HN_58     0.768331

80k_KWPN1.WB_3  670k_Mik1.SWB_69        0.713182
80k_KWPN1.WB_16 670k_Mik1.SWB_106       0.712522
80k_KWPN1.WB_25 670k_Mik1.SWB_69        0.703135
80k_KWPN1.WB_47 670k_Mik1.SWB_279       0.70379
80k_KWPN1.WB_48 670k_Mik1.SWB_31        0.709263
80k_KWPN1.WB_110        670k_Mik1.SWB_15        0.715808
80k_KWPN1.WB_110        670k_Mik1.SWB_51        0.708853
80k_KWPN1.WB_110        670k_Mik1.SWB_180       0.716722
80k_KWPN1.WB_110        670k_Mik1.SWB_371       0.714386
80k_KWPN1.WB_119        670k_Mik1.SWB_8 0.700746
80k_KWPN1.WB_119        670k_Mik1.SWB_131       0.717112
80k_KWPN1.WB_128        670k_Mik1.SWB_184       0.722059
80k_KWPN1.WB_139        670k_Mik1.SWB_184       0.702088
80k_KWPN1.WB_159        670k_Mik1.SWB_288       0.701447
80k_KWPN1.WB_208        670k_Mik1.SWB_26        0.700869
80k_KWPN1.WB_259        670k_Mik1.SWB_11        0.716175

80k_KWPN1.WB_134        670k_Bel3.FR_7  0.700554
80k_KWPN1.WB_134        670k_Bel3.FR_15 0.722026
80k_KWPN1.WB_134        670k_Bel3.FR_18 0.702643
80k_KWPN1.WB_134        670k_Bel3.FR_19 0.705674
80k_KWPN1.WB_134        670k_Bel3.FR_135        0.704366
80k_KWPN1.WB_134        670k_Bel3.FR_143        0.715208
80k_KWPN1.WB_134        670k_Bel3.FR_144        0.707372
80k_KWPN1.WB_134        670k_Bel3.FR_145        0.704992
80k_KWPN1.WB_134        670k_Bel3.FR_148        0.70215
80k_KWPN1.WB_134        670k_Bel3.FR_151        0.702216
80k_KWPN1.WB_134        670k_Bel3.FR_153        0.70502
80k_KWPN1.WB_134        670k_Bel3.FR_154        0.701942

80k_KWPN1.WB_151        670k_Hill3.ISH_12       0.702335

80k_KWPN1.WB_159        MNEc2M_McCue2.FT_50     0.7539
80k_KWPN1.WB_269        MNEc2M_McCue2.FT_50     0.709823
80k_KWPN1.WB_395        MNEc2M_McCue2.FT_50     0.701614
80k_KWPN1.WB_430        MNEc2M_McCue2.FT_50     0.716877
80k_KWPN1.WB_491        MNEc2M_McCue2.FT_50     0.706944

80k_KWPN1.WB_248        80k_Vit1.OL_13  0.727776
80k_KWPN1.WB_333        80k_Vit1.OL_13  0.710572

80k_KWPN1.WB_309        80k_Vit1.OS_6   0.70546
80k_KWPN1.WB_660        80k_Vit1.OS_6   0.703473

670k_Ant1.AR.Per_1      670k_Bro1.AR_266        0.727385
670k_Ant1.AR.Per_1      670k_Bro1.AR_267        0.721267
670k_Ant1.AR.Per_1      670k_Bro1.AR_271        0.700527
670k_Ant1.AR.Per_1      670k_Bro1.AR_272        0.707513

##########################
## temp

## New trial
## exclude variants with missing calling rate > 0.05 (i.g. genotyping rate > 0.95 to keep SNPs across arrays) or minor allele frequency < 0.3
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --geno 0.05 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.cand1
#605271 variants removed due to missing genotype data (--geno).
#3692 variants removed due to Hardy-Weinberg exact test.
#15124 variants removed due to minor allele threshold(s)
#10455 variants and 8471 samples pass filters and QC.

# calc the freq and generate a histogram
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq -out cand1  ## cand1.frq
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 cand1.frq | awk '{print $5}') > cand1.frq.hist

# calc the count and merge with the freq output
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq counts -out cand1 ## cand1.frq.counts
paste cand1.frq cand1.frq.counts | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$11,$6}' > cand1.frq.merged

# add hwe info
#plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freqx -out cand1 ## cand1.frqx
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --hardy -out cand1 ## cand1.hwe
paste cand1.frq.merged cand1.hwe | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$13,$14,$15,$16}' > cand1.frq.merged2
#tail -n+2 cand1.frq.merged2 | awk '{if(($11+0)<1e-50)print}' | wc -l ## filter like --hwe

# add info per study
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq --family -out cand1 ## cand1.frq.strat
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1 OFS $2]=$3 OFS $4 OFS $5 OFS $6 OFS $7;next}{if(a[$1 OFS $2]){print $1,$2,"allStudies",a[$1 OFS $2];a[$1 OFS $2]=0} print $1,$2,$3,$4,$5,$6,$7,$8}' <(tail -n+2 cand1.frq.merged) cand1.frq.strat > cand1.frq.strat.merged
#tail -n+2 cand1.frq.strat | awk 'BEGIN{OFS="\t"}{inp[$1 FS $2]=1;if($8)a[$1 FS $2]+=1;if($6>0.1)b[$1 FS $2]+=1;if($6>0.2)c[$1 FS $2]+=1;if($6>0.3)d[$1 FS $2]+=1;if($6>0.4)e[$1 FS $2]+=1;}END{for (i in inp)print i,a[i],b[i],c[i],d[i],e[i];}' > cand1.frq.strat.summary
echo $(head -n1 cand1.frq.merged2) "st" "st_0.1" "st_0.2" "st_0.3" "st_0.4" | tr ' ' '\t' > cand1.frq.merged3
awk 'BEGIN{OFS="\t"}FNR==NR{x=$1 FS $2;if($8)a[x]+=1;if($6>0.1)b[x]+=1;if($6>0.2)c[x]+=1;if($6>0.3)d[x]+=1;if($6>0.4)e[x]+=1;next}{i=$1 FS $2;print $0,(a[i]>0?a[i]:0),(b[i]>0?b[i]:0),(c[i]>0?c[i]:0),(d[i]>0?d[i]:0),(e[i]>0?e[i]:0);}' <(tail -n+2 cand1.frq.strat) <(tail -n+2 cand1.frq.merged2) >> cand1.frq.merged3

# add legacy info
cat cand1.frq.strat | awk '{if($3~"^50k_" && $8>0)print $2}' | sort | uniq > 50k.list
cat cand1.frq.strat | awk '{if($3~"^70k_" && $8>0)print $2}' | sort | uniq > 70k.list
echo "SNP legacy" | tr ' ' '\t' > legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]="50k";next}{if(a[$1])a[$1]="50_70k";else a[$1]="70k"}END{for (i in a)print i,a[i]}' 50k.list 70k.list >> legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"none";}' legacy.list cand1.frq.merged3 > cand1.frq.merged4

# Pruning
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --indep 10 2 2 --out allStudies.nonAmb.dedup.cand1.toPrune
#Pruning complete.  3062 of 10455 variants removed.
#plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --extract allStudies.nonAmb.dedup.cand1.toPrune.prune.in --make-bed --out allStudies.nonAmb.dedup.cand1.pruned
echo $(head -n1 cand1.frq.merged4) "select" | tr ' ' '\t' > cand1.frq.merged5
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$2])print $0,"exclude";else print $0,"keep"}' allStudies.nonAmb.dedup.cand1.toPrune.prune.out <(tail -n+2 cand1.frq.merged4) >> cand1.frq.merged5

# select candidates
#cat cand1.frq.merged4 | awk '{if(($10-$9)<-0.05)print}' > excessHet.list ## 46(<-0.1) 143(<-0.05) 2023(<0) ## excess heterozygosity has almost no effect
#cat cand1.frq.merged4 | awk '{if($12==54 && $15>25 && ($11+0)>1e-50 && (($10-$9)>0) && $5>0.4)print}' > sel.list ## 4215
#echo $(head -n1 cand1.frq.merged4) "select" | tr ' ' '\t' > cand1.frq.merged5
#tail -n+2 cand1.frq.merged4 | awk 'BEGIN{OFS="\t"}{if($12==54 && $15>25 && ($11+0)>1e-50 && (($10-$9)>0) && $5>0.4)print $0,1;else print $0,0}' >> cand1.frq.merged5

# Rank the candidates
tail -n+2 cand1.frq.merged5 | awk '{if($18=="keep")print $0}' | sort -k15,15nr | awk '{print $2,NR}' > cand1.rank
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])$18=a[$2];print $0;}' cand1.rank cand1.frq.merged5 > cand1.frq.merged6

# sel 1000 targets
head -n1 cand1.frq.merged6 > sel.list
tail -n+2 cand1.frq.merged6 | awk -F"\t" '{if($18<=1000)print}' | sort -k18,18n >> sel.list

#########################
## first trial
## exclude variants with missing calling rate > 0.7 (i.g. genotyping rate > 0.3) or minor allele frequency < 0.05
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --geno 0.7 --maf 0.05 --make-bed --out allStudies.nonAmb.dedup.cand1
#105126 variants removed due to missing genotype data (--geno).
#123315 variants removed due to minor allele threshold(s)
#406101 variants and 8471 samples pass filters and QC.

# calc the freq and generate a histogram
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq -out cand1  ## cand1.frq
awk -v size=0.05 '{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 cand1.frq | awk '{print $5}') > cand1.frq.hist

# calc the count and merge with the freq output
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq counts -out cand1 ## cand1.frq.counts
paste cand1.frq cand1.frq.counts | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$11,$6}' > cand1.frq.merged

# add hwe info
#plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freqx -out cand1 ## cand1.frqx
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --hardy -out cand1 ## cand1.hwe
paste cand1.frq.merged cand1.hwe | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$13,$14,$15,$16}' > cand1.frq.merged2
#tail -n+2 cand1.frq.merged2 | awk '{if(($11+0)<1e-10)print}' | wc -l ## filter like --hwe

# add info per study
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq --family -out cand1 ## cand1.frq.strat
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1 OFS $2]=$3 OFS $4 OFS $5 OFS $6 OFS $7;next}{if(a[$1 OFS $2]){print $1,$2,"allStudies",a[$1 OFS $2];a[$1 OFS $2]=0} print $1,$2,$3,$4,$5,$6,$7,$8}' <(tail -n+2 cand1.frq.merged) cand1.frq.strat > cand1.frq.strat.merged
#tail -n+2 cand1.frq.strat | awk 'BEGIN{OFS="\t"}{inp[$1 FS $2]=1;if($8)a[$1 FS $2]+=1;if($6>0.1)b[$1 FS $2]+=1;if($6>0.2)c[$1 FS $2]+=1;if($6>0.3)d[$1 FS $2]+=1;if($6>0.4)e[$1 FS $2]+=1;}END{for (i in inp)print i,a[i],b[i],c[i],d[i],e[i];}' > cand1.frq.strat.summary
echo $(head -n1 cand1.frq.merged2) "st" "st_0.1" "st_0.2" "st_0.3" "st_0.4" | tr ' ' '\t' > cand1.frq.merged3
awk 'BEGIN{OFS="\t"}FNR==NR{x=$1 FS $2;if($8)a[x]+=1;if($6>0.1)b[x]+=1;if($6>0.2)c[x]+=1;if($6>0.3)d[x]+=1;if($6>0.4)e[x]+=1;next}{i=$1 FS $2;print $0,(a[i]>0?a[i]:0),(b[i]>0?b[i]:0),(c[i]>0?c[i]:0),(d[i]>0?d[i]:0),(e[i]>0?e[i]:0);}' <(tail -n+2 cand1.frq.strat) <(tail -n+2 cand1.frq.merged2) >> cand1.frq.merged3

# add legacy info
cat cand1.frq.strat | awk '{if($3~"^50k_" && $8>0)print $2}' | sort | uniq > 50k.list
cat cand1.frq.strat | awk '{if($3~"^70k_" && $8>0)print $2}' | sort | uniq > 70k.list
echo "SNP legacy" | tr ' ' '\t' > legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]="50k";next}{if(a[$1])a[$1]="50_70k";else a[$1]="70k"}END{for (i in a)print i,a[i]}' 50k.list 70k.list >> legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"none";}' legacy.list cand1.frq.merged3 > cand1.frq.merged4

# select candidates
#cat cand1.frq.merged4 | awk '{if(($10-$9)<-0.05)print}' > excessHet.list ## 46(<-0.1) 143(<-0.05) 2023(<0) ## excess heterozygosity has almost no effect
cat cand1.frq.merged4 | awk '{if($12==54 && $15>25 && ($11+0)>1e-50 && (($10-$9)>0) && $5>0.4)print}' > sel.list ## 4215
echo $(head -n1 cand1.frq.merged4) "select" | tr ' ' '\t' > cand1.frq.merged5
tail -n+2 cand1.frq.merged4 | awk 'BEGIN{OFS="\t"}{if($12==54 && $15>25 && ($11+0)>1e-50 && (($10-$9)>0) && $5>0.4)print $0,1;else print $0,0}' >> cand1.frq.merged5

# compare to previous lists
rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs/media-1.txt $HOME/Horse_parentage_SNPs/backup_original
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8 FS $9]=$1;next}{if(a[$5 FS $6])$17=a[$5 FS $6];else if($3)$17=$3;else $17="N/A"; print $0;}' $map670 <(cat ../backup_original/media-1.txt | sed -e "s/\r//g" | sed '/^[[:space:]]*$/d') | cut -f1-17 | sed 's/SNP Name - Affymetrix 670K$/SNP/'> media-1-TM.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$0;next}{if(a[$17])print $0,a[$17];else print $0,"not_in_1ry_list";}' cand1.frq.merged5 media-1-TM.txt > media-1-TM_ext.txt
tail -n+2 media-1-TM_ext.txt | rev | cut -f 1 | rev | sort | uniq -c
cat media-1-TM_ext.txt | awk -F"\t" '{if($18=="not_in_1ry_list" && $17!="N/A")print $17}' | sed 's/_DUP.*//' | grep -Fwf - $map670
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34!="50_70k")print $35}' | uniq -c
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $22<=0.4)print $35}' | uniq -c
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $22>0.4 && $29<54)print $35}' | uniq -c
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $22>0.4 && $29==54 && $35!=1)print $0}'

## 288 are among my selected SNPs
## 59 marker failed to remap & 91 does not show up in both 50k and 70k arrays & 67 has MAF < 0.4 (all of them are > 0.3) & 62 failed to show up in all studies we have (most of them did not show up in a study or two) & 10 failed the HWE testing

## Task: upload cand1.frq.merged5
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/cand1.frq.merged5 remote_UCDavis_GoogleDr:Horse_parentage_share
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/media-1-TM_ext.txt remote_UCDavis_GoogleDr:Horse_parentage_share

###########################
## temp
## rate of heterozygosity
tail -n+2 cand1.frq.merged3 | awk 'BEGIN{FS=OFS="\t"}{split($8,a,"/");print a[2]/$7}' | head #awk '{if($10==54 && $13>25 && ($9+0)>0.1 && $5>0.4)print}' | head
tail -n+2 cand1.frq.merged3 | awk 'BEGIN{FS=OFS="\t"}{split($8,a,"/");print $0,a[2]/$7}' | awk '{if($15>0.25)print}' | head
###########################
## temp
## exclude variants with Hardy-Weinberg equilibrium exact test p-value below 1e-10.
#plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --hwe 1e-10  midp --make-bed --out allStudies.nonAmb.dedup.cand2
##187983 variants removed due to Hardy-Weinberg exact test.
##218118 variants and 8471 samples pass filters and QC.
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --hwe 1e-10 --make-bed --out allStudies.nonAmb.dedup.cand2
#186632 variants removed due to Hardy-Weinberg exact test.
#219469 variants and 8471 samples pass filters and QC.

# calc the freq and generate a histogram
plink --bfile allStudies.nonAmb.dedup.cand2 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq -out cand2  ## cand2.frq
awk -v size=0.05 '{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 cand2.frq | awk '{print $5}') > cand2.frq.hist
# calc the count and merge with the freq output
plink --bfile allStudies.nonAmb.dedup.cand2 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq counts -out cand2 ## cand2.frq.counts
paste cand2.frq cand2.frq.counts | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$11,$6}' > cand2.frq.merged
# calc the freq and count per study & merge with freq and count in all studies
plink --bfile allStudies.nonAmb.dedup.cand2 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq --family -out cand2 ## cand2.frq.strat
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1 OFS $2]=$3 OFS $4 OFS $5 OFS $6 OFS $7;next}{if(a[$1 OFS $2]){print $1,$2,"allStudies",a[$1 OFS $2];a[$1 OFS $2]=0} print $1,$2,$3,$4,$5,$6,$7,$8}' <(tail -n+2 cand2.frq.merged) cand2.frq.strat > cand2.frq.strat.merged

paste cand1.frq.hist cand2.frq.hist | awk '{if($6)print $3/$6}'



#########################################
## Temp
plink --bfile allStudies.nonAmb --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr 19 --rel-cutoff 0.95

64297 MB RAM detected; reserving 32148 MB for main workspace.
14871 out of 634542 variants loaded from .bim file.
8855 samples (550 males, 855 females, 7450 ambiguous) loaded from .fam.
Ambiguous sex IDs written to plink.nosex .
620 phenotype values loaded from .fam.
Using up to 31 threads (change this with --threads).
Before main variant filters, 8765 founders and 90 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate is 0.393133.
14871 variants and 8855 samples pass filters and QC (before --rel-cutoff).
Among remaining phenotypes, 142 are cases and 478 are controls.  (8235
phenotypes are missing.)
Relationship matrix calculation complete.
688 samples excluded by --rel-cutoff.
Remaining sample IDs written to plink.rel.id .



plink --bfile allStudies.nonAmb --distance square allele-ct ibs 1-ibs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr 19 --out allStudies.nonAmb.distance
#14871 variants and 8855 samples pass filters and QC.
cat allStudies.nonAmb.distance.mibs | awk 'BEGIN{FS="\t";OFS="\n";}{for (i = 1; i <= NF; i++)print $i}' > temp.1
awk 'BEGIN{OFS="\t"}FNR==NR{a[FNR]=$1;next}{for (i in a)print $1,a[i]}'  <(cat allStudies.nonAmb.distance.mdist.id | tr '\t' '.') <(cat allStudies.nonAmb.distance.mdist.id | tr '\t' '.') > temp.2
paste temp.2 temp.1 | awk 'BEGIN{FS=OFS="\t"}{if($3>0.95 && $1!=$2)print}' > ibs.95_temp
cat ibs.95_temp | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2"."$1] && $3!="nan"){a[$1"."$2]=1;print}}' > ibs.95


plink --bfile allStudies.nonAmb --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --geno 0.05 --indep 50 5 2 --out allStudies.nonAmb.toPrune
#605940 variants removed due to missing genotype data (--geno).
#395 variants removed due to minor allele threshold(s)
#Pruning complete.  14516 of 28207 variants removed.
plink --bfile allStudies.nonAmb --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --extract allStudies.nonAmb.toPrune.prune.in --make-bed --out allStudies.nonAmb.pruned
plink --bfile allStudies.nonAmb.pruned --distance square allele-ct ibs 1-ibs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --out allStudies.nonAmb.pruned.distance
cat allStudies.nonAmb.pruned.distance.mibs | awk 'BEGIN{FS="\t";OFS="\n";}{for (i = 1; i <= NF; i++)print $i}' > temp.pr.1
awk 'BEGIN{OFS="\t"}FNR==NR{a[FNR]=$1;next}{for (i in a)print $1,a[i]}'  <(cat allStudies.nonAmb.pruned.distance.mdist.id | tr '\t' '.') <(cat allStudies.nonAmb.pruned.distance.mdist.id | tr '\t' '.') > temp.pr.2
paste temp.pr.2 temp.pr.1 | awk 'BEGIN{FS=OFS="\t"}{if($3>0.95 && $1!=$2)print}' > ibs.pr.95_temp
cat ibs.pr.95_temp | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2"."$1] && $3!="nan"){a[$1"."$2]=1;print}}' > ibs.pr.95


plink --bfile allStudies.nonAmb --distance square allele-ct ibs 1-ibs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --geno 0.01 --out allStudies.nonAmb.geno.distance2
#15753 variants and 8855 samples pass filters and QC.
cat allStudies.nonAmb.geno.distance2.mibs | awk 'BEGIN{FS="\t";OFS="\n";}{for (i = 1; i <= NF; i++)print $i}' > temp.gen2.1
awk 'BEGIN{OFS="\t"}FNR==NR{a[FNR]=$1;next}{for (i in a)print $1,a[i]}'  <(cat allStudies.nonAmb.geno.distance2.mdist.id | tr '\t' '.') <(cat allStudies.nonAmb.geno.distance2.mdist.id | tr '\t' '.') > temp.gen2.2
paste temp.gen2.2 temp.gen2.1 | awk 'BEGIN{FS=OFS="\t"}{if($3>0.95 && $1!=$2)print}' > ibs.gen2.95_temp
cat ibs.gen2.95_temp | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2"."$1] && $3!="nan"){a[$1"."$2]=1;print}}' > ibs.gen2.95



## Temp
# convert you plink files to VCF
plink2 --bfile "$suffix".swab --chr-set 31 no-xy --allow-extra-chr --ref-from-fa --fa $equCab3_ref --recode vcf id-delim="." --out "$suffix"

# check for reference matching ## exclude unmatched records
bcftools norm -c ws -f $equCab3_ref $suffix.vcf 1> $suffix.bcftools.check.vcf 2> $suffix.bcftools.log ## 65155 out of 65157 variants loaded from 70K_Equcab2_Haflinger_Bellone.bim.

f1=670K_Equcab3_Belgian_Haflinger_Bellone.newIds.update.swab
f2=70K_Equcab2_PeurtoRicanPasoFino_Bellone.newIds.update.swab
plink --bfile $f1 --chr-set 31 no-xy --allow-extra-chr --bmerge $f2 --make-bed --out f1f2 --merge-mode 1


