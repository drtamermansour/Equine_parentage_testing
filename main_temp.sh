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

