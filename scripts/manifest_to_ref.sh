#!/bin/sh

#  manifest_to_ref.sh
#  
#
#  Created by Tamer Mansour on 3/31/22.
#  

## Input manifest is tab separated file with at least 8 columns
## col1="SNP ID", col3="chr", col4="pos", col7&8="2 alleles"
manifest="$1"
vcfContigs="$2"
ref="$3"
out="$4"


# create VCF template using the positions of the array SNPs
echo "##fileformat=VCFv4.3" > ${out}_pos.vcf
cat $vcfContigs >> ${out}_pos.vcf
echo "#CHROM POS ID REF ALT QUAL FILTER INFO" | tr ' ' '\t' >> ${out}_pos.vcf
tail -n+2 "$manifest" | awk 'BEGIN{FS=OFS="\t"}{if($3!="Y")print $3,$4,$1,"N",".",".",".","."}' >> ${out}_pos.vcf

# Use bcftools to obtain the reference alleles of the array SNPs
bcftools norm -c ws -f $ref ${out}_pos.vcf 1> ${out}_ref.vcf 2> ${out}_ref.vcf.log

# align the vcf carring reference alleles with the vcf carring positions (which already matches the manifest)
cp ${out}_ref.vcf ${out}_ref.vcf_temp
grep "^#" ${out}_ref.vcf_temp > ${out}_ref.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$3]=$4;next;}{$4=a[$3];print}' <(grep -v "^#" ${out}_ref.vcf_temp) <(grep -v "^#" ${out}_pos.vcf) >> ${out}_ref.vcf
rm ${out}_ref.vcf_temp

# Transform the array SNP alleles into reference and alternative alleles
grep -v "^#" ${out}_ref.vcf | awk 'BEGIN{FS=OFS="\t"}{print $4}' | tr 'tcga' 'TCGA' > ${out}.ref_alleles
tail -n+2 "$manifest" | awk 'BEGIN{FS=OFS="\t"}{if($3!="Y")print $7,$8}' > ${out}.chip_alleles
cat ${out}.chip_alleles | tr 'TCGA' 'AGCT' > ${out}.chip_alleles_oppStrand
paste "${out}.ref_alleles" "${out}.chip_alleles" "${out}.chip_alleles_oppStrand" | awk 'BEGIN{FS=OFS="\t";a=b=c=d=e=0}{if($1==$2)a+=1;else if($1==$3)b+=1;else if($1==$4)c+=1;else if($1==$5)d+=1;else e+=1;}END{print a,b,c,d,e;}'
paste "${out}.ref_alleles" "${out}.chip_alleles" "${out}.chip_alleles_oppStrand" | awk 'BEGIN{FS=OFS="\t"}{if($1==$2)print $2,$3;else if($1==$3)print $3,$2;else if($1==$4)print $4,$5;else if($1==$5)print $5,$4;else print $1,"."}' > ${out}.refalt_alleles
rm ${out}.ref_alleles ${out}.chip_alleles ${out}.chip_alleles_oppStrand

# Generate an output VCF
grep "^#" ${out}_ref.vcf > ${out}.vcf
paste <(grep -v "^#" ${out}_ref.vcf | cut -f1-3) ${out}.refalt_alleles <(grep -v "^#" ${out}_ref.vcf | cut -f6-8) >> ${out}.vcf

# Generate an output tab separated map for ALT/REF alleles on +ve and -ve strands
# Note that the columns in the output represent the SNP_id, Alt allele, Ref allele, compAlt, compRef
grep -v "^#" ${out}.vcf | awk 'BEGIN{FS=OFS="\t"}{print $3,$5,$4}' > ${out}.alt_ref.tab
paste ${out}.alt_ref.tab <(cut -f2-3 ${out}.alt_ref.tab | tr 'ACGT' 'TGCA') > ${out}.alt_ref_comp.tab
