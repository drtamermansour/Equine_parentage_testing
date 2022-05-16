#!/bin/sh

#  main.sh
#  
#
#  Created by Tamer Mansour on 3/2/22.
#  Project name (Horse_parentage_SNPs)
##### setup the stage
work_dir=$(pwd)
scripts="$work_dir"/scripts
## Create specific conda env
conda update conda ## conda info ==> conda version : 4.11.0
conda create -n equineSNP
conda activate equineSNP
conda install -c bioconda plink
conda install -c bioconda plink2
conda install -c bioconda bcftools

##### Download expermintal data
mkdir -p $work_dir/backup_original && cd $work_dir/backup_original
## 1. cp files from google drive to Farm
module load rclone #On Feb 2022: Module rclone/1.53.3 loaded
rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs $HOME/Horse_parentage_SNPs/backup_original
rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs_private $HOME/Horse_parentage_SNPs/backup_original

## File naming issues ## Fix and update GoogleDr
#1
mv FM670K_Metadata.txt 670K_Equcab3_FranchesMontagnes_Gmel.gender.txt
mv Metadata_FM50K.txt 50K_Equcab3_FranchesMontagnes_Gmel.gender.txt
mv FM_duplicates.txt 50K.670K.duplicates_Equcab3_FranchesMontagnes_Gmel.txt
#2
unzip GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.zip
rm GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.zip
#3
mv SNP_Map.csv GGP65K_SNP_Map.csv
#4
mv 670K_Equcab3__Mixedbreed_Bellone.vcf 670K_Equcab3_Mixedbreed_Bellone.vcf
mv 670K_Equcab3__Mixedbreed_Bellone.xlsx 670K_Equcab3_Mixedbreed_Bellone.xlsx
#5
mv 670K_Equcab3_Belgian_Haflinger_Bellone_breed_log.xlsx 670K_Equcab3_Belgian_Haflinger_Bellone.xlsx

## 2. Other public data
# Our Array data
# 1. 670K data: Arabian: https://www.sciencedirect.com/science/article/pii/S0890850820302760?via%3Dihub
wget -O 670k_Equcab2_Arabian_Patterson.zip https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/zs5h47xsxc-1.zip
unzip 670k_Equcab2_Arabian_Patterson.zip
for f in EMS_2020_Patterson_et_al.{ped,map,fam};do
 newF=$(echo $f | sed 's/EMS_2020_Patterson_et_al/670k_Equcab2_Arabian_Patterson/')
 mv $f $newF
done


# 2. 670K data: Thoroughbred (Japan): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6655603/
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ936/ERZ936110/ThoroughbredsJpn_670K_Fawcett2019.vcf.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/ERZ936/ERZ936110/ThoroughbredsJpn_670K_Fawcett2019.vcf.gz.tbi
mv ThoroughbredsJpn_670K_Fawcett2019.vcf.gz 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.vcf.gz
mv ThoroughbredsJpn_670K_Fawcett2019.vcf.gz.tbi 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.vcf.gz.tbi

# 3. 670K data: Arabian (Persia)    https://academic.oup.com/jhered/article/110/2/173/5264250
# download locally, unzip, upload (IR.paper.heredity.bed/bim) to the drive, then rclone

# 4. 670K data: Arabian    https://www.nature.com/articles/s41598-020-66232-1?fbclid=IwAR0ILYf_eIXrHcEhEpExRrgo7mnUfXwVo7Kiu8Xow2qlDS7FdmMthTs67Ms#data-availability
# download locally, upload (AHSdataset.bed/bim) to the drive, then rclone


# 5. 670K data: Kurdish/Arabian (Persia)/Thoroughbred (US): https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0247123
# download the largest unpruned dataset (58TBD_52KRD_24PAR_272724_Unpruned.*) locally, upload to the drive, then rclone


# Tosso Leeb <tosso.leeb@vetsuisse.unibe.ch>: WGS of 88 horses - PRJEB28306 - Publication: https://onlinelibrary.wiley.com/doi/10.1111/age.12753
# download links: http://ftp.ebi.ac.uk/pub/databases/eva/PRJEB28306/
# rename to WGS_EquCab3_MulitpleBreeds_TossoLeeb
wget http://ftp.ebi.ac.uk/pub/databases/eva/PRJEB28306/horses.88.vars.flt.pass.ebi.vcf.gz
wget http://ftp.ebi.ac.uk/pub/databases/eva/PRJEB28306/horses.88.vars.flt.pass.ebi.vcf.gz.tbi

# Molly and Sian's publication
wget http://ftp.ebi.ac.uk/pub/databases/eva/PRJEB47918/thesis_intersect_pub.decomposed.vcf.gz
wget http://ftp.ebi.ac.uk/pub/databases/eva/PRJEB47918/thesis_intersect_pub.decomposed.vcf.gz.tbi
#############
## Transform all files to map/ped
#1
#zcat 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.vcf.gz --double-id --allow-extra-chr --recode --out "670k_Equcab2_ThoroughbredsJpn_Fawcett2019"
#2
#cat 670K_EquCab2_Standardbreds_Bellone.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_EquCab2_Standardbreds_Bellone.vcf --double-id --chr-set 35 --recode --out "670K_EquCab2_Standardbreds_Bellone"

#3
#cat 670K_Equcab3_Mixedbreed_Bellone.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab3_Mixedbreed_Bellone.vcf --double-id --chr-set 31 --recode --out "670K_Equcab3_Mixedbreed_Bellone"

#4
#cat 670K_Equcab3_Belgian_Haflinger_Bellone.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab3_Belgian_Haflinger_Bellone.vcf --double-id --chr-set 31 --recode --out "670K_Equcab3_Belgian_Haflinger_Bellone"

#5
suffix=AHSdataset
plink --bfile "$suffix" --chr-set 34 no-y no-xy no-mt --allow-extra-chr --recode --out "$suffix"

#6
suffix=IR.paper.heredity
plink --bfile "$suffix" --chr-set 34 no-y no-xy no-mt --allow-extra-chr --recode --out "$suffix"

#7
#cat 70K_Equcab3_Bardigiano_Ablondi.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 70K_Equcab3_Bardigiano_Ablondi.vcf --double-id --chr-set 72 no-y no-xy no-mt --allow-extra-chr --recode --out "70K_Equcab3_Bardigiano_Ablondi"

#8
#cat 670K_Equcab2_AbagaBlack_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_AbagaBlack_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_AbagaBlack_Han"


#9
#cat 670K_Equcab2_BaichaIronHoof_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_BaichaIronHoof_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_BaichaIronHoof_Han"


#10
#cat 670K_Equcab2_Sanhe_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_Sanhe_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_Sanhe_Han"

#11
#cat 670K_Equcab2_Wushen_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_Wushen_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_Wushen_Han"


#12
#cat 670K_Equcab2_Wuzhumuqin_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_Wuzhumuqin_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_Wuzhumuqin_Han"

#13
#cat 70K_Equcab3_MongolianRacing_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 70K_Equcab3_MongolianRacing_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "70K_Equcab3_MongolianRacing_Hill"


#14
#cat 70K_Equcab3_ArabianHorse_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 70K_Equcab3_ArabianHorse_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "70K_Equcab3_ArabianHorse_Hill"

#15
#cat 670K_Equcab2_IrishDraught_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_IrishDraught_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_IrishDraught_Hill"

#16
#cat 670K_Equcab2_ConnemaraPony_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_ConnemaraPony_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_ConnemaraPony_Hill"

#17
#cat 670K_Equcab2_IrishSportHorse_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_IrishSportHorse_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_IrishSportHorse_Hill"

#-------
#18
#zcat QH_IlluminaSNP70_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf QH_IlluminaSNP70_EquCab2.vcf.gz --double-id --chr-set 31 no-xy no-mt --recode --out "QH_IlluminaSNP70_EquCab2"

#19
#zcat STDB_IlluminaSNP70_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf STDB_IlluminaSNP70_EquCab2.vcf.gz --double-id --chr-set 31 no-xy no-mt --recode --out "STDB_IlluminaSNP70_EquCab2"

#20
#zcat QH_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf QH_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "QH_IlluminaSNP50_EquCab2"


#21
#zcat STDB_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf STDB_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "STDB_IlluminaSNP50_EquCab2"


#22
#zcat TB_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf TB_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "TB_IlluminaSNP50_EquCab2"

#23
#zcat Belgian_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf Belgian_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "Belgian_IlluminaSNP50_EquCab2"

#24
#zcat Morgan_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf Morgan_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "Morgan_IlluminaSNP50_EquCab2"

#24
#zcat STDB_AffyMNEc670K_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf STDB_AffyMNEc670K_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --allow-extra-chr --recode --out "STDB_AffyMNEc670K_EquCab2"

#25
#zcat MNEc2M_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf MNEc2M_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --allow-extra-chr --recode --out "MNEc2M_EquCab2"


## make a backup copy of all newly derived files
mkdir -p backup_drived
cp *.{map,ped,log,nosex} backup/.
#############
## Move all ped/map files to the work_dir
cd $work_dir/.
cp backup_{original,drived}/*.{ped,map} $work_dir/.

cp backup_drived/Belgian_IlluminaSNP50_EquCab2.{ped,map} $work_dir/.
cp backup_drived/Morgan_IlluminaSNP50_EquCab2.{ped,map} $work_dir/.
cp backup_drived/QH_IlluminaSNP50_EquCab2.{ped,map} $work_dir/.
cp backup_drived/QH_IlluminaSNP70_EquCab2.{ped,map} $work_dir/.
cp backup_drived/STDB_AffyMNEc670K_EquCab2.{ped,map} $work_dir/.
cp backup_drived/STDB_IlluminaSNP50_EquCab2.{ped,map} $work_dir/.
cp backup_drived/STDB_IlluminaSNP70_EquCab2.{ped,map} $work_dir/.
cp backup_drived/TB_IlluminaSNP50_EquCab2.{ped,map} $work_dir/.
cp backup_drived/MNEc2M_EquCab2.{ped,map} $work_dir/.

############
#### download genomes
# Y chromosome: According to Terje Raudsepp, The Y chromosome assembly is still the one from 2018 Janecka et al, thus eMSYv3. The single copy assembly region is the one suitable for your SNP panel. Barbara Wallner would be a better contact for useful Y SNPs.
mkdir -p "$work_dir"/eMSY && cd "$work_dir"/eMSY ## download locally from Terje's email, upload to the drive, then rclone
rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs/eMSYv3.fasta $HOME/Horse_parentage_SNPs/eMSY/

# equCab2:
mkdir -p "$work_dir"/equCab2/download && cd "$work_dir"/equCab2/download
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab2/bigZips/chromFa.tar.gz' -O chromFa.tar.gz
tar xvzf chromFa.tar.gz
cd "$work_dir"
cat $(ls -1v equCab2/download/*.fa) > equCab2/equCab2_genome.fa
sed -i 's/>chr/>/' equCab2/equCab2_genome.fa
equCab2_ref="$work_dir"/equCab2/equCab2_genome.fa
cat $equCab2_ref | grep -v "^#" | awk -F">" '{if(!$1){if(NR!=1)print "##contig=<ID="ref",length="reflen">";ref=$2;reflen=0}else reflen+=length($1)}END{print "##contig=<ID="ref",length="reflen">";}' > equCab2/vcf_contigs.txt
equCab2_vcfContigs="$work_dir"/equCab2/vcf_contigs.txt

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/305/GCF_000002305.2_EquCab2.0/GCF_000002305.2_EquCab2.0_assembly_report.txt
grep -v "^#" GCF_000002305.2_EquCab2.0_assembly_report.txt | sed -e "s/\r//g" | awk 'BEGIN{FS=OFS="\t"}{print $10,$5}' | sed 's/^chr//' > equCab2_chrMap
equCab2_chrMap="$work_dir/equCab2/equCab2_chrMap"

# equCab3
mkdir -p "$work_dir"/equCab3/download && cd "$work_dir"/equCab3/download
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/equCab3/bigZips/equCab3.fa.gz' -O equCab3.fa.gz
gunzip equCab3.fa.gz
cd "$work_dir"
sed 's/>chr/>/' equCab3/download/equCab3.fa > equCab3/equCab3_genome.fa
equCab3_ref="$work_dir"/equCab3/equCab3_genome.fa

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/863/925/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_assembly_report.txt
grep -v "^#" GCF_002863925.1_EquCab3.0_assembly_report.txt | sed -e "s/\r//g" | awk 'BEGIN{FS=OFS="\t"}{print
$10,$5}' | sed 's/^chr//' > equCab3_chrMap
equCab3_chrMap="$work_dir/equCab3/equCab3_chrMap"


##### array data
## 70k
mkdir -p $work_dir/SNParrays/70k && cd $work_dir/SNParrays/70k
# GGP Equine Manifest File (CSV Format)
# https://sapac.illumina.com/products/by-type/microarray-kits/ggp-equine.html
# https://sapac.support.illumina.com/downloads/geneseek-ggp-equine-product-files.html
wget https://sapac.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/geneseek-ggp/geneseek-ggp-equine-manifest-file-csv.zip
unzip geneseek-ggp-equine-manifest-file-csv.zip
tail -n+8 GGP_Equine.csv | cut -d, -f10 | sort | uniq -c
head -n8 GGP_Equine.csv | tail -n1 | awk 'BEGIN{FS=",";OFS="\t"}{print $2,$1,$10,$11,$3,$4,"Allele_A","Allele_B"}' > GGP_Equine.simple.csv
tail -n+9 GGP_Equine.csv | head -n-24 | awk 'BEGIN{FS=",";OFS="\t"}{print $2,$1,$10,$11,$3,$4,substr($4, 2, 1),substr($4, 4, 1)}' >> GGP_Equine.simple.csv
#cat GGP_Equine.simple.csv | cut -f3 | sort | uniq -c  ## 31 X(3409) Y(2)
## QC
cat GGP_Equine.simple.csv | awk -F"\t" '{print $3"_"$4}' | sort | uniq -c | awk '{if($1>1)print $2}' > dup_loci ## 55 SNPs (different IDs but same position)
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$3"_"$4]==1){a[$3"_"$4]++;b[$3"_"$4]=$1;}else if(a[$3"_"$4]>1)b[$3"_"$4]=b[$3"_"$4] FS $1;}END{for(i in b)print i,b[i]}' dup_loci GGP_Equine.simple.csv > dup_snps


# tranform the array manifest into VCF & generate ALT/REF map
bash "$scripts"/manifest_to_ref.sh "GGP_Equine.simple.csv" "$equCab2_vcfContigs" "$equCab2_ref" "GGP_Equine"
arr70k="$work_dir/SNParrays/70k/GGP_Equine.alt_ref_comp.tab"

## get the minimal required data from the annotation file including Allele A and B in both senses
cat GGP_Equine.simple.csv | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$7,$8;}' > ann_70_temp.tab
tail -n+2 ann_70_temp.tab | awk 'BEGIN{FS=OFS="\t"}{print $5,$6}' | tr 'TCGA' 'AGCT' > ann_70.chip_alleles_oppStrand
echo -e "$(head -n1 ann_70_temp.tab)\tAllele_A_opp\tAllele_B_opp" > ann_70.tab
paste <(tail -n+2 ann_70_temp.tab) ann_70.chip_alleles_oppStrand >> ann_70.tab
rm ann_70_temp.tab ann_70.chip_alleles_oppStrand

# Molly files
#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/SNP70.NCBI.remap.csv.gz
wget https://www.animalgenome.org/repository/pub/UMN2018.1003/SNP70.unique_remap.FINAL.csv.gz
gunzip *.gz
#tail -n+2 SNP70.NCBI.remap.csv | cut -d"," -f4 | sort | uniq -c  > eq2.chr
#tail -n+2 SNP70.NCBI.remap.csv | cut -d"," -f5 | sort | uniq -c  > eq3.chr

## merge the annotation to the map file which has the reference alleles but not for the SNPs that failed to remap
cat SNP70.unique_remap.FINAL.csv | sed 's/chrUn_ref|/Un_/' | sed 's/\.1|,/v1,/' | sed 's/,chr/,/g' | tr ',' '\t' > SNP70_remap.tab
echo -e "$(head -n1 ann_70.tab)\t$(head -n1 SNP70_remap.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$6,$8,$9,$10}')" > map_70.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$5"_"$7]=$1 FS $3 FS $6 FS $8 FS $9 FS $10;next}{if(a[$3"_"$4])print $0,a[$3"_"$4];else print $0,0,0,0,0,0,0;}' SNP70_remap.tab <(tail -n+2 ann_70.tab) >> map_70.tab

## (Optional) add reference alleles for all SNPs in Equcab2 (except Un2 chromosome bc we do not have it sequence)
echo -e "$(head -n1 map_70.tab)\tSNP_id\tALT\tREF\tcALT\tcREF" > map_70_ext.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next}{if(a[$1])print $0,a[$1];else print $0,0,0,0,0,0;}' "$arr70k" <(tail -n+2 map_70.tab) >> map_70_ext.tab ## Note that column $10 should match column $17 except for the unmapped or Un2 SNPs
#cat map_70_ext.tab | awk -F"\t" '{if($10!=$17)print}' > QC.log
#cat QC.log | awk -F"\t" '{if($9=="0" && $15!="0")print}' | wc -l
#cat QC.log | awk -F"\t" '{if($9!="0" && $15=="0")print}' | wc -l
#cat QC.log | awk -F"\t" '{if($9=="0" && $15=="0")print}' | wc -l

str=$(head -n1 map_70_ext.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$9,$3,$4}')
end=$(head -n1 map_70_ext.tab | awk 'BEGIN{FS=OFS="\t"}{print $11,$12,$13,$14}')
echo -e "$str\tEC2_REF\tEC2_ALT\t$end" > map_70_complete.tab
tail -n+2 map_70_ext.tab | awk 'BEGIN{FS=OFS="\t"}{x=$1 FS $2 FS $9 FS $3 FS $4;y=$11 FS $12 FS $13 FS $14;\
if($10==$17)print x,$17,$16,y;\
else if($9=="0" && $15!="0")print x,$17,$16,y;\
else if($9!="0" && $15=="0"){\
     if($10==$5)print x,$5,$6,y;\
     else if($10==$6)print x,$6,$5,y;\
     else if($10==$7)print x,$7,$8,y;\
     else if($10==$8)print x,$8,$7,y;\
     else print x,$10,"0",y;}\
else if($9=="0" && $15=="0")print $5,$6;}' >> map_70_complete.tab

# find duplicate SNPs that should be excluded
tail -n+2 map_70_complete.tab | awk -F"\t" '{if($3 && $1!=$3)print $1}' > dup_snps_exclude
cat dup_snps_exclude | grep -vFwf - dup_snps | awk -F"\t" '{print $2}' > dup_snps_exclude.temp
cat dup_snps_exclude.temp >> dup_snps_exclude; rm dup_snps_exclude.temp;
# remove duplicate SNPs
cp map_70_complete.tab map_70_complete.tab.temp
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(!a[$1])print}' dup_snps_exclude map_70_complete.tab.temp > map_70_complete.tab
rm map_70_complete.tab.temp
map70="$work_dir/SNParrays/70k/map_70_complete.tab"

## 50k
mkdir -p $work_dir/SNParrays/50k && cd $work_dir/SNParrays/50k
# Molly files
wget https://www.animalgenome.org/repository/pub/UMN2018.1003/SNP50.NCBI.remap.csv.gz
gunzip *.gz
tail -n+2 SNP50.NCBI.remap.csv | cut -d"," -f4 | sort | uniq -c  > eq2.chr
tail -n+2 SNP50.NCBI.remap.csv | cut -d"," -f5 | sort | uniq -c  > eq3.chr

## 670k
mkdir -p $work_dir/SNParrays/670k && cd $work_dir/SNParrays/670k
## https://www.thermofisher.com/order/catalog/product/550583
## http://media.affymetrix.com/analysis/downloads/lf/genotyping/Axiom_MNEc670/
wget http://media.affymetrix.com/analysis/downloads/lf/genotyping/Axiom_MNEc670/r2/Axiom_MNEc670_Annotation.r2.csv.zip
unzip Axiom_MNEc670_Annotation.r2.csv.zip
sed -i 's/","/"|"/g' Axiom_MNEc670_Annotation.r2.csv
cat Axiom_MNEc670_Annotation.r2.csv | grep -v "^#" | awk 'BEGIN{FS="|";OFS="\t"}{print $1,$2,$4,$5,$6,$7,$8,$9,$11}' | sed 's/"//g' > Axiom_MNEc670_Annotation.r2.simple.csv
#cat Axiom_MNEc670_Annotation.r2.simple.csv | cut -f3 | sort | uniq -c  ## 31 Un(21468) Un2(9395) X(28017) Y(1)
## QC
cat Axiom_MNEc670_Annotation.r2.simple.csv | awk -F"\t" '{print $3"_"$4}' | sort | uniq -c | awk '{if($1>1)print $2}' > dup_loci  ## 0 SNPs (different IDs but same position)


# tranform the array manifest into VCF & generate ALT/REF map
grep -vw "Un2" Axiom_MNEc670_Annotation.r2.simple.csv > Axiom_MNEc670_Annotation.r2.simple_excludeUN2.csv
bash "$scripts"/manifest_to_ref.sh "Axiom_MNEc670_Annotation.r2.simple_excludeUN2.csv" "$equCab2_vcfContigs" "$equCab2_ref" "Axiom_MNEc670_r2"
arr670k="$work_dir/SNParrays/670k/Axiom_MNEc670_r2.alt_ref_comp.tab"


## get the minimal required data from the annotation file including Allele A and B in both senses
cat Axiom_MNEc670_Annotation.r2.simple.csv | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$7,$8;}' > ann_670_temp.tab
tail -n+2 ann_670_temp.tab | awk 'BEGIN{FS=OFS="\t"}{print $5,$6}' | tr 'TCGA' 'AGCT' > ann_670.chip_alleles_oppStrand
echo -e "$(head -n1 ann_670_temp.tab)\tAllele_A_opp\tAllele_B_opp" > ann_670.tab
paste <(tail -n+2 ann_670_temp.tab) ann_670.chip_alleles_oppStrand >> ann_670.tab
rm ann_670_temp.tab ann_670.chip_alleles_oppStrand

# Molly files
#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc2M.blast.EquCab3.txt.gz
#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc2M.NCBI.remap.csv.gz
#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc2M.probe_sequence.fa.gz
#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc2M.remap_failure.csv.gz
#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc2M.unique_remap.FINAL.csv.gz
#wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc670k.remap_failure.csv.gz
wget https://www.animalgenome.org/repository/pub/UMN2018.1003/MNEc670k.unique_remap.FINAL.csv.gz
gunzip *.gz
#tail -n+2 MNEc2M.NCBI.remap.csv | cut -d"," -f4 | sort | uniq -c  > eq2.chr
#tail -n+2 MNEc2M.NCBI.remap.csv | cut -d"," -f5 | sort | uniq -c  > eq3.chr

## merge the annotation to the map file which has the reference alleles but not for the SNPs that failed to remap
cat MNEc670k.unique_remap.FINAL.csv | sed 's/chrUn_ref|/Un_/' | sed 's/\.1|,/v1,/' | sed 's/,chr/,/g' | tr ',' '\t' > MNEc670k_remap.tab
echo -e "$(head -n1 ann_670.tab)\t$(head -n1 MNEc670k_remap.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$6,$8,$9,$10}')" > map_670.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$5"_"$7]=$1 FS $3 FS $6 FS $8 FS $9 FS $10;next}{if(a[$3"_"$4])print $0,a[$3"_"$4];else print $0,0,0,0,0,0,0;}' MNEc670k_remap.tab <(tail -n+2 ann_670.tab) >> map_670.tab

## (Optional) add reference alleles for all SNPs in Equcab2 (except Un2 chromosome bc we do not have it sequence)
echo -e "$(head -n1 map_670.tab)\tSNP_id\tALT\tREF\tcALT\tcREF" > map_670_ext.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next}{if(a[$1])print $0,a[$1];else print $0,0,0,0,0,0;}' "$arr670k" <(tail -n+2 map_670.tab) >> map_670_ext.tab ## Note that column $10 should match column $17 except for the unmapped or Un2 SNPs
#cat map_670_ext.tab | awk -F"\t" '{if($10!=$17)print}' > QC.log
#cat QC.log | awk -F"\t" '{if($9=="0" && $15!="0")print}' | wc -l
#cat QC.log | awk -F"\t" '{if($9!="0" && $15=="0")print}' | wc -l
#cat QC.log | awk -F"\t" '{if($9=="0" && $15=="0")print}' | wc -l

str=$(head -n1 map_670_ext.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$9,$3,$4}')
end=$(head -n1 map_670_ext.tab | awk 'BEGIN{FS=OFS="\t"}{print $11,$12,$13,$14}')
echo -e "$str\tEC2_REF\tEC2_ALT\t$end" > map_670_complete.tab
tail -n+2 map_670_ext.tab | awk 'BEGIN{FS=OFS="\t"}{x=$1 FS $2 FS $9 FS $3 FS $4;y=$11 FS $12 FS $13 FS $14;\
if($10==$17)print x,$17,$16,y;\
else if($9=="0" && $15!="0")print x,$17,$16,y;\
else if($9!="0" && $15=="0"){\
     if($10==$5)print x,$5,$6,y;\
     else if($10==$6)print x,$6,$5,y;\
     else if($10==$7)print x,$7,$8,y;\
     else if($10==$8)print x,$8,$7,y;\
     else print x,$10,"0",y;}\
else if($9=="0" && $15=="0")print $5,$6;}' >> map_670_complete.tab
map670="$work_dir/SNParrays/670k/map_670_complete.tab"

## compare to the Equcab3 manifest file
wget http://media.affymetrix.com/analysis/downloads/lf/genotyping/Axiom_MNEc670/r3/Axiom_MNEc670.na35.r3.a2.annot.csv.zip
unzip Axiom_MNEc670.na35.r3.a2.annot.csv.zip
sed -i 's/","/"|"/g' Axiom_MNEc670.na35.r3.a2.annot.csv
cat Axiom_MNEc670.na35.r3.a2.annot.csv | grep -v "^#" | awk 'BEGIN{FS="|";OFS="\t"}{print $1,$2,$5,$6,$8,$11,$12,$13,$14,$15,$33,$34,$36}' | sed 's/"//g' > Axiom_MNEc670.na35.r3.a2.annot.simple.csv
#cat Axiom_MNEc670.na35.r3.a2.annot.simple.csv | cut -f3 | sort | uniq -c  ## 0(39127) 31 X(27327) Y(1) + many unplaced contigs

cat map_670_complete.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$8,$9,$10,$11}' > molly.tab
cat Axiom_MNEc670.na35.r3.a2.annot.simple.csv | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$4,$9,$10}' > Axiom.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next}{if(a[$1])print $0,a[$1];else print $0,0,0,0,0,0;}' molly.tab Axiom.tab > molly_Axiom_compare.tab
cat molly_Axiom_compare.tab | awk -F"\t" '{if($6=="0")print}' ## should be nothing (i.e. both files have same SNPs)
cat molly_Axiom_compare.tab | awk -F"\t" '{if($2=="---" && $7=="0")print}' | wc -l ## 26959  failed to remap in both
cat molly_Axiom_compare.tab | awk -F"\t" '{if($2!="---" && $7=="0")print}' | wc -l ## 1919   mapped by Axiom only
cat molly_Axiom_compare.tab | awk -F"\t" '{if($2=="---" && $7!="0")print}' | wc -l ## 12168  mapped by Molly only
cat molly_Axiom_compare.tab | awk -F"\t" '{if($2!="---" && $7!="0" && $3!=$8)print}' | wc -l ## 94  different mapping

## 80k
## There is no official manifest for the 80k array
## instead, I will use SNP_ids from all the map files I have for this array, then extract those SNPs from 70k & 670k maps
## Please note that this code assumes that the map files have been downloaded & passed the first 3 fixing steps
mkdir -p $work_dir/SNParrays/80k && cd $work_dir/SNParrays/80k
cat $work_dir/Clyde_70K.map \
$work_dir/EQ75v03_EquCab2_Falabella_DiazS.map \
$work_dir/70K_Equcab2_Maremmano_Capomaccio.map \
$work_dir/GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map \
$work_dir/80K_Equcab3_App_Knab_Bellone.map \
$work_dir/80K_Equcab3_MultipleBreeds_vit.map | awk -F"\t" '{print $2}' | sort | uniq > candidate_80snps

head -n1 $map70 > map_80_complete.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$1])print $0;}' candidate_80snps $map70 >> map_80_complete.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$2]){t=$1; $1=$2; $2=t;print;}}' candidate_80snps $map670 >> map_80_complete.tab
map80="$work_dir/SNParrays/80k/map_80_complete.tab"  ## 71947

## QC
cat map_80_complete.tab | awk -F"\t" '{print $4"_"$5}' | sort | uniq -c | awk '{if($1>1)print $2}' > dup_loci ## empty bc the 55 duplicate SNPs were excluded from map70
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$4"_"$5]==1){a[$4"_"$5]++;b[$4"_"$5]=$1;}else if(a[$4"_"$5]>1)b[$4"_"$5]=b[$4"_"$5] FS $1;}END{for(i in b)print i,b[i]}' dup_loci map_80_complete.tab > dup_snps


## https://github.com/UMN-EGGL/MNEc2to3
mkdir $work_dir/SNParrays/MNEc2to3 && cd $work_dir/SNParrays/MNEc2to3
wget https://raw.githubusercontent.com/UMN-EGGL/MNEc2to3/master/data/MNEc2M.probe_blast_counts.csv.gz
tail -n+2 MNEc2M.probe_blast_counts.csv | awk -F"," '{if($11=="True")print $5}' | sort | uniq -c
# Note: There are multiple records with the same coordinates but different SNP ID e.g.
# 24940,MNEc.2.1.29763.UKUL2,Affx-101643281,AX-104566380,chr1,29763,T,C,C,T,True,True,chr1.29763,1,1,1
# 24941,MNEc.2.1.29763.UKUL2,Affx-101643281,AX-103622354,chr1,29763,T,C,C,T,True,True,chr1.29763,1,1,1

## map between 70k and 670k
cd $work_dir/SNParrays
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1"_"$2]=$3 FS $4 FS $5;next}{if(a[$1"_"$2])print $1,$2,a[$1"_"$2],$3,$4,$5}'  <(grep -v "^#" 70k/GGP_Equine.vcf) <(grep -v "^#" 670k/Axiom_MNEc670_r2.vcf) > 70k_670k_map.txt ## 62276
#cat 70k_670k_map.txt | awk -F"\t" '{if($4!=$7)print}'
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8"_"$9]=$1 FS $10 FS $11;next}{if(a[$8"_"$9])print $8"_"$9,a[$8"_"$9],$1,$10,$11}'  $map70 $map670 > 70k_670k_map.txt ## 59051
####################
### Fix input map/ped files format
## FIX: 1a. fix the files with abonormal end characters
for f in *.{map,ped};do
 sed -i -e "s/\r//g" $f;
done

## FIX: 1b. unify the file separator
for f in *.{map,ped};do if [ -f $f ];then cp $f $f.temp; cat $f.temp | grep -v "^#" | tr ' ' '\t' > $f;rm $f.temp;fi;done

####################
## CHECK: check map files format
for f in *.map;do echo $(wc -l $f) "--" $(cat -v $f | grep -v "^#" | head -n1) "--" $(head -n1 $f | awk -F"\t" '{print NF}');done > checkMapFormat.txt

## Results:
## 80K_Equcab3_MultipleBreeds_vit.map has 3 columns only

## Fix the missing columns in the map files
cp 80K_Equcab3_MultipleBreeds_vit.map 80K_Equcab3_MultipleBreeds_vit.map.temp
cat 80K_Equcab3_MultipleBreeds_vit.map.temp | awk 'BEGIN{FS=OFS="\t"}{$3=0 FS $3;print}' > 80K_Equcab3_MultipleBreeds_vit.map
rm 80K_Equcab3_MultipleBreeds_vit.map.temp

####################
## CHECK: check ped files format
for f in *.ped;do echo $f $(head -n1 $f | cut -f1-10);done > checkPedFormat.txt

## Results:
## Some ped files are Irregularly-formatted PLINK text files
##    - 670K_Equcab2_Exmoor_Lindgren-Velie.ped
##    - 670K_Equcab2_Icelandic_Lindgren-Velie.ped
## Some ped files are numberically encoded
##    - 670K_Equcab3_FranchesMontagnes_Gmel.ped
## Some ped files are using A/B call
##    - Clyde_70K.ped

## Fix the Irregularly-formatted PLINK text files
mv 670K_Equcab2_Exmoor_Lindgren-Velie.ped 670K_Equcab2_Exmoor_Lindgren-Velie.ped.toBedeleted
cat 670K_Equcab2_Exmoor_Lindgren-Velie.ped.toBedeleted | awk 'BEGIN{FS=OFS="\t"}{$1=1 FS $1 FS 0 FS 0 FS 0 FS "-9";print $0}' > 670K_Equcab2_Exmoor_Lindgren-Velie.ped
rm 670K_Equcab2_Exmoor_Lindgren-Velie.ped.toBedeleted

mv 670K_Equcab2_Icelandic_Lindgren-Velie.ped 670K_Equcab2_Icelandic_Lindgren-Velie.ped.toBedeleted
cat 670K_Equcab2_Icelandic_Lindgren-Velie.ped.toBedeleted | awk 'BEGIN{FS=OFS="\t"}{$1=1 FS $1 FS 0 FS 0 FS 0 FS "-9";print $0}' > 670K_Equcab2_Icelandic_Lindgren-Velie.ped
rm 670K_Equcab2_Icelandic_Lindgren-Velie.ped.toBedeleted

## Fix the numberically encoded ped files
mv 670K_Equcab3_FranchesMontagnes_Gmel.map 670K_Equcab3_FranchesMontagnes_Gmel.temp.map
mv 670K_Equcab3_FranchesMontagnes_Gmel.ped 670K_Equcab3_FranchesMontagnes_Gmel.temp.ped
plink --file 670K_Equcab3_FranchesMontagnes_Gmel.temp --chr-set 36 --recode --alleleACGT --out 670K_Equcab3_FranchesMontagnes_Gmel
rm 670K_Equcab3_FranchesMontagnes_Gmel.temp.* 670K_Equcab3_FranchesMontagnes_Gmel.{log,nosex}
mv 670K_Equcab3_FranchesMontagnes_Gmel.ped 670K_Equcab3_FranchesMontagnes_Gmel.temp.ped
cat 670K_Equcab3_FranchesMontagnes_Gmel.temp.ped | grep -v "^#" | tr ' ' '\t' > 670K_Equcab3_FranchesMontagnes_Gmel.ped
rm 670K_Equcab3_FranchesMontagnes_Gmel.temp.ped

## Fix: Remove datasets using A/B calls
rm Clyde_70K.{map,ped}

#################
##CHECK for the No of chromosomes per map file
#cat 80K_Equcab3_App_Knab_Bellone.map  | grep -v "^#" | awk '{print $1}' | sort | uniq -c | sort -k2,2n
for f in *.map;do echo $f; cat $f  | grep -v "^#" | awk '{print $1}' | sort | uniq -c | sort -k2,2n;done > checkMapChrNo.txt

## Note2: chromosome IDs
# 80K_Equcab3_App_Knab_Bellone.map has PJAA* chromosomes. These need to be mapped to the UCSC format
# 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map has CM* chr. These need to be mapped to the UCSC format
# Several map files have inconsistant encoding for X, Y , Un, and Un2 chromosomes
# experiments remapped to Equcab3 are using chr & pos = 0 for unmapped features
  # 1. 670K_Equcab3_FranchesMontagnes_Gmel.map
  # 2. 80K_Equcab3_MultipleBreeds_vit.map
  # 3. 80K_Equcab3_App_Knab_Bellone.map
## default annotation in Plink
#--horse = --chr-set 31 no-xy no-mt
# * Note that, when there are n autosome pairs, the X chromosome is assigned numeric code n+1, Y is n+2, XY (pseudo-autosomal region of X) is n+3, and MT (mitochondria) is n+4.
# * Normally, PLINK reports an error if the input data contains unrecognized chromosome codes (such as unplaced contigs). If none of the additional codes start with a digit, you can permit them with the --allow-extra-chr flag. (These contigs are ignored by most analyses which skip unplaced regions.)

## Fix chromosome symbols
# the only file that has CM* chromosomes
cp 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map.temp
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$1;next}{if(a[$1])$1=a[$1];print $0;}' "$equCab2_chrMap" 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map.temp > 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map
rm 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map.temp

# the only file that has PJAA* chromosomes
cp 80K_Equcab3_App_Knab_Bellone.map 80K_Equcab3_App_Knab_Bellone.map.temp
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$1;next}{if(a[$1])$1=a[$1];print $0;}' "$equCab3_chrMap" 80K_Equcab3_App_Knab_Bellone.map.temp | sed 's/^32/X/' > 80K_Equcab3_App_Knab_Bellone.map
rm 80K_Equcab3_App_Knab_Bellone.map.temp

# the files that have chrUn and chrUn2 instead of Un and Un2
cp STDB_AffyMNEc670K_EquCab2.map STDB_AffyMNEc670K_EquCab2.map.temp
cat STDB_AffyMNEc670K_EquCab2.map.temp | sed 's/^chrUn/Un/; s/^chrUn2/Un2/; s/^32/X/;' > STDB_AffyMNEc670K_EquCab2.map
rm STDB_AffyMNEc670K_EquCab2.map.temp

cp MNEc2M_EquCab2.map MNEc2M_EquCab2.map.temp
cat MNEc2M_EquCab2.map.temp | sed 's/^chrUn/Un/; s/^chrUn2/Un2/; s/^32/X/;' > MNEc2M_EquCab2.map
rm MNEc2M_EquCab2.map.temp


# the files missing X (and Y) chromosome (X encoded as 32 & Y encoded as 33)
for f in \
"50K_Equcab2_FranchesMontagnes_Gmel.map" \
"QH_IlluminaSNP50_EquCab2.map" \
"STDB_IlluminaSNP50_EquCab2.map" \
"TB_IlluminaSNP50_EquCab2.map" \
"Belgian_IlluminaSNP50_EquCab2.map" \
"Morgan_IlluminaSNP50_EquCab2.map" \
"QH_IlluminaSNP70_EquCab2.map" \
"STDB_IlluminaSNP70_EquCab2.map" \
"EQN65KV02_EquCab2_ArgentineanCreole_GiovambattistaG-DiazS.map" \
"EQ75v03_EquCab2_Falabella_DiazS.map" \
"70K_Equcab2_Maremmano_Capomaccio.map" \
"670K_Equcab3_FranchesMontagnes_Gmel.map" \
"670K_Equcab3_Belgian_Haflinger_Bellone.map" \
"670K_Equcab3_Mixedbreed_Bellone.map";do
cp $f $f.temp; cat $f.temp | sed 's/^32/X/; s/^33/Y/' > $f; rm $f.temp;done

# the files missing X (and Y) and Un/Un2 chromosomes (Un, Un2, X, and Y are encoded as 33,34,35,36)
for f in \
"670k_Equcab2_Arabian_Patterson.map" \
"670K_Equcab2_Exmoor_Lindgren-Velie.map" \
"670K_Equcab2_Icelandic_Lindgren-Velie.map" \
"670K_EquCab2_SWB_Mikko.map" \
"58TBD_52KRD_24PAR_272724_Unpruned.map" \
"AHSdataset.map" \
"IR.paper.heredity.map" \
"670K_EquCab2_Standardbreds_Bellone.map";do
cp $f $f.temp; cat $f.temp | sed 's/^33/Un/; s/^34/Un2/; s/^35/X/; s/^36/Y/' > $f; rm $f.temp;done


##CHECK 3. check for SNP ids
> unrecognized_SNP_ids.txt
#arr="70"
for input in \
"70K_Equcab2_Haflinger_Bellone.map" \
"70K_Equcab2_PeurtoRicanPasoFino_Bellone.map" \
"GGP65K_Equcab_haflinger_agrotis.map" \
"GGP65K_Equcab_Italian-heavy-draft-horse_agrotis.map" \
"GGP65K_Equcab_noriker_agrotis.map" \
"EQN65KV02_EquCab2_ArgentineanCreole_GiovambattistaG-DiazS.map" \
"70K_Equcab3_ArabianHorse_Hill.map" \
"70K_Equcab3_Bardigiano_Ablondi.map" \
"70K_Equcab3_MongolianRacing_Hill.map" \
"QH_IlluminaSNP70_EquCab2.map" \
"STDB_IlluminaSNP70_EquCab2.map";do
echo $input "--" $(comm -23 <(awk -F"\t" '{print $2}' "$input" | sort) <(tail -n+2 "$map70" | awk -F"\t" '{print $1}' | sort) | wc -l);
done >> unrecognized_SNP_ids.txt


#arr="670"
for input in \
"670k_Equcab2_Arabian_Patterson.map" \
"670K_Equcab2_Exmoor_Lindgren-Velie.map" \
"670K_Equcab2_Icelandic_Lindgren-Velie.map" \
"670K_EquCab2_SWB_Mikko.map" \
"58TBD_52KRD_24PAR_272724_Unpruned.map" \
"AHSdataset.map" \
"IR.paper.heredity.map" \
"670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map" \
"670K_EquCab2_Standardbreds_Bellone.map" \
"670K_Equcab2_AbagaBlack_Han.map" \
"670K_Equcab2_BaichaIronHoof_Han.map" \
"670K_Equcab2_ConnemaraPony_Hill.map" \
"670K_Equcab2_IrishDraught_Hill.map" \
"670K_Equcab2_IrishSportHorse_Hill.map" \
"670K_Equcab2_Sanhe_Han.map" \
"670K_Equcab2_Wushen_Han.map" \
"670K_Equcab2_Wuzhumuqin_Han.map" \
"670K_Equcab3_FranchesMontagnes_Gmel.map" \
"670K_Equcab3_Mixedbreed_Bellone.map" \
"670K_Equcab3_Belgian_Haflinger_Bellone.map" \
"STDB_AffyMNEc670K_EquCab2.map" \
"MNEc2M_EquCab2.map";do
echo $input "--" $(comm -23 <(awk -F"\t" '{print $2}' "$input" | sort) <(tail -n+2 "$map670" | awk -F"\t" '{print $1}' | sort) | wc -l);
done >> unrecognized_SNP_ids.txt

#arr="80"
for input in \
"EQ75v03_EquCab2_Falabella_DiazS.map" \
"70K_Equcab2_Maremmano_Capomaccio.map" \
"GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map" \
"80K_Equcab3_App_Knab_Bellone.map" \
"80K_Equcab3_MultipleBreeds_vit.map";do
echo $input "--" $(comm -23 <(awk -F"\t" '{print $2}' "$input" | sort) <(tail -n+2 "$map80" | awk -F"\t" '{print $1}' | sort) | wc -l);
done >> unrecognized_SNP_ids.txt

## Results: unrecognized SNP ids
# All 70k & 80k will show the duplicate SNPs as unrecognized SNP ids bc we removed from map70
# Some datasets has many unrecognized SNP ids due to having extra-chromosomes or mapping to unrecognized Equcab3 position
#670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map -- 639909
#670K_Equcab3_Belgian_Haflinger_Bellone.map -- 619466
#GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map -- 1876 ==> has 275 duplicate SNPs (labeled *_dup)
#80K_Equcab3_App_Knab_Bellone.map -- 71947
#80K_Equcab3_MultipleBreeds_vit.map -- 16143
# postpone fixing unrecognized SNP ids after removing duplicates


##CHECK 4. check for duplicate SNP positions (This step edit all ids to make them unique)
for f in *.map;do
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($1)a[$1"_"$4]=$2;next}{if(a[$1"_"$4])a[$1"_"$4]=0;else print $2;}' "$f" "$f" > $f.exclude.lst
x=$(cat $f.exclude.lst | wc -l)
if [ $x -gt 0 ];then
  wc -l $f.exclude.lst
fi;done > dupSNPs.txt
rm *.exclude.lst

## Results:
## All 70k and 80k have built in 55 SNP duplicates
## other map files have their own duplicates


## FIX deduplicate SNP positions
# 1. change SNP IDs to be unique
for f in *.map;do
cp $f $f.temp; cat $f.temp | awk 'BEGIN{FS=OFS="\t"}{print $1,$2"."NR,$3,$4}' > $f; rm $f.temp;
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($1)a[$1"_"$4]=$2;next}{if(a[$1"_"$4])a[$1"_"$4]=0;else print $2;}' "$f" "$f" > $f.exclude.lst
done

# 2. exclude duplicated SNPs (This command excludes the unmapped SNPs as well)
for f in *.map;do
x=$(cat $f.exclude.lst | wc -l)
if [ $x -gt 0 ];then
  echo $f
  t=${f%.map}.ped
  dupf=${f%.map}.dup.map
  dupt=${f%.map}.dup.ped
  mv $f $dupf; mv $t $dupt;
  plink --file "${f%.map}.dup" --chr-set 31 no-xy --allow-extra-chr --recode --exclude $f.exclude.lst --out "${f%.map}"
  cp $f $f.temp; cat $f.temp | sed 's/^32/X/; s/^33/Y/' > $f; rm $f.temp;
  rm $dupf $dupt;
fi;done &> dedup.txt
rm *.exclude.lst
rm *.{hh,log,nosex}

##########
## Fix: remap (if needed) & fix the bad SNP ids & use 670 IDs only
mkdir -p tempbackup
cp *.{map,ped} tempbackup/.
##asm="2"
##arr="70"
#for f in \
#"70K_Equcab2_Haflinger_Bellone.map" \
#"70K_Equcab2_PeurtoRicanPasoFino_Bellone.map" \
#"GGP65K_Equcab_haflinger_agrotis.map" \
#"GGP65K_Equcab_Italian-heavy-draft-horse_agrotis.map" \
#"GGP65K_Equcab_noriker_agrotis.map" \
#"EQN65KV02_EquCab2_ArgentineanCreole_GiovambattistaG-DiazS.map";do
#cp $f $f.temp;
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$4"_"$5]=$8 FS $1 FS 0 FS $9;next}{if(a[$1"_"$4])print a[$1"_"$4];else print "0",$2,"0","0";}' "$map70" $f.temp > $f
#rm $f.temp;done

##asm="2"
##arr="80"
#for f in \
#"EQ75v03_EquCab2_Falabella_DiazS.map" \
#"70K_Equcab2_Maremmano_Capomaccio.map" \
#"GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map";do
#cp $f $f.temp;
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$4"_"$5]=$8 FS $1 FS 0 FS $9;next}{if(a[$1"_"$4])print a[$1"_"$4];else print "0",$2,"0","0";}' "$map80" $f.temp > $f
#rm $f.temp;done

##asm="2"
##arr="670"
#for f in \
#"670k_Equcab2_Arabian_Patterson.map" \
#"670K_Equcab2_Exmoor_Lindgren-Velie.map" \
#"670K_Equcab2_Icelandic_Lindgren-Velie.map" \
#"670K_EquCab2_SWB_Mikko.map" \
#"58TBD_52KRD_24PAR_272724_Unpruned.map" \
#"AHSdataset.map" \
#"IR.paper.heredity.map" \
#"670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map" \
#"670K_EquCab2_Standardbreds_Bellone.map" \
#"670K_Equcab2_AbagaBlack_Han.map" \
#"670K_Equcab2_BaichaIronHoof_Han.map" \
#"670K_Equcab2_ConnemaraPony_Hill.map" \
#"670K_Equcab2_IrishDraught_Hill.map" \
#"670K_Equcab2_IrishSportHorse_Hill.map" \
#"670K_Equcab2_Sanhe_Han.map" \
#"670K_Equcab2_Wushen_Han.map" \
#"670K_Equcab2_Wuzhumuqin_Han.map";do
#cp $f $f.temp;
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$4"_"$5]=$8 FS $1 FS 0 FS $9;next}{if(a[$1"_"$4])print a[$1"_"$4];else print "0",$2,"0","0";}' "$map670" $f.temp > $f
#rm $f.temp;done

##asm="3"
##arr="80"
#for f in \
#"80K_Equcab3_App_Knab_Bellone.map" \
#"80K_Equcab3_MultipleBreeds_vit.map" \
#"70K_Equcab3_ArabianHorse_Hill.map" \
#"70K_Equcab3_Bardigiano_Ablondi.map" \
#"70K_Equcab3_MongolianRacing_Hill.map";do
#cp $f $f.temp;
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8"_"$9]=$8 FS $1 FS 0 FS $9;next}{if(a[$1"_"$4])print a[$1"_"$4];else print "0",$2,"0","0";}' "$map80" $f.temp > $f
#rm $f.temp;done

##asm="3"
##arr="670"
#for f in \
#"670K_Equcab3_FranchesMontagnes_Gmel.map" \
#"670K_Equcab3_Mixedbreed_Bellone.map" \
#"670K_Equcab3_Belgian_Haflinger_Bellone.map";do
#cp $f $f.temp;
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8"_"$9]=$8 FS $1 FS 0 FS $9;next}{if(a[$1"_"$4])print a[$1"_"$4];else print "0",$2,"0","0";}' "$map670" $f.temp > $f
#rm $f.temp;done

for f in *.map;do
  equ2=$(comm -12 <(cat $f | awk '{print $1"."$4}' | sort) <(cat $map670 |  awk '{print $4"."$5}' | sort) | wc -l)
  equ3=$(comm -12 <(cat $f | awk '{print $1"."$4}' | sort) <(cat $map670 |  awk '{print $8"."$9}' | sort) | wc -l)
  cp $f $f.temp;
  if [ $equ2 -gt $equ3 ];then
    echo "equ2=$equ2 & equ3=$equ3 ... $f is equ2"
    awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$4"_"$5]=$8 FS $1 FS 0 FS $9;next}{if(a[$1"_"$4])print a[$1"_"$4];else print "0",$2,"0","0";}' "$map670" $f.temp > $f
  else
    echo "equ2=$equ2 & equ3=$equ3 ... $f is equ3"
    awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8"_"$9]=$8 FS $1 FS 0 FS $9;next}{if(a[$1"_"$4])print a[$1"_"$4];else print "0",$2,"0","0";}' "$map670" $f.temp > $f
  fi
done
rm *.map.temp

## Check the no of unmapped SNPs
for f in tempbackup/*.map;do echo $(grep ^0 $f | wc -l) $(basename $f);done > 1.temp
for f in *.map;do echo $(grep ^0 $f | wc -l) $f;done > 2.temp
paste 1.temp 2.temp > final_unmappedSNPs.txt
rm 1.temp 2.temp


## Remove unmapped & convert to binary files
for f in *.map;do
cat $f | awk -F"\t" '{if($1==0)print $2}' > $f.exclude.lst
x=$(cat $f.exclude.lst | wc -l)
if [ $x -gt 0 ];then
  echo $f
  t=${f%.map}.ped
  unmapf=${f%.map}.unmap.map
  unmapt=${f%.map}.unmap.ped
  mv $f $unmapf; mv $t $unmapt;
  plink --file "${f%.map}.unmap" --chr-set 31 no-xy --allow-extra-chr --make-bed --exclude $f.exclude.lst --out "${f%.map}"
  #plink --file "$suffix" --chr-set 31 no-xy --allow-extra-chr --make-bed --exclude $f.exclude.lst --out "$suffix.fixed"
else
  plink --file "${f%.map}" --chr-set 31 no-xy --allow-extra-chr --make-bed --out "${f%.map}"
fi;done &> unmap.txt
rm *.exclude.lst
rm *.{hh,log,nosex}

##CHECK for the No of chromosomes per bim file
for f in *.bim;do echo $f; cat $f  | grep -v "^#" | awk '{print $1}' | sort | uniq -c | sort -k2,2n;done > checkBimChrNo.txt
