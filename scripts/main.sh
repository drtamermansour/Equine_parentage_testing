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
#6
mv 80K_Equcab3_App_Knab_Bellone_Breed_Key.xlsx 80K_Equcab3_App_Knab_Bellone.xlsx

#7
unzip Kurdish_Arabian_670.zip
mv Plink/Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.{map,ped} .
rmdir Plink

#8
#mv New2M_UMN.csv MNEc2M_EquCab2.csv

#change xlsx to csv files (keep original copies in backup_original/overwritten_xlsx)
#670K_Equcab3_Belgian_Haflinger_Bellone.xlsx  
#670K_Equcab3_Mixedbreed_Bellone.xlsx  
#80K_Equcab3_App_Knab_Bellone.xlsx
#New2M_UMN.xlsx

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

# 6. 70K data: Appaloosa (Bellone Lab Rockwell et al. 2020 https://onlinelibrary.wiley.com/doi/epdf/10.1111/age.12883)
# download locally, upload (ERU_GWAS_70K.ped/map) to the drive, then rclone

# 7. 70K data: Icelandic Horses (Bellone Lab 80 K data )
# download locally, upload (ERU_Icelandics_07292021.ped/map) to the drive, then rclone

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
## and move newly derived files to new backup folder
mkdir -p ../backup_drived

#1
#zcat 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.vcf.gz --double-id --allow-extra-chr --recode --out "670k_Equcab2_ThoroughbredsJpn_Fawcett2019"
mv 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.{map,ped,log,nosex} ../backup_drived/.

#2
#cat 670K_EquCab2_Standardbreds_Bellone.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_EquCab2_Standardbreds_Bellone.vcf --double-id --chr-set 35 --recode --out "670K_EquCab2_Standardbreds_Bellone"
mv 670K_EquCab2_Standardbreds_Bellone.{map,ped,log,nosex} ../backup_drived/.

#3
#cat 670K_Equcab3_Mixedbreed_Bellone.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab3_Mixedbreed_Bellone.vcf --double-id --chr-set 31 --recode --out "670K_Equcab3_Mixedbreed_Bellone"
mv 670K_Equcab3_Mixedbreed_Bellone.{map,ped,log,nosex} ../backup_drived/.

#4
#cat 670K_Equcab3_Belgian_Haflinger_Bellone.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab3_Belgian_Haflinger_Bellone.vcf --double-id --chr-set 31 --recode --out "670K_Equcab3_Belgian_Haflinger_Bellone"
mv 670K_Equcab3_Belgian_Haflinger_Bellone.{map,ped,log,nosex} ../backup_drived/.

#5
suffix=AHSdataset
plink --bfile "$suffix" --chr-set 34 no-y no-xy no-mt --allow-extra-chr --recode --out "$suffix"
mv AHSdataset.{map,ped,log,nosex} ../backup_drived/.

#6
suffix=IR.paper.heredity
plink --bfile "$suffix" --chr-set 34 no-y no-xy no-mt --allow-extra-chr --recode --out "$suffix"
mv IR.paper.heredity.{map,ped,log,nosex} ../backup_drived/.

#7
#cat 70K_Equcab3_Bardigiano_Ablondi.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 70K_Equcab3_Bardigiano_Ablondi.vcf --double-id --chr-set 72 no-y no-xy no-mt --allow-extra-chr --recode --out "70K_Equcab3_Bardigiano_Ablondi"
mv 70K_Equcab3_Bardigiano_Ablondi.{map,ped,log,nosex} ../backup_drived/.

#8
#cat 670K_Equcab2_AbagaBlack_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_AbagaBlack_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_AbagaBlack_Han"
mv 670K_Equcab2_AbagaBlack_Han.{map,ped,log,nosex} ../backup_drived/.

#9
#cat 670K_Equcab2_BaichaIronHoof_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_BaichaIronHoof_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_BaichaIronHoof_Han"
mv 670K_Equcab2_BaichaIronHoof_Han.{map,ped,log,nosex} ../backup_drived/.

#10
#cat 670K_Equcab2_Sanhe_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_Sanhe_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_Sanhe_Han"
mv 670K_Equcab2_Sanhe_Han.{map,ped,log,nosex} ../backup_drived/.

#11
#cat 670K_Equcab2_Wushen_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_Wushen_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_Wushen_Han"
mv 670K_Equcab2_Wushen_Han.{map,ped,log,nosex} ../backup_drived/.

#12
#cat 670K_Equcab2_Wuzhumuqin_Han.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_Wuzhumuqin_Han.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_Wuzhumuqin_Han"
mv 670K_Equcab2_Wuzhumuqin_Han.{map,ped,log,nosex} ../backup_drived/.

#13
#cat 70K_Equcab3_MongolianRacing_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 70K_Equcab3_MongolianRacing_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "70K_Equcab3_MongolianRacing_Hill"
mv 70K_Equcab3_MongolianRacing_Hill.{map,ped,log,nosex} ../backup_drived/.

#14
#cat 70K_Equcab3_ArabianHorse_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 70K_Equcab3_ArabianHorse_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "70K_Equcab3_ArabianHorse_Hill"
mv 70K_Equcab3_ArabianHorse_Hill.{map,ped,log,nosex} ../backup_drived/.

#15
#cat 670K_Equcab2_IrishDraught_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_IrishDraught_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_IrishDraught_Hill"
mv 670K_Equcab2_IrishDraught_Hill.{map,ped,log,nosex} ../backup_drived/.

#16
#cat 670K_Equcab2_ConnemaraPony_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_ConnemaraPony_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_ConnemaraPony_Hill"
mv 670K_Equcab2_ConnemaraPony_Hill.{map,ped,log,nosex} ../backup_drived/.

#17
#cat 670K_Equcab2_IrishSportHorse_Hill.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_Equcab2_IrishSportHorse_Hill.vcf --double-id --chr-set 31 no-x no-y no-xy no-mt --recode --out "670K_Equcab2_IrishSportHorse_Hill"
mv 670K_Equcab2_IrishSportHorse_Hill.{map,ped,log,nosex} ../backup_drived/.

#-------
#18
#zcat QH_IlluminaSNP70_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf QH_IlluminaSNP70_EquCab2.vcf.gz --double-id --chr-set 31 no-xy no-mt --recode --out "QH_IlluminaSNP70_EquCab2"
mv QH_IlluminaSNP70_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#19
#zcat STDB_IlluminaSNP70_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf STDB_IlluminaSNP70_EquCab2.vcf.gz --double-id --chr-set 31 no-xy no-mt --recode --out "STDB_IlluminaSNP70_EquCab2"
mv STDB_IlluminaSNP70_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#20
#zcat QH_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf QH_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "QH_IlluminaSNP50_EquCab2"
mv QH_IlluminaSNP50_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#21
#zcat STDB_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf STDB_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "STDB_IlluminaSNP50_EquCab2"
mv STDB_IlluminaSNP50_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#22
#zcat TB_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf TB_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "TB_IlluminaSNP50_EquCab2"
mv TB_IlluminaSNP50_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#23
#zcat Belgian_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf Belgian_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "Belgian_IlluminaSNP50_EquCab2"
mv Belgian_IlluminaSNP50_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#24
#zcat Morgan_IlluminaSNP50_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf Morgan_IlluminaSNP50_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --recode --out "Morgan_IlluminaSNP50_EquCab2"
mv Morgan_IlluminaSNP50_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#24
#zcat STDB_AffyMNEc670K_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf STDB_AffyMNEc670K_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --allow-extra-chr --recode --out "STDB_AffyMNEc670K_EquCab2"
mv STDB_AffyMNEc670K_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#25
#zcat MNEc2M_EquCab2.vcf.gz | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf MNEc2M_EquCab2.vcf.gz --double-id --chr-set 31 no-y no-xy no-mt --allow-extra-chr --recode --out "MNEc2M_EquCab2"
mv MNEc2M_EquCab2.{map,ped,log,nosex} ../backup_drived/.

#26
cd Tozaki
wget https://sapac.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/geneseek-ggp/geneseek-ggp-equine-manifest-file-csv.zip
unzip geneseek-ggp-equine-manifest-file-csv.zip
tail -n+9 GGP_Equine.csv | head -n-24 | awk 'BEGIN{FS=",";OFS="\t"}{print $10,$2,"0",$11}' > GGP_Equine.design
for f in *.txt;do out=${f%.txt}; echo $out;
tail -n+10 $f | awk 'BEGIN{FS=OFS="\t"}{if($2)print $2,$2,$1,$5,$6}' > $out.lgen
tail -n+10 $f | awk 'BEGIN{FS=OFS="\t"}{if($2)print $2,$2,"0","0","0","0"}' | sort | uniq > $out.fam
cp GGP_Equine.design $out.map
plink --horse --nonfounders --allow-no-sex --lfile $out --missing-genotype - --output-missing-genotype 0 --recode --out $out.drived
done
mv *.drived.{map,ped,log,nosex} ../../backup_drived/.

#27
#cat 670K_192WP_32WP_Combined.vcf | grep -v "^#" | awk '{print $1}' | uniq
plink --vcf 670K_192WP_32WP_Combined.vcf --double-id --chr-set 31 --allow-extra-chr --recode --out "670K_192WP_32WP_Combined"
mv 670K_192WP_32WP_Combined.{map,ped,log,nosex} ../backup_drived/.

#############
## Move all ped/map files to the work_dir
cd $work_dir/.
cp backup_{original,drived}/*.{ped,map} $work_dir/.

############
#### download genomes
# Y chromosome: According to Terje Raudsepp, The Y chromosome assembly is still the one from 2018 Janecka et al, thus eMSYv3. The single copy assembly region is the one suitable for your SNP panel. Barbara Wallner would be a better contact for useful Y SNPs.
mkdir -p "$work_dir"/eMSY && cd "$work_dir"/eMSY ## download locally from Terje's email, upload to the drive, then rclone
rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs/eMSYv3.fasta $HOME/Horse_parentage_SNPs/eMSY/
# Janecka et al. 2018 (data availability section and Description of Additional Supplementary Files) states that this asembly is eMSYv3.1, in which  two vector sequences (18,381-bp in position 2,505,950 and 1396-bp in position 8,066,811) have been trimmed from assembly with respect to v3.0. Therefore, to obtain v3, you must insert a 1396‐bp gap sequence (Ns) after position 8,048,430 and another 18,381‐bp gap sequence (Ns) after position 2,506,950 into the eMSYv3.1 sequence. The total sequence length after making these two insertions should be 9,497,449 bp.
seq=$(head -n2 $HOME/Horse_parentage_SNPs/eMSY/eMSYv3.fasta | tail -n1)
point1=2506950 ## it was mentioned as 2505950 in the data availability section but both positions should not affect our SNPs
point2=8048430
interval_len=$((point2-point1))
echo ${seq:0:$point1} | tr -d '\n' > seq1.part
echo ${seq:$point1:$interval_len} | tr -d '\n' > seq2.part
echo ${seq:$point2} > seq3.part
for i in {1..18381}; do echo -n "N";done > rep1.part
for i in {1..1396}; do echo -n "N";done > rep2.part
head -n1 eMSYv3.fasta > eMSYv3_nontrimmed.fasta
cat seq1.part rep1.part seq2.part rep2.part seq3.part >> eMSYv3_nontrimmed.fasta

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
head -n1 "$work_dir"/eMSY/eMSYv3_nontrimmed.fasta >> equCab3/equCab3_genome.fa
head -n2 "$work_dir"/eMSY/eMSYv3_nontrimmed.fasta | tail -n1 | fold -w50 >> equCab3/equCab3_genome.fa
equCab3_ref="$work_dir"/equCab3/equCab3_genome.fa
cat $equCab3_ref | grep -v "^#" | awk -F">" '{if(!$1){if(NR!=1)print "##contig=<ID="ref",length="reflen">";ref=$2;reflen=0}else reflen+=length($1)}END{print "##contig=<ID="ref",length="reflen">";}' > equCab3/vcf_contigs.txt
#head -n2 "$work_dir"/eMSY/eMSYv3_nontrimmed.fasta | awk -F">" '{if(!$1){if(NR!=1)print "##contig=<ID="ref",length="reflen">";ref=$2;reflen=0}else reflen+=length($1)}END{print "##contig=<ID="ref",length="reflen">";}' >> equCab3/vcf_contigs.txt
equCab3_vcfContigs="$work_dir"/equCab3/vcf_contigs.txt

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/863/925/GCF_002863925.1_EquCab3.0/GCF_002863925.1_EquCab3.0_assembly_report.txt
grep -v "^#" GCF_002863925.1_EquCab3.0_assembly_report.txt | sed -e "s/\r//g" | awk 'BEGIN{FS=OFS="\t"}{print
$10,$5}' | sed 's/^chr//' > equCab3_chrMap
equCab3_chrMap="$work_dir/equCab3/equCab3_chrMap"


cat $equCab3_ref | awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' > "$work_dir"/equCab3/equCab3_genome_unwrap.fa
equCab3_ref_unwrap="$work_dir"/equCab3/equCab3_genome_unwrap.fa

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
#echo -e "[A/C]\n[A/G]\n[T/C]\n[T/G]" | grep -vFwf - <(tail -n+2 GGP_Equine.simple.csv) ## all SNPs are unambiguous and the A or T alleles always comes first
#cat GGP_Equine.simple.csv | cut -f3 | sort | uniq -c  ## 31 X(3409) Y(2)
## QC
cat GGP_Equine.simple.csv | awk -F"\t" '{print $3"_"$4}' | sort | uniq -c | awk '{if($1>1)print $2}' > dup_loci ## 55 SNPs (different IDs but same position)
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$3"_"$4]==1){a[$3"_"$4]++;b[$3"_"$4]=$1;}else if(a[$3"_"$4]>1)b[$3"_"$4]=b[$3"_"$4] FS $1;}END{for(i in b)print i,b[i]}' dup_loci GGP_Equine.simple.csv > dup_snps


# tranform the array manifest into VCF & generate ALT/REF map
bash "$scripts"/manifest_to_ref.sh "GGP_Equine.simple.csv" "$equCab2_vcfContigs" "$equCab2_ref" "GGP_Equine" ## This script excludes the Y SNPs
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

## merge the annotation to the map file which has the reference alleles but not for the SNPs that failed to remap
cat SNP70.unique_remap.FINAL.csv | sed 's/chrUn_ref|/Un_/' | sed 's/\.1|,/v1,/' | sed 's/,chr/,/g' | tr ',' '\t' > SNP70_remap.tab
echo -e "$(head -n1 ann_70.tab)\t$(head -n1 SNP70_remap.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$6,$8,$9,$10}')" > map_70.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$5"_"$7]=$1 FS $3 FS $6 FS $8 FS $9 FS $10;next}{if(a[$3"_"$4])print $0,a[$3"_"$4];else print $0,0,0,0,0,0,0;}' SNP70_remap.tab <(tail -n+2 ann_70.tab) >> map_70.tab

## (Optional) add reference alleles for all SNPs in Equcab2 (except Un2 chromosome bc we do not have its sequence)
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
#echo -e "[A/C]\n[A/G]\n[T/C]\n[T/G]" | grep -vFwf - <(tail -n+2 Axiom_MNEc670_Annotation.r2.simple.csv | awk -F"\t" '{print "["$7"/"$8"]"}') | wc -l ## we have 8790 ambiguous SNPs (one of them is [CT/GA])
#echo -e "[A/T]\n[C/G]\n[T/A]\n[G/C]\n[CT/GA]" | grep -Ff - Axiom_MNEc670_Annotation.r2.simple.csv | grep -w "-" | wc -l ## All of ambiguous SNPs are on +ve strand
#cat Axiom_MNEc670_Annotation.r2.simple.csv | grep -w "-" | wc -l ## actually all the array SNPs are on +ve strand but my QC (during tranformation of the array manifest into VCF) shows 22 SNPs on Chr14 that I found coming from -ve strand)
echo -e "[A/T]\n[C/G]\n[T/A]\n[G/C]\n[CT/GA]" | grep -Ff - Axiom_MNEc670_Annotation.r2.simple.csv > ambiguous_670 ## 8790
ambiguous_670="$work_dir/SNParrays/670k/ambiguous_670"


#cat Axiom_MNEc670_Annotation.r2.simple.csv | cut -f3 | sort | uniq -c  ## 31 Un(21468) Un2(9395) X(28017) Y(1)
## QC
cat Axiom_MNEc670_Annotation.r2.simple.csv | awk -F"\t" '{print $3"_"$4}' | sort | uniq -c | awk '{if($1>1)print $2}' > dup_loci  ## 0 SNPs (different IDs but same position)


# tranform the array manifest into VCF & generate ALT/REF map (except Un2 chromosome bc we do not have its sequence)
grep -vw "Un2" Axiom_MNEc670_Annotation.r2.simple.csv > Axiom_MNEc670_Annotation.r2.simple_excludeUN2.csv
bash "$scripts"/manifest_to_ref.sh "Axiom_MNEc670_Annotation.r2.simple_excludeUN2.csv" "$equCab2_vcfContigs" "$equCab2_ref" "Axiom_MNEc670_r2" ## This script excludes the Y SNPs
arr670k="$work_dir/SNParrays/670k/Axiom_MNEc670_r2.alt_ref_comp.tab"

#gf="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
#cat Axiom_MNEc670_r2.vcf | awk -v gf=$gf 'BEGIN{FS=OFS="\t"}/##/{print;next}/#/{print gf; print $0,"FORMAT","simFake";next}{print $0,"GT","0/1"}' > Axiom_MNEc670_r2.simFake.vcf
#plink --vcf Axiom_MNEc670_r2.simFake.vcf --chr-set 31 no-xy --allow-extra-chr --out Axiom_MNEc670_r2.simFake
#plinkref670=$(pwd)/Axiom_MNEc670_r2.simFake

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

## merge the annotation to the map file which has the reference alleles but not for the SNPs that failed to remap
cat MNEc670k.unique_remap.FINAL.csv | awk -F"," '{if($10)print}'  | sed 's/chrUn_ref|/Un_/' | sed 's/\.1|,/v1,/' | sed 's/,chr/,/g' | tr ',' '\t' > MNEc670k_remap.tab
echo -e "$(head -n1 ann_670.tab)\t$(head -n1 MNEc670k_remap.tab | awk 'BEGIN{FS=OFS="\t"}{print $1,$3,$6,$8,$9,$10}')" > map_670.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$5"_"$7]=$1 FS $3 FS $6 FS $8 FS $9 FS $10;next}{if(a[$3"_"$4])print $0,a[$3"_"$4];else print $0,0,0,0,0,0,0;}' MNEc670k_remap.tab <(tail -n+2 ann_670.tab) >> map_670.tab

## (Optional) add reference alleles for all SNPs in Equcab2 (except Un2 chromosome bc we do not have its sequence)
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

paste <(tail -n+2 $map670 | awk 'BEGIN{FS=OFS="\t"}{print $1,$11,$10}') \
<(tail -n+2 $map670 | awk 'BEGIN{FS=OFS="\t"}{print $11,$10}' | tr 'TCGA' 'AGCT') > Axiom_MNEc670_r2.alt_ref_comp.EC3.tab
arr670k_EC3="$work_dir/SNParrays/670k/Axiom_MNEc670_r2.alt_ref_comp.EC3.tab"

## compare to the Equcab3 manifest file
wget http://media.affymetrix.com/analysis/downloads/lf/genotyping/Axiom_MNEc670/r3/Axiom_MNEc670.na35.r3.a2.annot.csv.zip
unzip Axiom_MNEc670.na35.r3.a2.annot.csv.zip
sed -i 's/","/"|"/g' Axiom_MNEc670.na35.r3.a2.annot.csv
cat Axiom_MNEc670.na35.r3.a2.annot.csv | grep -v "^#" | awk 'BEGIN{FS="|";OFS="\t"}{gsub("^\"NW_", "\"Un_NW_", $5);gsub("\\.1\"", "v1\"", $5);print $1,$2,$5,$6,$8,$11,$12,$13,$14,$15,$33,$34,$36}' | sed 's/"//g' > Axiom_MNEc670.na35.r3.a2.annot.simple.csv

#cat Axiom_MNEc670.na35.r3.a2.annot.simple.csv | cut -f3 | sort | uniq -c  ## 0(39127) 31 X(27327) Y(1) + many unplaced contigs
echo -e "$(head -n1 $map670)\tarrEC3_chrom\tarrEC3_pos\tarrEC3_REF\tarrEC3_ALT" > map_670_complete2.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$3 FS $4 FS $9 FS $10;next}{if(a[$1])print $0,a[$1];else print $0,0,0,0,0;}' "Axiom_MNEc670.na35.r3.a2.annot.simple.csv" <(tail -n+2 $map670) | sed 's/---/0/g' >> map_670_complete2.tab

tail -n+2 map_670_complete2.tab | awk -F"\t" '{if($8==$12 && $8==0)print}' > unmapped_670 ## 26962
tail -n+2 map_670_complete2.tab | awk -F"\t" '{if($8!=$12 || $9!=$13 || $10!=$14 || $11!=$15)print}' | wc -l ## 14184 ## inconsistant remapping bettween illumina & Molly
tail -n+2 map_670_complete2.tab | awk -F"\t" '{if($8!=$12 && $12==0)print}' > remap_IlluminaFail_670 ## 12165
tail -n+2 map_670_complete2.tab | awk -F"\t" '{if($8!=$12 && $8==0)print}' > remap_MollyFail_670 ## 1920
tail -n+2 map_670_complete2.tab | awk -F"\t" '{if($8!=$12 && $8!=0 && $12!=0)print}' > remap_IlluminaMollyIncons1_670 ## 19
tail -n+2 map_670_complete2.tab | awk -F"\t" '{if($9!=$13 && $8==$12 && $8!=0 && $12!=0)print}' > remap_IlluminaMollyIncons2_670 ## 73
tail -n+2 map_670_complete2.tab | awk -F"\t" '{if(($10!=$14 || $11!=$15) && $9==$13 && $8==$12 && $8!=0 && $12!=0)print}' > remap_IlluminaMollyIncons3_670 ## 7 ## Illumina reports the position correct but report the alleleles as 0

## 80k
## There is no official manifest for the 80k array
## instead, I will use SNP_ids from all the map files I have for this array, then extract those SNPs from 70k & 670k maps
## Please note that this code assumes that the map files have been downloaded & passed the first 3 fixing steps
## Note: We asked Jiansheng Qiu for the 80k manifest. He sent us the map file "GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map"
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


## Remapping files from University of Minnesota Equine Genetics and Genomics Lab
## https://github.com/UMN-EGGL/MNEc2to3
mkdir $work_dir/SNParrays/MNEc2to3 && cd $work_dir/SNParrays/MNEc2to3
wget https://raw.githubusercontent.com/UMN-EGGL/MNEc2to3/master/data/MNEc2M.probe_blast_counts.csv.gz
tail -n+2 MNEc2M.probe_blast_counts.csv | awk -F"," '{if($11=="True")print $5}' | sort | uniq -c
# Note: There are multiple records with the same coordinates but different SNP ID e.g.
# 24940,MNEc.2.1.29763.UKUL2,Affx-101643281,AX-104566380,chr1,29763,T,C,C,T,True,True,chr1.29763,1,1,1
# 24941,MNEc.2.1.29763.UKUL2,Affx-101643281,AX-103622354,chr1,29763,T,C,C,T,True,True,chr1.29763,1,1,1

cd $work_dir/SNParrays
## assess overlap in Equcap2
tail -n+2 50k/SNP50.NCBI.remap.csv | awk -F"," '{print $4"."$6}' | sed 's/chr//' | sed 's/|/_/' | sort | uniq > 50.pos
tail -n+2 $map70 | awk -F"\t" '{print $4"."$5}' | sed 's/chr//' | sed 's/v1././' | sort | uniq > 70.pos
tail -n+2 $map670 | awk -F"\t" '{print $4"."$5}' | sed 's/chr//' | sed 's/v1././' | sort | uniq > 670.pos
comm -12 670.pos 50.pos > 670-50.pos
comm -12 670.pos 70.pos > 670-70.pos
comm -12 70.pos 50.pos > 70-50.pos
comm -12 670-50.pos 70.pos > 670-50-70.pos
wc -l *.pos
#   54602 50.pos
#   44098 670-50-70.pos
#   50936 670-50.pos
#   62277 670-70.pos
#  670796 670.pos
#   45986 70-50.pos
#   65102 70.pos

## assess overlap in Equcap3
tail -n+2 50k/SNP50.NCBI.remap.csv | awk -F"," '{if($5)print $5"-"$7}' | cut -d"." -f1 | sed 's/chr//' | sed 's/-/./' | sed 's/|/_/' | sort | uniq > 50.pos2
tail -n+2 $map70 | awk -F"\t" '{if($8)print $8"."$9}' | sed 's/chr//' | sed 's/v1././' | sort | uniq > 70.pos2
tail -n+2 $map670 | awk -F"\t" '{if($8)print $8"."$9}' | sed 's/chr//' | sed 's/v1././' | sort | uniq > 670.pos2
comm -12 670.pos2 50.pos2 > 670-50.pos2
comm -12 670.pos2 70.pos2 > 670-70.pos2
comm -12 70.pos2 50.pos2 > 70-50.pos2
comm -12 670-50.pos2 70.pos2 > 670-50-70.pos2
wc -l *.pos2
#   55571 50.pos2
#   42223 670-50-70.pos2
#   49582 670-50.pos2
#   59050 670-70.pos2
#  641918 670.pos2
#   44174 70-50.pos2
#   61939 70.pos2

## map between 70k and 670k (Equcap3)
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1"_"$2]=$3 FS $4 FS $5;next}{if(a[$1"_"$2])print $1,$2,a[$1"_"$2],$3,$4,$5}'  <(grep -v "^#" 70k/GGP_Equine.vcf) <(grep -v "^#" 670k/Axiom_MNEc670_r2.vcf) > 70k_670k_map.txt ## 62276
#cat 70k_670k_map.txt | awk -F"\t" '{if($4!=$7)print}'
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8"_"$9]=$1 FS $10 FS $11;next}{if(a[$8"_"$9])print $8"_"$9,a[$8"_"$9],$1,$10,$11}'  $map70 $map670 > 70k_670k_map.txt ## 59051

## inconsistent Molly remapping between 70k and 670k
awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$5"."$7]=$6"."$8;next}{if(a[$5"."$7] && a[$5"."$7] != $6"."$8)print $0,a[$5"."$7]}' 70k/SNP70.unique_remap.FINAL.csv 670k/MNEc670k.unique_remap.FINAL.csv > molly.remap.inconsistent ## Same Equicab2 position remap to different Equicab3 positions
#awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$5"."$7]=$6"."$8;next}{if(!a[$5"."$7])print $0}' 670k/MNEc670k.unique_remap.FINAL.csv 70k/SNP70.unique_remap.FINAL.csv > molly.remap.missed
awk 'BEGIN{FS=",";OFS="\t"}FNR==NR{a[$5"."$7]=1;b[$6"."$8]=$0;next}{if(!a[$5"."$7] && b[$6"."$8])print b[$6"."$8],$5"."$7}' 670k/MNEc670k.unique_remap.FINAL.csv 70k/SNP70.unique_remap.FINAL.csv > molly.remap.inconsistent2 ## Two different Equicab2 positions to the remap to the same Equicab3 position

## a QC to differentiate 80k from 70k arrays
cd ~/Horse_parentage_SNPs/SNParrays/70k
cut -f3,4 GGP_Equine.simple.csv | tr '\t' '.' | sort | uniq > temp.eq2.70Pos
cut -f3,4,11,12 ../670k/map_670.tab | tr '\t' '.' | sort | uniq > temp.eq2_eq3.670Pos
grep -v -Ff temp.eq2.70Pos temp.eq2_eq3.670Pos | cut -d"." -f3,4 | grep -v ^0 | sort | uniq > temp.eq3.not70Pos
not70="$HOME/Horse_parentage_SNPs/SNParrays/70k/temp.eq3.not70Pos"

cd $work_dir
for f in *.bim;do echo $f $(cat $f | wc -l) $(grep -Fwf $not70 <(cut -f1,4 $f | tr '\t' '.' | sort) | wc -l);done > temp.not70
####################
## CHECK: check map files format
for f in *.map;do echo $(wc -l $f) "--" $(cat -v $f | grep -v "^#" | head -n1) "--" $(head -n1 $f | awk -F"\t" '{print NF}');done > checkMapFormat.txt

## CHECK: check ped files format
for f in *.ped;do echo $f $(grep -v "^#" $f | head -n1 | cut -f1-10);done > checkPedFormat.txt
####################
### Fix input map/ped files format
## Fix the files with abonormal end characters
for f in *.{map,ped};do
 sed -i -e "s/\r//g" $f;
done

####################
## Map file check results:

## FIX: Unify the file separator
for f in *.map;do if [ -f $f ];then cp $f $f.temp; cat $f.temp | grep -v "^#" | tr ' ' '\t' > $f;rm $f.temp;fi;done

## 80K_Equcab3_MultipleBreeds_vit.map has 3 columns only
## Fix the missing columns in the map files
cp 80K_Equcab3_MultipleBreeds_vit.map 80K_Equcab3_MultipleBreeds_vit.map.temp
cat 80K_Equcab3_MultipleBreeds_vit.map.temp | awk 'BEGIN{FS=OFS="\t"}{$3=0 FS $3;print}' > 80K_Equcab3_MultipleBreeds_vit.map
rm 80K_Equcab3_MultipleBreeds_vit.map.temp

####################
## Ped file check results:
## Some ped files are Irregularly-formatted PLINK text files
##    - 670K_Equcab2_Exmoor_Lindgren-Velie.ped
##    - 670K_Equcab2_Icelandic_Lindgren-Velie.ped
## Some ped files are numberically encoded
##    - 670K_Equcab3_FranchesMontagnes_Gmel.ped
## Some ped files are using A/B call
##    - Clyde_70K.ped


## Fix the Irregularly-formatted PLINK text files (part1)
cp Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.ped Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.ped.toBedeleted
cat Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.ped.toBedeleted | grep -v "^#" | awk 'BEGIN{FS=OFS="\t"}{gsub(" ","_",$1);$1=1 FS $1 FS 0 FS 0 FS 0 FS "-9";print $0}' > Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.ped
rm Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.ped.toBedeleted

## FIX: unify the file separator
for f in *.ped;do if [ -f $f ];then cp $f $f.temp; cat $f.temp | grep -v "^#" | tr ' ' '\t' > $f;rm $f.temp;fi;done

## Fix the Irregularly-formatted PLINK text files (part2)
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
#rm Clyde_70K.{map,ped} ## It was updated by new dataset
#rm ERU_GWAS_70K.{map,ped} ## It was replaced by ERU_65K_09102016
#################
##CHECK for the No of chromosomes per map file
#cat 80K_Equcab3_App_Knab_Bellone.map  | grep -v "^#" | awk '{print $1}' | sort | uniq -c | sort -k2,2n
for f in *.map;do echo $f; cat $f  | grep -v "^#" | awk '{print $1}' | sort | uniq -c | sort -k2,2n;done > checkMapChrNo.txt

## Note2: chromosome IDs
# 80K_Equcab3_App_Knab_Bellone.map has PJAA* chromosomes. These need to be mapped to the UCSC format
# 670k_Equcab2_ThoroughbredsJpn_Fawcett2019.map has CM* chr. These need to be mapped to the UCSC format
# ERU_Icelandics_07292021.map has many unknown contigs
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

cp 670K_192WP_32WP_Combined.map 670K_192WP_32WP_Combined.map.temp
cat 670K_192WP_32WP_Combined.map.temp | sed 's/^chrUn/Un/; s/^chrUn2/Un2/; s/^32/X/; s/^33/Y/;' > 670K_192WP_32WP_Combined.map
rm 670K_192WP_32WP_Combined.map.temp
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
"670K_Equcab3_Mixedbreed_Bellone.map" \
"Clyde_70K.map" \
"JapaneseThoroughbred.drived.map" \
"Misaki.drived.map" \
"Noma.drived.map" \
"Hokkaido.drived.map" \
"Kiso.drived.map" \
"Miyako.drived.map" \
"Taishu.drived.map" \
"ERU_65K_09102016.map";do
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
"670K_EquCab2_Standardbreds_Bellone.map" \
"Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.map";do
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
"STDB_IlluminaSNP70_EquCab2.map" \
"JapaneseThoroughbred.drived.map" \
"Misaki.drived.map" \
"Noma.drived.map" \
"Hokkaido.drived.map" \
"Kiso.drived.map" \
"Miyako.drived.map" \
"Taishu.drived.map" \
"ERU_65K_09102016.map";do
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
"MNEc2M_EquCab2.map" \
"670K_192WP_32WP_Combined.map" \
"Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink.map";do
echo $input "--" $(comm -23 <(awk -F"\t" '{print $2}' "$input" | sort) <(tail -n+2 "$map670" | awk -F"\t" '{print $1}' | sort) | wc -l);
done >> unrecognized_SNP_ids.txt

#arr="80"
for input in \
"EQ75v03_EquCab2_Falabella_DiazS.map" \
"70K_Equcab2_Maremmano_Capomaccio.map" \
"GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map" \
"80K_Equcab3_App_Knab_Bellone.map" \
"80K_Equcab3_MultipleBreeds_vit.map" \
"Clyde_70K.map";do
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
#ERU_Icelandics_07292021.map -- 6359
# postpone fixing unrecognized SNP ids after removing duplicates


##CHECK 4. check for duplicate SNP positions
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
sleep 3
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

############################
mkdir -p tempbackup2
cp *.fam tempbackup2/.

### Update sample IDs to reflect the study and breed
cd $work_dir
cp backup_original/{670K_Equcab3_Belgian_Haflinger_Bellone.csv,\
670K_Equcab3_Mixedbreed_Bellone.csv,\
80K_Equcab3_App_Knab_Bellone.csv,\
80K_Equcab3_MultipleBreeds_vit.csv,\
New2M_UMN.csv} .

## fix the end of line characters
for f in 670K_Equcab3_Belgian_Haflinger_Bellone.csv \
670K_Equcab3_Mixedbreed_Bellone.csv \
80K_Equcab3_App_Knab_Bellone.csv \
80K_Equcab3_MultipleBreeds_vit.csv \
New2M_UMN.csv;do
 sed -i -e "s/\r//g" $f;
done

## remove header and reorder the columns 
cp 80K_Equcab3_MultipleBreeds_vit.csv 80K_Equcab3_MultipleBreeds_vit.csv_temp
tail -n+2 80K_Equcab3_MultipleBreeds_vit.csv_temp | awk 'BEGIN{FS=OFS=","}{print $2,$1,$3}' > 80K_Equcab3_MultipleBreeds_vit.csv
rm 80K_Equcab3_MultipleBreeds_vit.csv_temp

cp New2M_UMN.csv New2M_UMN.csv_temp
tail -n+2 New2M_UMN.csv_temp | awk 'BEGIN{FS=OFS=","}{print $2,$1}' > New2M_UMN.csv
rm New2M_UMN.csv_temp

## edit the breed column to remove spaces
for f in 670K_Equcab3_Belgian_Haflinger_Bellone.csv \
670K_Equcab3_Mixedbreed_Bellone.csv \
80K_Equcab3_App_Knab_Bellone.csv \
New2M_UMN.csv \
80K_Equcab3_MultipleBreeds_vit.csv;do
cp $f $f.temp;
awk 'BEGIN{FS=OFS=","} {gsub(" ", "", $2)} 1' $f.temp > $f;
rm $f.temp
done

## generate newIds
f=670K_Equcab3_Belgian_Haflinger_Bellone
cut -d"," -f2 $f.csv | sort | uniq -c
cat $f.csv | sed 's/,Belgian.*$/ BE/; s/,Haflinger.*$/ HF/' > $f.temp
awk 'BEGIN{FS=" ";OFS="\t"}FNR==NR{a[$2]=$1 OFS $2;next}{print a[$1],"670k_Bel2",$2"_"FNR;}' $f.fam $f.temp > $f.newIds
rm $f.temp

f=670K_Equcab3_Mixedbreed_Bellone
cut -d"," -f2 $f.csv | sort | uniq -c
cat $f.csv | sed 's/,Appaloosa.*$/ AP/; s/,Friesian.*$/ FR/; s/,Knabstrupper.*$/ KNB/; s/,Knapstrupper.*$/ KNB/; s/,PonyoftheAmericas.*$/ POA/; s/,ShetlandPony.*$/ SP/;' > $f.temp
awk 'BEGIN{FS=" ";OFS="\t"}FNR==NR{a[$2]=$1 OFS $2;next}{print a[$1],"670k_Bel3",$2"_"FNR;}' $f.fam $f.temp > $f.newIds
rm $f.temp

f=80K_Equcab3_App_Knab_Bellone
cut -d"," -f2 $f.csv | sort | uniq -c
cat $f.csv | sed 's/,Appaloosa.*$/ AP/; s/,Knabstrupper.*$/ KNB/; s/,PonyoftheAmericas.*$/ POA/; s/,QuarterHorse.*$/ QH/;' > $f.temp
awk 'BEGIN{FS=" ";OFS="\t"}FNR==NR{a[$2]=$1 OFS $2;next}{print a[$1],"80k_Bel6",$2"_"FNR;}' $f.fam $f.temp > $f.newIds
rm $f.temp


f=80K_Equcab3_MultipleBreeds_vit
cut -d"," -f2 $f.csv | sort | uniq -c
cat $f.csv | sed 's/,Haflingerhorse.*$/ HF/; s/,HolsteinWarmblood(HOL).*$/ HO/; s/,OldenburgJumpingHorse(OS).*$/ OS/; s/,OldenburgWarmblood(OL).*$/ OL/; s/,SwissWarmblood.*$/ SWB/; s/,Trakehner(TRAK).*$/ TK/; s/,WestfalianWarmblood(WESTF).*$/ WF/;' > $f.temp
awk 'BEGIN{FS=" ";OFS="\t"}FNR==NR{a[$2]=$1 OFS $2;next}{print a[$1],"80k_Vit1",$2"_"FNR;}' $f.fam $f.temp > $f.newIds
rm $f.temp

cut -d" " -f2 MNEc2M_EquCab2.fam | grep -Fwf - New2M_UMN.csv > MNEc2M_EquCab2.csv
f=MNEc2M_EquCab2
cut -d"," -f2 $f.csv | sort | uniq -c
cat $f.csv | sed 's/,Arabian.*$/ AR/; s/,Belgian.*$/ BE/; s/,FrenchTrotter.*$/ FT/; s/,Hanoverian.*$/ HN/; s/,Icelandic.*$/ ICH/; s/,Lusitano.*$/ LU/; s/,Maremmano.*$/ MRM/; s/,Morgan.*$/ MH/; s/,QuarterHorse.*$/ QH/; s/,Standardbred.*$/ STB/; s/,Thoroughbred.*$/ TB/; s/,Unknownbreed.*$/ UN/; s/,WelshPony.*$/ WP/;' > $f.temp
awk 'BEGIN{FS=" ";OFS="\t"}FNR==NR{a[$2]=$1 OFS $2;next}{print a[$1],"MNEc2M_McCue2",$2"_"FNR;}' $f.fam $f.temp > $f.newIds
rm $f.temp

f=Univ_of_Kentucky_Bailey_EQNAFFY_20170105_Plink
cut -d" " -f1-2 $f.fam | awk 'BEGIN{OFS="\t"}{$3="670k_Pet2";$4=$2;gsub("^A.*","AR.Per",$4);gsub("^K.*","KD",$4);print $0"_"NR}' > $f.newIds

cp backup_original/id_map.csv .
while IFS=',' read st st_id br;do if [ ! -f "$st.newIds" ];then
echo $st $st_id $br
cat $st.fam | awk -v st_id=$st_id -v br=$br 'BEGIN{FS=" ";OFS="\t"}{print $1,$2,st_id,br"_"NR}' > $st.newIds
fi;done < id_map.csv #<(head -n1 id_map.csv)

> update-ids.log
for nst in *.newIds;do if [ ! -f "$nst.fam" ];then
sleep 3
st=${nst%.newIds}
plink --bfile "$st" --chr-set 31 no-xy --allow-extra-chr --update-ids $st.newIds --make-bed --out "$st.newIds"
fi;done &> update-ids.log

mkdir -p newIds
mv *.newIds* newIds/.
############################
## Fix the alleles to match the +ve strand of Equcab3
cd newIds
# fix the bim file where both alleles are abscent
for suffix in *.newIds;do
  cp "$suffix".bim "$suffix".bim_temp
  awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $3;next;}{if($5 == "0" && $6 == "0" && a[$2])print $1,$2,$3,$4,a[$2];else print $0}' $arr670k_EC3 "$suffix".bim_temp > "$suffix".bim
  rm "$suffix".bim_temp
done

# Prepare for the comming steps
for suffix in *.newIds;do
  awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{if(a[$2])print $0,a[$2]}' $arr670k_EC3  "$suffix".bim > "$suffix".bim.ext
done

# generate lists of bad SNPs to be excluded
for suffix in *.newIds;do
  cat "$suffix".bim.ext | awk 'BEGIN{FS=OFS="\t"}{
    if($5 != 0 && $6 != 0){
    if( ! (($5 == $8 && $6 == $9) || ($5 == $9 && $6 == $8) || ($5 == $10 && $6 == $11) || ($5 == $11 && $6 == $10))) print
    }}' > "$suffix".badSNPs
done

# exclude bad SNPs
for f in *.badSNPs;do
  if [ -s $f ];then echo $f;
  else rm $f;fi
done

mkdir -p failingStudies
mv 670K_Equcab3_FranchesMontagnes_Gmel.* failingStudies/.
mv 70K_Equcab3_Bardigiano_Ablondi.* failingStudies/.

suffix="EQN65KV02_EquCab2_ArgentineanCreole_GiovambattistaG-DiazS.newIds"
mv $suffix.bed $suffix.temp.bed
mv $suffix.bim $suffix.temp.bim
mv $suffix.fam $suffix.temp.fam
mv $suffix.log $suffix.temp.log
plink --bfile "$suffix".temp --chr-set 31 no-xy --allow-extra-chr --exclude  <(cut -f2 $suffix.badSNPs) --make-bed --out "$suffix"
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$0;next;}{if(a[$2])print $0,a[$2]}' $arr670k_EC3  "$suffix".bim > "$suffix".bim.ext

# generate the SNP list to be updated
for suffix in *.newIds;do
  cat "$suffix".bim.ext | awk 'BEGIN{FS="\t"}{
  if($5 != $8 && $5 != $9){
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
  }else if($6 == "0"){
   if($5 == $8)print $2,$5,$6,$8,$9;
   else if($5 == $9)print $2,$5,$6,$9,$8;
  }
  }' > "$suffix".SNP_list_update.txt
done
mkdir -p bim_ext
mv *.bim.ext bim_ext/.


# update the bim/bed files (complete the missing allele and fix strand issues)
for suffix in *.newIds;do
  if [ -s "$suffix".SNP_list_update.txt ];then
    sleep 3
    plink --bfile "$suffix" --chr-set 31 no-xy --allow-extra-chr --update-alleles "$suffix".SNP_list_update.txt --make-bed --out "$suffix".update
  else
    cp "$suffix".bed "$suffix".update.bed
    cp "$suffix".bim "$suffix".update.bim
    cp "$suffix".fam "$suffix".update.fam
fi;done &> strand_update.log
mkdir -p SNP_lists
mv *.SNP_list_update.txt SNP_lists/.

# Fix the REF/ALT alleles if needed
for suffix in *.newIds;do
  sleep 3
  plink --bfile "$suffix".update --chr-set 31 no-xy --allow-extra-chr --make-bed --a2-allele $arr670k_EC3 3 1 --out "$suffix".swab
done &> refAlt_update.log
rm *.update*
rm *.hh
#####################################
## merge all studies
for suffix in *.newIds;do
echo "$suffix".swab;done > studies_to_merge
plink --merge-list studies_to_merge --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --make-bed --out allStudies
#Total genotyping rate is 0.38255.
#641914 variants and 8855 samples pass filters and QC.

## exclude ambiguous SNPs
plink --bfile allStudies --chr-set 31 no-xy --allow-extra-chr --exclude  <(cut -f1 $ambiguous_670) --make-bed --out allStudies.nonAmb  ## There are 8790 ambiguous SNPs but only 7372 in the merged file
#Total genotyping rate is 0.383443.
#634542 variants and 8855 samples pass filters and QC.

mv allStudies.nonAmb.fam allStudies.nonAmb.fam_temp
cat allStudies.nonAmb.fam_temp | awk '{print $1,$2,"0","0","0","-9"}' > allStudies.nonAmb.fam
#########################################
## I am using IBS distance for pruning of sample duplicates
## Plink has --rel-cutoff to exclude one member of each pair of samples with observed genomic relatedness greater than the given cutoff value (default 0.025). However, I am not using this approach to control which sample to be excluded
## Instead, I am generating a distance matrix using --distance but we can not run this on all variants
## I tried using variants of 1 chromosome, pruning, filtering with MAF and geno. I found that filtration on low --geno cutoff is the easiest and the best which also make sens because it harmonize different arrays
plink --bfile allStudies.nonAmb --distance square allele-ct ibs 1-ibs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --geno 0.01 --out allStudies.nonAmb.geno.distance #15753 variants and 8855 samples pass filters and QC.
# tranform the distance matrix into one colum
cat allStudies.nonAmb.geno.distance.mibs | awk 'BEGIN{FS="\t";OFS="\n";}{for (i = 1; i <= NF; i++)print $i}' > temp.gen.1
# Prep the ids
awk 'BEGIN{OFS="\t"}FNR==NR{a[FNR]=$1;next}{for (i in a)print $1,a[i]}'  <(cat allStudies.nonAmb.geno.distance.mdist.id | tr '\t' '.') <(cat allStudies.nonAmb.geno.distance.mdist.id | tr '\t' '.') > temp.gen.2
# create a pairwise table for significantly similar samples
paste temp.gen.2 temp.gen.1 | awk 'BEGIN{FS=OFS="\t"}{if($3>0.95 && $1!=$2)print}' > ibs.gen.95_temp
cat ibs.gen.95_temp | awk 'BEGIN{FS=OFS="\t"}{if(!a[$2"."$1] && $3!="nan"){a[$1"."$2]=1;print}}' > ibs.gen.95

# generate a histogram
awk -v size=0.02 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i]/2 }' <(paste temp.gen.2 temp.gen.1 | awk 'BEGIN{FS=OFS="\t"}{if($1!=$2 && $3!="nan")print $3}') > allStudies.nonAmb.geno.distance.mibs.hist


# Extract the within study duplicates
cat ibs.gen.95 | awk 'BEGIN{FS=OFS="\t"}{gsub("\\.",FS,$1);gsub("\\.",FS,$2);print}' > ibs.gen.95.ids
cat ibs.gen.95.ids | awk -F"\t" '{if($1==$3)print}' > ibs.gen.95.sameST
awk 'BEGIN{OFS="\t"}FNR==NR{a[$3"."$4]=$1"|"$2;next}{if($1=="50k_Gmel1")print a[$1"."$2],a[$3"."$4],$5}' 50K_Equcab2_FranchesMontagnes_Gmel.newIds ibs.gen.95.sameST > 50k_Gmel1.duplicate
awk 'BEGIN{OFS="\t"}FNR==NR{a[$3"."$4]=$1"|"$2;next}{if($1=="50k_McCue1")print a[$1"."$2],a[$3"."$4],$5}' Belgian_IlluminaSNP50_EquCab2.newIds ibs.gen.95.sameST > 50k_McCue1.duplicate
awk 'BEGIN{OFS="\t"}FNR==NR{a[$3"."$4]=$1"|"$2;next}{if($1=="50k_McCue4")print a[$1"."$2],a[$3"."$4],$5}' QH_IlluminaSNP50_EquCab2.newIds ibs.gen.95.sameST > 50k_McCue4.duplicate

## remove duplicates
cat ibs.gen.95 | cut -f1 | sort | uniq | sed 's/\./|/' | tr '|' '\t' > sample_prune.list2
plink --bfile allStudies.nonAmb --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --remove sample_prune.list --make-bed --out allStudies.nonAmb.dedup #634542 variants and 8465 samples pass filters and QC.
##########################
## Annotation
# calc the freq and generate a histogram
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq -out allInfo  ## allInfo.frq
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 allInfo.frq | awk '{print $5}') > allInfo.frq.hist

# calc the count and merge with the freq output
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq counts -out allInfo ## allInfo.frq.counts
paste allInfo.frq allInfo.frq.counts | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$11,$6}' > allInfo.frq.merged

# add hwe info
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --hardy -out allInfo ## allInfo.hwe
paste allInfo.frq.merged allInfo.hwe | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$13,$14,$15,$16}' > allInfo.frq.merged2

# add info per study
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq --family -out allInfo ## allInfo.frq.strat
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 allInfo.frq.strat | awk '{print $6}') > allInfo.frq.strat.hist
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 allInfo.frq.strat | awk '{if($6>0.5)$6=1-$6;print $6}') > allInfo.frq.strat.hist2

echo $(head -n1 allInfo.frq.merged2) "st" "st_0.1" "st_0.2" "st_0.3" "st_0.4" | tr ' ' '\t' > allInfo.frq.merged3
awk 'BEGIN{OFS="\t"}FNR==NR{x=$1 FS $2; if($8)a[x]+=1; if($6>0.5)$6=1-$6; if($6>0.1)b[x]+=1;if($6>0.2)c[x]+=1;if($6>0.3)d[x]+=1;if($6>0.4)e[x]+=1;next}{i=$1 FS $2;print $0,(a[i]>0?a[i]:0),(b[i]>0?b[i]:0),(c[i]>0?c[i]:0),(d[i]>0?d[i]:0),(e[i]>0?e[i]:0);}' <(tail -n+2 allInfo.frq.strat) <(tail -n+2 allInfo.frq.merged2) >> allInfo.frq.merged3

# add legacy info
cat allInfo.frq.strat | awk '{if($3~"^50k_" && $8>0)print $2}' | sort | uniq > all50k.list
cat allInfo.frq.strat | awk '{if($3~"^70k_" && $8>0)print $2}' | sort | uniq > all70k.list
echo "SNP legacy" | tr ' ' '\t' > all_legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]="50k";next}{if(a[$1])a[$1]="50_70k";else a[$1]="70k"}END{for (i in a)print i,a[i]}' all50k.list all70k.list >> all_legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"none";}' all_legacy.list allInfo.frq.merged3 > allInfo.frq.merged4

# add breed MAF info
# TB and TB.Jpn are treated as 1 breed. Also, AR and AR.Per are treated as 1 breed. Moreover, 54 animals in 670k_Pet1 with unknown breed labels are treated as TB
cat allStudies.nonAmb.dedup.fam | awk '{fam=$2;gsub("\\..*_","_",fam);gsub("_.*","",fam);if($1=="670k_Pet1")fam="TB";ind=$1"-"$2;$1=fam;$2=ind;print}' > allStudies.nonAmb.dedup.br.fam
cp allStudies.nonAmb.dedup.bim allStudies.nonAmb.dedup.br.bim
cp allStudies.nonAmb.dedup.bed allStudies.nonAmb.dedup.br.bed
plink --bfile allStudies.nonAmb.dedup.br --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq --family -out allInfo.br ## allInfo.br.frq.strat
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 allInfo.br.frq.strat | awk '{print $6}') > allInfo.br.frq.strat.hist
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 allInfo.br.frq.strat | awk '{if($6>0.5)$6=1-$6;print $6}') > allInfo.br.frq.strat.hist2

echo $(head -n1 allInfo.frq.merged4) "br" "br_0.1" "br_0.2" "br_0.3" "br_0.4" | tr ' ' '\t' > allInfo.frq.merged4br
awk 'BEGIN{OFS="\t"}FNR==NR{x=$1 FS $2; if($8)a[x]+=1; if($6>0.5)$6=1-$6; if($6>0.1)b[x]+=1;if($6>0.2)c[x]+=1;if($6>0.3)d[x]+=1;if($6>0.4)e[x]+=1;next}{i=$1 FS $2;print $0,(a[i]>0?a[i]:0),(b[i]>0?b[i]:0),(c[i]>0?c[i]:0),(d[i]>0?d[i]:0),(e[i]>0?e[i]:0);}' <(tail -n+2 allInfo.br.frq.strat) <(tail -n+2 allInfo.frq.merged4) >> allInfo.frq.merged4br
##########################
# select candidates
## 1) exclude variants with missing calling rate > 0.05 (i.g. genotyping rate > 0.95 to keep SNPs across arrays) or minor allele frequency > 0.3
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --geno 0.05 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.cand1
#605271 variants removed due to missing genotype data (--geno).
#3687 variants removed due to Hardy-Weinberg exact test.
#15125 variants removed due to minor allele threshold(s)
#10459 variants and 8465 samples pass filters and QC.


## 2) Pruning
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --indep 10 2 2 --out allStudies.nonAmb.dedup.cand1.toPrune
#Pruning complete.  3055 of 10459 variants removed.

## 3) Rank & label the candidates
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$2])print $0;}' allStudies.nonAmb.dedup.cand1.toPrune.prune.in allInfo.frq.merged4br | sort -k15,15nr | awk '{print $2,NR}' > cand1.rank
cat allStudies.nonAmb.dedup.cand1.toPrune.prune.out | awk '{print $1,"prune";}' >> cand1.rank

echo $(head -n1 allInfo.frq.merged4br) "select" | tr ' ' '\t' > allInfo.frq.merged5
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"not_in_cand_list"}' cand1.rank <(tail -n+2 allInfo.frq.merged4br) >> allInfo.frq.merged5
## file for the candidate SNPs only even if excluded by pruning
cat allInfo.frq.merged5 | awk 'BEGIN{FS=OFS="\t"}{if($23!="not_in_cand_list")print}' > allInfo.frq.merged5.cand
##########################
# compare to previous lists
rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs/media-1.txt $HOME/Horse_parentage_SNPs/backup_original
awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8 FS $9]=$1;next}{if(a[$5 FS $6])$17=a[$5 FS $6];else if($3)$17=$3;else $17="N/A"; print $0;}' $map670 <(cat ../backup_original/media-1.txt | sed -e "s/\r//g" | sed '/^[[:space:]]*$/d') | cut -f1-17 | sed 's/SNP Name - Affymetrix 670K$/SNP/'> media-1-TM.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$0;next}{if(a[$17])print $0,a[$17];else print $0,"not_in_1ry_list";}' allInfo.frq.merged5 media-1-TM.txt > media-1-TM_ext.txt
tail -n+2 media-1-TM_ext.txt | rev | cut -f 1 | rev | sort | uniq -c > select.list ## selected SNPs
awk -v size=1000 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(head -n-3 select.list | awk '{print $2}') > select.list.hist ## distribution of selected SNPs

tail -n+2 media-1-TM_ext.txt | cut -f40 | grep prune | wc -l
tail -n+2 media-1-TM_ext.txt | cut -f18 | sort | uniq -c
cat media-1-TM_ext.txt | awk -F"\t" '{if($18=="not_in_1ry_list" && $17!="N/A")print $17}' | sed 's/_DUP.*//' | grep -Fwf - $map670
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34!="50_70k")print $40}' | uniq -c
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $24<16095)print $40}' | uniq -c
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $24>16095 && ($28+0)<=1e-50)print $40}' | uniq -c
tail -n+2 media-1-TM_ext.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $24>16095 && ($28+0)>1e-50 && $40=="not_in_cand_list")print $0}'

## 282 are among my selected SNPs (85 in the 1st 1000, and 69,50,36,23,13,5,1 in the 2nd,3rd,... and 8th 1000)
## 92 are among my selected SNPs but excluded for pruning
## 59 not in 670k array (including 13 markers that exist but failed to remap to EquCab3)
## 91 does not show up in both 50k and 70k arrays
## 28 have genotyping rates < 95%
## 25 failed the HWE testing
########################################
# sel 1500 targets
head -n1 allInfo.frq.merged5 > sel_1500.list
tail -n+2 allInfo.frq.merged5 | awk -F"\t" '{if($23<=1500)print}' | sort -k23,23n >> sel_1500.list

## Task: upload cand1.frq.merged5
#rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/allInfo.frq.merged5.format remote_UCDavis_GoogleDr:Horse_parentage_share ## made manually
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/allInfo.frq.merged5 remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/allInfo.frq.merged5.cand remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/media-1-TM_ext.txt remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/sel_1500.list remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning
#######################################
## Assess similarities
## prep data set after removal of pruned SNPs
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --extract allStudies.nonAmb.dedup.cand1.toPrune.prune.in --make-bed --out allStudies.nonAmb.dedup.cand1.pruned

## remove extra-chr & make indvidual id unique
plink --bfile allStudies.nonAmb.dedup.cand1.pruned --chr-set 31 no-xy --allow-extra-chr 0 --allow-no-sex --recode --out allStudies.nonAmb.dedup.cand1.pruned.forEigensoft
cp allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped_temp
#cat allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped_temp | awk '{$2=$1"."$2;$6=1;print $0}' > allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped
cat allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped_temp | awk '{br=$2;gsub("_.*","",br);$2=$1"."$2;$1=br;$6=1;print $0}' > allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped
rm allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped_temp

## generate stats of breed frequency
cat allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped | cut -d" " -f1 | sort | uniq -c | sort -k1,1nr > br.freq
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/br.freq remote_UCDavis_GoogleDr:Horse_parentage_share

## pca again
#plink --file allStudies.nonAmb.dedup.cand1.pruned.forEigensoft  --pca header tabs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --out allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.pca ## keep failing

## try eigensoft
conda install -c bioconda eigensoft
cp $scripts/spca_parfile .
smartpca -p spca_parfile > smartpca.log
#ploteig -i out.evec -c 1:2 -x 1vs2.xtxt

##visualize
cp $scripts/plot.R .
cat allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.ped | \
    awk '{if(!a[$1]){n+=1;a[$1]=n;} \
    if(a[$1]<=18)x=a[$1];else x=(a[$1]-1)%18+1; \
    c=int((a[$1]-1)/18);if(c==0)col="black";else if(c==1)col="red";else if(c==2)col="green";else col="blue";
    print $1,$2,x,col}' > st_sample.map
cat st_sample.map | awk '{print $1,$3,$4}' | sort | uniq > st.dat

Rscript plot.R "PC1" "PC2"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC2.out.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

## Top left
#FM      ==> 1,circle                        black
## Top right (merge into AR)
#TB      ==> 6,triangle point down           black
#TB.Jap  ==> 14,square and triangle down     black
#Pet1.UN ==> 12,square plus                  green
## Bottom
#STB     ==> 5,diamond                       black

tail -n+2 out.evec | awk 'BEGIN{OFS="\t"}{print $1,$2}' | sort -k2,2gr | grep -v "\.TB_" | grep -v "\.TB.Jpn_" | grep -v "\Pet1.UN_" | head
#80k_KWPN1.WB_29         0.0210          9,diamond plus                          Red
#50k_McCue7.STB_405      0.0209          5,diamond                               Black
#50k_McCue1.BE_70        0.0205          2,triangle point up                     Black

tail -n+2 out.evec | awk 'BEGIN{OFS="\t"}{print $1,$3}' | sort -k2,2g | grep -v "\.STB_" | head
#670k_Bro1.AR_380        -0.0217         13,circle cross                         Black
#50k_McCue1.BE_55        -0.0192         2,triangle point up                     Black
#670k_Bro1.AR_379        -0.0187         13,circle cross (perfect overlap)       Black
#80k_Bel6.AP_174         -0.0121         11,11,triangles up and down             Black

Rscript plot.R "PC1" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC3.out.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

## Middle left
#BE      ==> 2,triangle point up             Black
#CL      ==> 10,circle plus                  Red

#ICH     ==> 5,diamond                       Red

## Bottom left
#EX      ==> 10,circle plus                  Green

## Bottom (merge into KD & )
#AR      ==> 13,circle cross                 Black
#AR.Per  ==> 17, filled triangle point-up    Red

Rscript plot.R "PC1" "PC4"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC4.out.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

Rscript plot.R "PC1" "PC5"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC5.out.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

Rscript plot.R "PC2" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC2.PC3.out.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

#########
## pairwise distance
plink --bfile allStudies.nonAmb.dedup.cand1.pruned --distance square ibs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --out allStudies.nonAmb.dedup.cand1.pruned.dist

dist="allStudies.nonAmb.dedup.cand1.pruned.dist.mibs"
cat $dist.id | awk 'BEGIN{FS="\t";OFS="|"}{gsub("\\..*_","_",$2);print $1,$2}' > $dist.id2
paste <(cat $dist.id | awk 'BEGIN{FS=OFS="\t";}{gsub("\\..*_","_",$2);print $1,$2}') $dist | grep -v "UN_" > $dist.comp

echo "study sample score match" | tr ' ' '\t' > br_bestMatch
while read study sample;do
 id="$study|$sample"
 found=$(grep -wn "$id" $dist.id2)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 newBreed=$(cut -f 1,2,$((index+2)) $dist.comp | grep -v $study$'\t'$sample | sort -k3,3nr | head -n10 | cut -f2 | sed 's/_.*//' | sort | uniq -c | sort -k1,1nr | head -n1 | sed 's/ *//')
 echo $study $sample "$newBreed" | tr ' ' '\t'
 fi
done < <(cat *.newIds | awk 'BEGIN{FS=OFS="\t";}{gsub("\\..*_","_",$4);print $3,$4}') >> br_bestMatch

echo "Breed expMatch expMatchAveScore difMatch difMatchAveScore" | tr ' ' '\t' > br_bestMatch_Summary
tail -n+2 br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);\
if($2==$4){a[$2]+=1;aSc[$2]+=$3;}else{b[$2]+=1;bSc[$2]+=$3;}}
END{for(i in a){aM=sprintf("%.2f",aSc[i]/a[i]); if(b[i]){bC=b[i];bM=sprintf("%.2f",bSc[i]/b[i]);}else{bC=bM=0;} print i,a[i],aM,bC,bM;};\
    for(x in b){if(!a[x]){aC=aM=0;bM=sprintf("%.2f",bSc[x]/b[x]);print x,aC,aM,b[x],bM;}}
}' | sort -k2,2nr >> br_bestMatch_Summary

echo "count query match" | tr ' ' '\t' > br_bestMatch_difSummary
tail -n+2 br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);if($2!=$4)print $2,$4}' | sort | uniq -c >> br_bestMatch_difSummary

head -n1 br_bestMatch > br_bestMatch_AP
cat br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{q=$2;gsub("_.*","",q);m=$4;if(q!=m && q=="AP")print $0}' >> br_bestMatch_AP

echo "study breed matchCount AveScore" | tr ' ' '\t' > br_samMatch_aveScore
cat br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);if($2==$4){a[$1 FS $2]+=$3;b[$1 FS $2]+=1}}END{for(i in a)print i,b[i],a[i]/b[i];}' | sort -k2,2 >> br_samMatch_aveScore

#echo "study qurey_breed matching_breed matchCount AveScore" | tr ' ' '\t' > br_difMatch_aveScore
#cat br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);if($2!=$4){a[$1 FS $2 FS $4]+=$3;b[$1 FS $2 FS $4]+=1}}END{for(i in a)print i,b[i],a[i]/b[i];}' | sort -k2,2 -k1,1 > br_difMatch_aveScore

echo "study qurey_breed matching_breed matchCount AveScore" | tr ' ' '\t' > br_bestMatch_aveScore
cat br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);a[$1 FS $2 FS $4]+=$3;b[$1 FS $2 FS $4]+=1}END{for(i in a)print i,b[i],a[i]/b[i];}' | sort -k2,2 -k1,1 >> br_bestMatch_aveScore

rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/br_bestMatch_Summary remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/br_bestMatch_aveScore remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

## bad/good ratios; hi > 7; Lo < 4; med 4-7
## Many breeds have no problems: MH, ICH, HF, MRM, FR, HOK, IHD, NOR, TAI, ARC, CO, MIS, NOM, SP, KIS, LU, PF, FB, ABB,
## Error labeling (hi/hi): STB(1/1174), FM(2/1100), AR(20/567), EX(5/279), BE(6/73),
## hetero animals (Lo_med/hi): TB(2med/1003hi), QH(30var/381hi), WP(4var/255hi), KNB(13var,87hi), CL(24med/51hi), MIY(1Lo/35hi), POA(7Lo/30hi),
## hetero breeds (Lo_med/Lo_med): AP(81var/171med), IRD(1Lo,20med), MOR(6Lo/18med), BIH(2Lo/17med), WU(6Lo/16med), WUZ(7Lo/14med), SN(14Lo/9Lo),
## Small sample size: HN(5med/1Lo), HO(4Hi/0), OL(5med/0), OS(4med/0), TK(5med/0), WF(5med/0),
## Special cases:
# a) WB(45med/1395hi) & SWB(57med/324hi)
# SWB/WB
#670k_Mik1       SWB     SWB     324     7.57407
#670k_Mik1       SWB     TB      7       6
#670k_Mik1       SWB     WB      49      6.81633
#80k_Vit1        SWB     WB      1       8
#80k_KWPN1       WB      AR      1       7
#80k_KWPN1       WB      FR      7       6.85714
#80k_KWPN1       WB      HF      1       5
#80k_KWPN1       WB      KNB     1       3
#80k_KWPN1       WB      SWB     28      5.78571
#80k_KWPN1       WB      TB      6       6.16667
#80k_KWPN1       WB      WB      1395    8.98638
#80k_KWPN1       WB      WP      1       7

# b) KD(57med/1med), FT(11hi/0), ISH (15med/0),
#670k_Pet2       KD      AR      57      6.80702
#670k_Pet2       KD      KD      1       5

#MNEc2M_McCue2   FT      AR      6       7.5
#MNEc2M_McCue2   FT      STB     3       7.66667
#MNEc2M_McCue2   FT      WB      2       8.5

#670k_Hill3      ISH     FM      1       3
#670k_Hill3      ISH     IRD     2       4
#670k_Hill3      ISH     SWB     1       5
#670k_Hill3      ISH     TB      5       5.8
#670k_Hill3      ISH     WB      6       6.66667


# c) UN
#670k_Pet1       UN      TB      54      9.88889
#MNEc2M_McCue2   UN      HF      6       10
#MNEc2M_McCue2   UN      NOR     4       9.25
#MNEc2M_McCue2   UN      BE      4       7.75
#MNEc2M_McCue2   UN      LU      2       7.5
#MNEc2M_McCue2   UN      WB      1       5
#MNEc2M_McCue2   UN      FM      1       3
#MNEc2M_McCue2   UN      ICH     1       3
#MNEc2M_McCue2   UN      WP      1       3
#MNEc2M_McCue2   UN      WUZ     1       3

#############################################################
## Replicate with aggressive pruning
## 2) Pruning
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --indep 100 5 1.5 --out allStudies.nonAmb.dedup.cand1.toPrune2
## Pruning complete.  5571 of 10459 variants removed.

## 3) Rank & label the candidates
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$2])print $0;}' allStudies.nonAmb.dedup.cand1.toPrune2.prune.in allInfo.frq.merged4br | sort -k15,15nr | awk '{print $2,NR}' > cand2.rank
cat allStudies.nonAmb.dedup.cand1.toPrune2.prune.out | awk '{print $1,"prune";}' >> cand2.rank

echo $(head -n1 allInfo.frq.merged4br) "select" | tr ' ' '\t' > allInfo.frq.merged5v2
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"not_in_cand_list"}' cand2.rank <(tail -n+2 allInfo.frq.merged4br) >> allInfo.frq.merged5v2
## file for the candidate SNPs only even if excluded by pruning
cat allInfo.frq.merged5v2 | awk 'BEGIN{FS=OFS="\t"}{if($23!="not_in_cand_list")print}' > allInfo.frq.merged5v2.cand
tail -n+2 allInfo.frq.merged5v2.cand | awk -F"\t" '{if($23!="prune")print $1}' | sort | uniq -c | sort -k2,2n > allInfo.frq.merged5v2.cand.chr.dist
##########################
# compare to previous lists
#rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs/media-1.txt $HOME/Horse_parentage_SNPs/backup_original
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{if($8)a[$8 FS $9]=$1;next}{if(a[$5 FS $6])$17=a[$5 FS $6];else if($3)$17=$3;else $17="N/A"; print $0;}' $map670 <(cat ../backup_original/media-1.txt | sed -e "s/\r//g" | sed '/^[[:space:]]*$/d') | cut -f1-17 | sed 's/SNP Name - Affymetrix 670K$/SNP/'> media-1-TM.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$0;next}{if(a[$17])print $0,a[$17];else print $0,"not_in_1ry_list";}' allInfo.frq.merged5v2 media-1-TM.txt > media-1-TM_extv2.txt
tail -n+2 media-1-TM_extv2.txt | rev | cut -f 1 | rev | sort | uniq -c > select.listv2 ## selected SNPs
awk -v size=1000 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(head -n-3 select.listv2 | awk '{print $2}') > select.listv2.hist ## distribution of selected SNPs

tail -n+2 media-1-TM_extv2.txt | cut -f40 | grep prune | wc -l
tail -n+2 media-1-TM_extv2.txt | cut -f18 | sort | uniq -c
cat media-1-TM_extv2.txt | awk -F"\t" '{if($18=="not_in_1ry_list" && $17!="N/A")print $17}' | sed 's/_DUP.*//' | grep -Fwf - $map670
tail -n+2 media-1-TM_extv2.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34!="50_70k")print $40}' | uniq -c
tail -n+2 media-1-TM_extv2.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $24<16095)print $40}' | uniq -c
tail -n+2 media-1-TM_extv2.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $24>16095 && ($28+0)<=1e-50)print $40}' | uniq -c
tail -n+2 media-1-TM_extv2.txt | awk -F"\t" '{if($18!="not_in_1ry_list" && $34=="50_70k" && $24>16095 && ($28+0)>1e-50 && $40=="not_in_cand_list")print $0}'

## 217 are among my selected SNPs (89 in the 1st 1000, and 65,36,19,5 in the 2nd,3rd,4th and 5th 1000)
## 158 are among my selected SNPs but excluded for pruning
## 59 not in 670k array (including 13 markers that exist but failed to remap to EquCab3)
## 91 does not show up in both 50k and 70k arrays
## 28 have genotyping rates < 95%
## 25 failed the HWE testing

########################################
# sel 1500 targets
head -n1 allInfo.frq.merged5v2 > sel_1500.listv2
tail -n+2 allInfo.frq.merged5v2 | awk -F"\t" '{if($23<=1500)print}' | sort -k23,23n >> sel_1500.listv2

# sel 1300 targets
head -n1 allInfo.frq.merged5v2 > sel_1300.listv2
tail -n+2 allInfo.frq.merged5v2 | awk -F"\t" '{if($23<=1300)print}' | sort -k23,23n >> sel_1300.listv2
#tail -n+2 sel_1300.listv2 | awk -F"\t" '{if($23!="prune")print $1}' | sort | uniq -c | sort -k2,2n > sel_1300.listv2.chr.dist


## Task: upload cand1.frq.merged5
#rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/allInfo.frq.merged5.format remote_UCDavis_GoogleDr:Horse_parentage_share ## made manually
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/allInfo.frq.merged5v2 remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/allInfo.frq.merged5v2.cand remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/allInfo.frq.merged5v2.cand.chr.dist remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/media-1-TM_extv2.txt remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/sel_1500.listv2 remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/sel_1300.listv2 remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
#rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/sel_1300.listv2.chr.dist remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning


########################################
## Assess similarities
## prep data set after removal of pruned SNPs
plink --bfile allStudies.nonAmb.dedup.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --extract allStudies.nonAmb.dedup.cand1.toPrune2.prune.in --make-bed --out allStudies.nonAmb.dedup.cand1.pruned2

## remove extra-chr & make indvidual id unique
plink --bfile allStudies.nonAmb.dedup.cand1.pruned2 --chr-set 31 no-xy --allow-extra-chr 0 --allow-no-sex --recode --out allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft
cp allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped_temp
#cat allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped_temp | awk '{$2=$1"."$2;$6=1;print $0}' > allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped
cat allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped_temp | awk '{br=$2;gsub("_.*","",br);$2=$1"."$2;$1=br;$6=1;print $0}' > allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped
rm allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped_temp

## generate stats of breed frequency
#cat allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped | cut -d" " -f1 | sort | uniq -c | sort -k1,1nr > br.freq2
#rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/br.freq2 remote_UCDavis_GoogleDr:Horse_parentage_share

## pca again
#plink --file allStudies.nonAmb.dedup.cand1.pruned.forEigensoft  --pca header tabs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --out allStudies.nonAmb.dedup.cand1.pruned.forEigensoft.pca ## keep failing

## try eigensoft
conda install -c bioconda eigensoft
cp $scripts/spca_parfile_pruned2 .
smartpca -p spca_parfile_pruned2 > smartpca_pruned2.log
#ploteig -i out.evec -c 1:2 -x 1vs2.xtxt

##visualize
cp $scripts/plot_pruned2.R .
cat allStudies.nonAmb.dedup.cand1.pruned2.forEigensoft.ped | \
    awk '{if(!a[$1]){n+=1;a[$1]=n;} \
    if(a[$1]<=18)x=a[$1];else x=(a[$1]-1)%18+1; \
    c=int((a[$1]-1)/18);if(c==0)col="black";else if(c==1)col="red";else if(c==2)col="green";else col="blue";
    print $1,$2,x,col}' > st_sample2.map  ## indifferent from st_sample.map
cat st_sample2.map | awk '{print $1,$3,$4}' | sort | uniq > st2.dat ## indifferent from st.dat

Rscript plot_pruned2.R "PC1" "PC2"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC2.out2.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
Rscript plot_pruned2.R "PC1" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC3.out2.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
Rscript plot_pruned2.R "PC1" "PC4"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC4.out2.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
Rscript plot_pruned2.R "PC1" "PC5"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC1.PC5.out2.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
Rscript plot_pruned2.R "PC2" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/PC2.PC3.out2.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning

#########
## pairwise distance
plink --bfile allStudies.nonAmb.dedup.cand1.pruned2 --distance square ibs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --out allStudies.nonAmb.dedup.cand1.pruned2.dist

dist="allStudies.nonAmb.dedup.cand1.pruned2.dist.mibs"
cat $dist.id | awk 'BEGIN{FS="\t";OFS="|"}{gsub("\\..*_","_",$2);print $1,$2}' > $dist.id2
paste <(cat $dist.id | awk 'BEGIN{FS=OFS="\t";}{gsub("\\..*_","_",$2);print $1,$2}') $dist | grep -v "UN_" > $dist.comp

echo "study sample score match" | tr ' ' '\t' > br_bestMatchv2
while read study sample;do
 id="$study|$sample"
 found=$(grep -wn "$id" $dist.id2)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 newBreed=$(cut -f 1,2,$((index+2)) $dist.comp | grep -v $study$'\t'$sample | sort -k3,3nr | head -n10 | cut -f2 | sed 's/_.*//' | sort | uniq -c | sort -k1,1nr | head -n1 | sed 's/ *//')
 echo $study $sample "$newBreed" | tr ' ' '\t'
 fi
done < <(cat *.newIds | awk 'BEGIN{FS=OFS="\t";}{gsub("\\..*_","_",$4);print $3,$4}') >> br_bestMatchv2

echo "Breed expMatch expMatchAveScore difMatch difMatchAveScore" | tr ' ' '\t' > br_bestMatchv2_Summary
tail -n+2 br_bestMatchv2 | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);\
if($2==$4){a[$2]+=1;aSc[$2]+=$3;}else{b[$2]+=1;bSc[$2]+=$3;}}
END{for(i in a){aM=sprintf("%.2f",aSc[i]/a[i]); if(b[i]){bC=b[i];bM=sprintf("%.2f",bSc[i]/b[i]);}else{bC=bM=0;} print i,a[i],aM,bC,bM;};\
    for(x in b){if(!a[x]){aC=aM=0;bM=sprintf("%.2f",bSc[x]/b[x]);print x,aC,aM,b[x],bM;}}
}' | sort -k2,2nr >> br_bestMatchv2_Summary


echo "count query match" | tr ' ' '\t' > br_bestMatchv2_difSummary
tail -n+2 br_bestMatchv2 | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);if($2!=$4)print $2,$4}' | sort | uniq -c >> br_bestMatchv2_difSummary

head -n1 br_bestMatchv2 > br_bestMatchv2_AP
cat br_bestMatchv2 | awk 'BEGIN{FS=OFS="\t"}{q=$2;gsub("_.*","",q);m=$4;if(q!=m && q=="AP")print $0}' >> br_bestMatchv2_AP

echo "study breed matchCount AveScore" | tr ' ' '\t' > br_samMatchv2_aveScore
cat br_bestMatchv2 | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);if($2==$4){a[$1 FS $2]+=$3;b[$1 FS $2]+=1}}END{for(i in a)print i,b[i],a[i]/b[i];}' | sort -k2,2 >> br_samMatchv2_aveScore

#echo "study qurey_breed matching_breed matchCount AveScore" | tr ' ' '\t' > br_difMatch_aveScore
#cat br_bestMatch | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);if($2!=$4){a[$1 FS $2 FS $4]+=$3;b[$1 FS $2 FS $4]+=1}}END{for(i in a)print i,b[i],a[i]/b[i];}' | sort -k2,2 -k1,1 > br_difMatch_aveScore

echo "study qurey_breed matching_breed matchCount AveScore" | tr ' ' '\t' > br_bestMatchv2_aveScore
cat br_bestMatchv2 | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);a[$1 FS $2 FS $4]+=$3;b[$1 FS $2 FS $4]+=1}END{for(i in a)print i,b[i],a[i]/b[i];}' | sort -k2,2 -k1,1 >> br_bestMatchv2_aveScore

rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/br_bestMatchv2_Summary remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/newIds/br_bestMatchv2_aveScore remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning


#############################################################
## Replicate with generic pruning (Useless. Therefore, I did not rerun after fixing the deduplication bug)
## 2) Pruning
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --indep 50 5 1.5 --out allStudies.nonAmb.dedup.toPrune
#Pruning complete.  459604 of 634542 variants removed.

## Assess similarities
## prep data set after removal of pruned SNPs
plink --bfile allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --extract allStudies.nonAmb.dedup.toPrune.prune.in --make-bed --out allStudies.nonAmb.dedup.pruned
#174938 variants and 8471 samples pass filters and QC.

## pairwise distance
plink --bfile allStudies.nonAmb.dedup.pruned --distance square ibs --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --out allStudies.nonAmb.dedup.pruned.dist

dist="allStudies.nonAmb.dedup.pruned.dist.mibs"
cat $dist.id | awk 'BEGIN{FS="\t";OFS="|"}{gsub("\\..*_","_",$2);print $1,$2}' > $dist.id2
paste <(cat $dist.id | awk 'BEGIN{FS=OFS="\t";}{gsub("\\..*_","_",$2);print $1,$2}') $dist | grep -v "UN_" > $dist.comp

echo "study sample score match" | tr ' ' '\t' > br_bestMatchv3
while read study sample;do
 id="$study|$sample"
 found=$(grep -wn "$id" $dist.id2)
 if [ "$found" ];then
 index=$(echo $found | cut -d":" -f1)
 newBreed=$(cut -f 1,2,$((index+2)) $dist.comp | grep -v $study$'\t'$sample | sort -k3,3nr | head -n10 | cut -f2 | sed 's/_.*//' | sort | uniq -c | sort -k1,1nr | head -n1 | sed 's/ *//')
 echo $study $sample "$newBreed" | tr ' ' '\t'
 fi
done < <(cat *.newIds | awk 'BEGIN{FS=OFS="\t";}{gsub("\\..*_","_",$4);print $3,$4}') >> br_bestMatchv3

echo "Breed expMatch expMatchAveScore difMatch difMatchAveScore" | tr ' ' '\t' > br_bestMatchv3_Summary
tail -n+2 br_bestMatchv3 | awk 'BEGIN{FS=OFS="\t"}{gsub("_.*","",$2);\
if($2==$4){a[$2]+=1;aSc[$2]+=$3;}else{b[$2]+=1;bSc[$2]+=$3;}}
END{for(i in a){aM=sprintf("%.2f",aSc[i]/a[i]); if(b[i]){bC=b[i];bM=sprintf("%.2f",bSc[i]/b[i]);}else{bC=bM=0;} print i,a[i],aM,bC,bM;};\
    for(x in b){if(!a[x]){aC=aM=0;bM=sprintf("%.2f",bSc[x]/b[x]);print x,aC,aM,b[x],bM;}}
}' | sort -k2,2nr >> br_bestMatchv3_Summary


##############
## Y markers
mkdir "$work_dir"/y_chr && cd "$work_dir"/y_chr
rclone -v --copy-links copy remote_UCDavis_GoogleDr:Horse_parentage_SNPs/Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality.txt $HOME/Horse_parentage_SNPs/backup_original/
cp $HOME/Horse_parentage_SNPs/backup_original/Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality.txt .
cp $HOME/Horse_parentage_SNPs/backup_original/GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map .

head -n3 Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality.txt | tail -n1 | cut -f1-16 > Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality_simplified.txt
tail -n+6 Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality.txt | cut -f1-16 >> Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality_simplified.txt

#tail -n+2 Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality_simplified.txt | awk 'BEGIN{FS=OFS="\t"}{print $13}' | sort | uniq -c
#tail -n+2 Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality_simplified.txt | awk 'BEGIN{FS=OFS="\t"}{print $14}' | sort | uniq -c
#tail -n+2 Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality_simplified.txt | awk 'BEGIN{FS=OFS="\t"}{if($13=="ok" && $14=="Yes_OK")print $2}' | sort > Y_SNP-List.ids
#grep ^EMSY GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map | cut -f2 | sort > 80k_emsy.ids
#comm -13 80k_emsy.ids Y_SNP-List.ids

cat <(echo "chr snpID mapPos coordinate" | tr ' ' '\t') <(head -n1 Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality_simplified.txt) | paste - - > Y_SNP-List.txt
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$0;next}{if(a[$2])print a[$2],$0}' <(grep ^EMSY GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.map | sed -e "s/\r//g") Y_SNP-List_paternityCHIP_with_KASP_Ill70K_quality_simplified.txt >> Y_SNP-List.txt
head -n1  Y_SNP-List.txt >  Y_SNP-List_sel.txt
cat Y_SNP-List.txt | awk 'BEGIN{FS=OFS="\t"}{if($17=="ok" && $18=="Yes_OK")print}' >>  Y_SNP-List_sel.txt
rclone -v --copy-links copy Y_SNP-List_sel.txt remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected

## summary: the map file has 170 variant & the original Y_SNP-List has 40 good variant
## merging the 2 files, we have 155, with only 39 from the good variants in the  Y_SNP-List


##############
## X markers
mkdir "$work_dir"/x_chr && cd "$work_dir"/x_chr


## A) PAR
plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 700000 --to-bp 2071050 --make-bed --out allStudies.nonAmb.dedup.parAll
#595 variants and 8465 samples pass filters and QC.

plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 700000 --to-bp 2071050 --geno 0.05 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.par
#10 variants and 8465 samples pass filters and QC.

plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 700000 --to-bp 2071050 --geno 0.1 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.par90
#13 variants and 8465 samples pass filters and QC.

plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 700000 --to-bp 2071050 --geno 0.2 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.par80
#14 variants and 8465 samples pass filters and QC.

plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 700000 --to-bp 2071050 --geno 0.3 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.par70
#14 variants and 8465 samples pass filters and QC.

plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 700000 --to-bp 2071050 --geno 0.4 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.par60
#18 variants and 8465 samples pass filters and QC.


plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 700000 --to-bp 2071050 --geno 0.4 --maf 0.2 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.par60maf2
#20 variants and 8465 samples pass filters and QC.


head -n1 ../newIds/allInfo.frq.merged5.cand > PAR_sel.txt
cat allStudies.nonAmb.dedup.par.bim | cut -f2 | grep -Fwf - ../newIds/allInfo.frq.merged5.cand >> PAR_sel.txt
rclone -v --copy-links copy PAR_sel.txt remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

head -n1 ../newIds/allInfo.frq.merged5v2.cand > PAR_selv2.txt
cat allStudies.nonAmb.dedup.par.bim | cut -f2 | grep -Fwf - ../newIds/allInfo.frq.merged5v2.cand >> PAR_selv2.txt
rclone -v --copy-links copy PAR_selv2.txt remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning


head -n1 ../newIds/allInfo.frq.merged5.cand > PAR60maf2_sel.txt
cat allStudies.nonAmb.dedup.par60maf2.bim | cut -f2 | grep -Fwf - ../newIds/allInfo.frq.merged5 >> PAR60maf2_sel.txt
rclone -v --copy-links copy PAR60maf2_sel.txt remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

head -n1 ../newIds/allInfo.frq.merged5v2.cand > PAR60maf2_selv2.txt
cat allStudies.nonAmb.dedup.par60maf2.bim | cut -f2 | grep -Fwf - ../newIds/allInfo.frq.merged5v2 >> PAR60maf2_selv2.txt
rclone -v --copy-links copy PAR60maf2_selv2.txt remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected


## B) X specific
plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 2130000 --to-bp 78677038 --make-bed --out allStudies.nonAmb.dedup.1.xSp ## 15461 variants and 8465 samples pass filters and QC (genotyping rate is 0.366897).
plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 78722677 --to-bp 95845791 --make-bed --out allStudies.nonAmb.dedup.2.xSp ## 3487 variants and 8465 samples pass filters and QC (genotyping rate is 0.374464)
plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --chr X --from-bp 95861454 --make-bed --out allStudies.nonAmb.dedup.3.xSp ## 6997 variants and 8465 samples pass filters and QC (genotyping rate is 0.370756)

for x in 1 2 3;do
echo allStudies.nonAmb.dedup.$x.xSp;done > x_studies_to_merge
plink --merge-list x_studies_to_merge --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --make-bed --out allStudies.nonAmb.dedup.xSp
## 25945 variants and 8465 samples pass filters and QC (genotyping rate is 0.368812).

# remove cases with very low X genotyping rate
plink --bfile allStudies.nonAmb.dedup.xSp --chr-set 32 --allow-no-sex --het --out xSp
awk -v size=259 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 xSp.het | awk '{print $5}') > xSp.het.hist ## the bin is 1% of the input variants
plink --bfile allStudies.nonAmb.dedup.xSp --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --mind 0.95 --make-bed --out allStudies.nonAmb.dedup.xSp.mind ## 25945 variants and 8191 samples pass filters and QC.

# select females based on the rate of heterozygosity
# The plink command computes observed and expected autosomal homozygous genotype counts for each sample, and reports method-of-moments F coefficient estimates (i.e. (<observed hom. count> - <expected count>) / (<total observations> - <expected count>))
plink --bfile allStudies.nonAmb.dedup.xSp.mind --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --check-sex --out xSp
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 xSp.sexcheck | awk '{print $6}') > xSp.sexcheck.hist
tail -n+2 xSp.sexcheck | awk 'BEGIN{OFS="\t"}{if($6<0.5)print $1,$2}' > female.list
plink --bfile allStudies.nonAmb.dedup.xSp.mind --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --keep female.list --make-bed --out allStudies.nonAmb.dedup.xSp.female ## 25945 variants and 4671 samples pass filters and QC.
##########################
## Annotation
# calc the freq and generate a histogram
plink --bfile allStudies.nonAmb.dedup.xSp.female --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq -out xSp.allInfo  ## allInfo.frq
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 xSp.allInfo.frq | awk '{print $5}') > xSp.allInfo.frq.hist

# calc the count and merge with the freq output
plink --bfile allStudies.nonAmb.dedup.xSp.female --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq counts -out xSp.allInfo ## allInfo.frq.counts
paste xSp.allInfo.frq xSp.allInfo.frq.counts | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$11,$6}' > xSp.allInfo.frq.merged

# add hwe info
plink --bfile allStudies.nonAmb.dedup.xSp.female --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --hardy -out xSp.allInfo ## allInfo.hwe
paste xSp.allInfo.frq.merged xSp.allInfo.hwe | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$13,$14,$15,$16}' > xSp.allInfo.frq.merged2

# add info per study
plink --bfile allStudies.nonAmb.dedup.xSp.female --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq --family -out xSp.allInfo ## allInfo.frq.strat
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 xSp.allInfo.frq.strat | awk '{print $6}') > xSp.allInfo.frq.strat.hist
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 xSp.allInfo.frq.strat | awk '{if($6>0.5)$6=1-$6;print $6}') > xSp.allInfo.frq.strat.hist2

echo $(head -n1 xSp.allInfo.frq.merged2) "st" "st_0.1" "st_0.2" "st_0.3" "st_0.4" | tr ' ' '\t' > xSp.allInfo.frq.merged3
awk 'BEGIN{OFS="\t"}FNR==NR{x=$1 FS $2; if($8)a[x]+=1; if($6>0.5)$6=1-$6; if($6>0.1)b[x]+=1;if($6>0.2)c[x]+=1;if($6>0.3)d[x]+=1;if($6>0.4)e[x]+=1;next}{i=$1 FS $2;print $0,(a[i]>0?a[i]:0),(b[i]>0?b[i]:0),(c[i]>0?c[i]:0),(d[i]>0?d[i]:0),(e[i]>0?e[i]:0);}' <(tail -n+2 xSp.allInfo.frq.strat) <(tail -n+2 xSp.allInfo.frq.merged2) >> xSp.allInfo.frq.merged3

# add legacy info
cat xSp.allInfo.frq.strat | awk '{if($3~"^50k_" && $8>0)print $2}' | sort | uniq > xSp.all50k.list
cat xSp.allInfo.frq.strat | awk '{if($3~"^70k_" && $8>0)print $2}' | sort | uniq > xSp.all70k.list
echo "SNP legacy" | tr ' ' '\t' > xSp.all_legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]="50k";next}{if(a[$1])a[$1]="50_70k";else a[$1]="70k"}END{for (i in a)print i,a[i]}' xSp.all50k.list xSp.all70k.list >> xSp.all_legacy.list
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"none";}' xSp.all_legacy.list xSp.allInfo.frq.merged3 > xSp.allInfo.frq.merged4

# add breed MAF info
# TB and TB.Jpn are treated as 1 breed. Also, AR and AR.Per are treated as 1 breed. Moreover, 54 animals in 670k_Pet1 with unknown breed labels are treated as TB
cat allStudies.nonAmb.dedup.xSp.female.fam | awk '{fam=$2;gsub("\\..*_","_",fam);gsub("_.*","",fam);if($1=="670k_Pet1")fam="TB";ind=$1"-"$2;$1=fam;$2=ind;print}' > allStudies.nonAmb.dedup.xSp.female.br.fam
cp allStudies.nonAmb.dedup.xSp.female.bim allStudies.nonAmb.dedup.xSp.female.br.bim
cp allStudies.nonAmb.dedup.xSp.female.bed allStudies.nonAmb.dedup.xSp.female.br.bed
plink --bfile allStudies.nonAmb.dedup.xSp.female.br --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --freq --family -out xSp.allInfo.br ## allInfo.br.frq.strat ## 40 clusters defined.
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 xSp.allInfo.br.frq.strat | awk '{print $6}') > xSp.allInfo.br.frq.strat.hist
awk -v size=0.05 'BEGIN{bmin=bmax=0}{ b=int($1/size); a[b]++; bmax=b>bmax?b:bmax; bmin=b<bmin?b:bmin } END { for(i=bmin;i<=bmax;++i) print i*size,(i+1)*size,a[i] }' <(tail -n+2 xSp.allInfo.br.frq.strat | awk '{if($6>0.5)$6=1-$6;print $6}') > xSp.allInfo.br.frq.strat.hist2

echo $(head -n1 xSp.allInfo.frq.merged4) "br" "br_0.1" "br_0.2" "br_0.3" "br_0.4" | tr ' ' '\t' > xSp.allInfo.frq.merged4br
awk 'BEGIN{OFS="\t"}FNR==NR{x=$1 FS $2; if($8)a[x]+=1; if($6>0.5)$6=1-$6; if($6>0.1)b[x]+=1;if($6>0.2)c[x]+=1;if($6>0.3)d[x]+=1;if($6>0.4)e[x]+=1;next}{i=$1 FS $2;print $0,(a[i]>0?a[i]:0),(b[i]>0?b[i]:0),(c[i]>0?c[i]:0),(d[i]>0?d[i]:0),(e[i]>0?e[i]:0);}' <(tail -n+2 xSp.allInfo.br.frq.strat) <(tail -n+2 xSp.allInfo.frq.merged4) >> xSp.allInfo.frq.merged4br
##########################
# select candidates
## 1) exclude variants with missing calling rate > 0.05 (i.g. genotyping rate > 0.95 to keep SNPs across arrays) or minor allele frequency > 0.3
plink --bfile allStudies.nonAmb.dedup.xSp.female --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --geno 0.05 --maf 0.3 --hwe 1e-50 --make-bed --out allStudies.nonAmb.dedup.xSp.female.cand1
## 24699 variants removed due to missing genotype data (--geno).
## 24 variants removed due to Hardy-Weinberg exact test.
## 841 variants removed due to minor allele threshold(s)
## 381 variants and 4671 samples pass filters and QC.

## 2) Pruning
plink --bfile allStudies.nonAmb.dedup.xSp.female.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --indep 10 2 2 --out allStudies.nonAmb.dedup.xSp.female.cand1.toPrune
#Pruning complete.  90 of 381 variants removed.(genotyping rate is 0.990579)

## 3) Rank & label the candidates
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$2])print $0;}' allStudies.nonAmb.dedup.xSp.female.cand1.toPrune.prune.in xSp.allInfo.frq.merged4br | sort -k15,15nr | awk '{print $2,NR}' > xSp.cand1.rank
cat allStudies.nonAmb.dedup.xSp.female.cand1.toPrune.prune.out | awk '{print $1,"prune";}' >> xSp.cand1.rank

echo $(head -n1 xSp.allInfo.frq.merged4br) "select" | tr ' ' '\t' > xSp.allInfo.frq.merged5
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"not_in_cand_list"}' xSp.cand1.rank <(tail -n+2 xSp.allInfo.frq.merged4br) >> xSp.allInfo.frq.merged5
## file for the candidate SNPs only even if excluded by pruning
cat xSp.allInfo.frq.merged5 | awk 'BEGIN{FS=OFS="\t"}{if($23!="not_in_cand_list")print}' > xSp.allInfo.frq.merged5.cand
########################################
# sel 150 targets
head -n1 xSp.allInfo.frq.merged5 > xSp.sel_150.list
tail -n+2 xSp.allInfo.frq.merged5 | awk -F"\t" '{if($23<=150)print}' | sort -k23,23n >> xSp.sel_150.list

## Task: upload cand1.frq.merged5
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/x_chr/xSp.allInfo.frq.merged5 remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/x_chr/xSp.allInfo.frq.merged5.cand remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/x_chr/xSp.sel_150.list remote_UCDavis_GoogleDr:Horse_parentage_share/Moderate_pruning

#############################################################
## Replicate with aggressive pruning
## 2) Pruning
plink --bfile allStudies.nonAmb.dedup.xSp.female.cand1 --chr-set 31 no-xy --allow-extra-chr --allow-no-sex --indep 100 5 1.5 --out allStudies.nonAmb.dedup.xSp.female.cand1.toPrune2
## Pruning complete.  168 of 381 variants removed.

## 3) Rank & label the candidates
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=1;next}{if(a[$2])print $0;}' allStudies.nonAmb.dedup.xSp.female.cand1.toPrune2.prune.in xSp.allInfo.frq.merged4br | sort -k15,15nr | awk '{print $2,NR}' > xSp.cand2.rank
cat allStudies.nonAmb.dedup.xSp.female.cand1.toPrune2.prune.out | awk '{print $1,"prune";}' >> xSp.cand2.rank

echo $(head -n1 xSp.allInfo.frq.merged4br) "select" | tr ' ' '\t' > xSp.allInfo.frq.merged5v2
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{if(a[$2])print $0,a[$2];else print $0,"not_in_cand_list"}' xSp.cand2.rank <(tail -n+2 xSp.allInfo.frq.merged4br) >> xSp.allInfo.frq.merged5v2
## file for the candidate SNPs only even if excluded by pruning
cat xSp.allInfo.frq.merged5v2 | awk 'BEGIN{FS=OFS="\t"}{if($18!="not_in_cand_list")print}' > xSp.allInfo.frq.merged5v2.cand
########################################
# sel 150 targets
head -n1 xSp.allInfo.frq.merged5v2 > xSp.sel_150.listv2
tail -n+2 xSp.allInfo.frq.merged5v2 | awk -F"\t" '{if($23<=150)print}' | sort -k23,23n >> xSp.sel_150.listv2

## Task: upload cand1.frq.merged5
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/x_chr/xSp.allInfo.frq.merged5v2 remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/x_chr/xSp.allInfo.frq.merged5v2.cand remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/x_chr/xSp.sel_150.listv2 remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
########################################
## Make one final list (using aggressively pruned data)
mkdir "$work_dir"/final && cd "$work_dir"/final
#echo "SNP_ID SNP_name Chr Pos Strand Flank REF ALT" | tr ' ' '\t' > Equ_Parentv2.manifest
#awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$1FS$2FS$8FS$9FS"."FS"."FS$10FS$11;next}{if(a[$2])print a[$2];}' $map670 <(cat "$work_dir"/newIds/sel_1300.listv2 "$work_dir"/x_chr/PAR_selv2.txt "$work_dir"/x_chr/xSp.sel_150.listv2) >> Equ_Parentv2.manifest

echo "##fileformat=VCFv4.3" > Equ_Parentv2.vcf
cat $equCab3_vcfContigs >> Equ_Parentv2.vcf
echo "#CHROM POS ID REF ALT QUAL FILTER INFO" | tr ' ' '\t' >> Equ_Parentv2.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$8FS$9FS$1FS$10FS$11FS"."FS"."FS".";next}{if(a[$2])print a[$2];}' $map670 <(cat "$work_dir"/newIds/sel_1300.listv2 "$work_dir"/x_chr/PAR_selv2.txt "$work_dir"/x_chr/xSp.sel_150.listv2) >> Equ_Parentv2.vcf
tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel.txt | awk 'BEGIN{FS=OFS="\t"}{print $8,$9}' | tr 'TCGA' 'AGCT' > temp_oppStrand
paste <(tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel.txt) temp_oppStrand | awk 'BEGIN{FS=OFS="\t"}{if($20=="norm")print "eMSYv3",$4,$7,$8,$9,".",".",".";else print "eMSYv3",$4,$7,$21,$22,".",".",".";}' >> Equ_Parentv2.vcf
bcftools sort -o Equ_Parentv2_sorted.vcf Equ_Parentv2.vcf
bcftools norm -c ws -f $equCab3_ref Equ_Parentv2_sorted.vcf 1> Equ_Parentv2.check.vcf 2> Equ_Parentv2.bcftools.log


grep -v "^#" Equ_Parentv2.check.vcf > temp_check
grep -v "^#" Equ_Parentv2_sorted.vcf > temp_orig
diff temp_check temp_orig
diff temp_check temp_orig | awk '{if($4)print $4}' | sort | uniq | grep -Fwf - "$work_dir"/y_chr/Y_SNP-List_sel.txt | tr '[' '\t' | tr ']' '\t' > dif1
cat dif1 | cut -f15 | tr 'TCGA' 'AGCT' | rev > dif2
paste dif1 dif2 | awk 'BEGIN{FS=OFS="\t"}{if($22=="norm")print $17;else print $23}' | while read marker;do
grep $marker ../eMSY/eMSYv3_nontrimmed.fasta | sed "s/$marker/\n/" | head -n1 | tr -d '\n' | wc -c;done > dif3
paste dif1 dif3 | cut -f7,15-17,21,22,23
#fYR, fTI, fWM (strand swap)
#fWQ also suffers from strand swap. Moreover, it should map to eMSYv3:5042449 instead of eMSYv3:5042499
marker=$(echo "AAGTGATGACATTGTGTTTTGCTCTAATTCCAACAGCCTTAGATACACTT" | tr 'TCGA' 'AGCT' | rev)
grep $marker ../eMSY/eMSYv3_nontrimmed.fasta | sed "s/$marker/\n/" | head -n1 | tr -d '\n' | wc -c
#rAL also suffers from strand swap. Moreover, it should be labeled "reco" instead of "norm"


## correct the position of fWQ and orintation of rAL (bcftools will swap the alleles of the 5 SNPs)
cat "$work_dir"/y_chr/Y_SNP-List_sel.txt | awk 'BEGIN{FS=OFS="\t"}{if($7=="fWQ"){$4="5042449";$19="eMSYv3:5042449"};if($7=="rAL")$20="reco";print;}' > "$work_dir"/y_chr/Y_SNP-List_sel_cor.txt
echo "##fileformat=VCFv4.3" > Equ_Parentv2cor.vcf
cat $equCab3_vcfContigs >> Equ_Parentv2cor.vcf
echo "#CHROM POS ID REF ALT QUAL FILTER INFO" | tr ' ' '\t' >> Equ_Parentv2cor.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$8FS$9FS$1FS$10FS$11FS"."FS"."FS".";next}{if(a[$2])print a[$2];}' $map670 <(cat "$work_dir"/newIds/sel_1300.listv2 "$work_dir"/x_chr/PAR_selv2.txt "$work_dir"/x_chr/xSp.sel_150.listv2) >> Equ_Parentv2cor.vcf
tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel_cor.txt | awk 'BEGIN{FS=OFS="\t"}{print $8,$9}' | tr 'TCGA' 'AGCT' > temp_oppStrand_cor
paste <(tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel_cor.txt) temp_oppStrand_cor | awk 'BEGIN{FS=OFS="\t"}{if($20=="norm")print "eMSYv3",$4,$7,$8,$9,".",".",".";else print "eMSYv3",$4,$7,$21,$22,".",".",".";}' >> Equ_Parentv2cor.vcf
bcftools sort -o Equ_Parentv2cor_sorted.vcf Equ_Parentv2cor.vcf
bcftools norm -c ws -f $equCab3_ref Equ_Parentv2cor_sorted.vcf 1> Equ_Parentv2cor.check.vcf 2> Equ_Parentv2cor.bcftools.log


rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final/Equ_Parentv2cor.check.vcf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning
#########
## check for distribution (using final list from aggressively pruned data)
grep -v "^#" Equ_Parentv2cor.check.vcf | awk -F"\t" '{print $1}' | sort | uniq -c | sort -k2,2n > Equ_Parentv2cor.check.chr.dist

rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final/Equ_Parentv2cor.check.chr.dist remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning

Rscript rainfall.R

##########
## Make one final list (using aggressively pruned data) but include the more relaxed X SNPs
mkdir "$work_dir"/final2 && cd "$work_dir"/final2

# define how many autosomal markers to keep
sex_markers=$(($(wc -l "$work_dir"/x_chr/PAR60maf2_selv2.txt "$work_dir"/x_chr/xSp.sel_150.listv2 "$work_dir"/y_chr/Y_SNP-List_sel.txt | tail -n1 | awk '{print $1}')-3))
keep=$((1500-sex_markers))
head -n1 "$work_dir"/newIds/allInfo.frq.merged5v2 > sel_auto.listv2
tail -n+2 "$work_dir"/newIds/allInfo.frq.merged5v2 | awk -v keep=$keep -F"\t" '{if($23<=keep)print}' | sort -k23,23n >> sel_auto.listv2 ## make sure this file does not have any X markers

# make the combined VCF
echo "##fileformat=VCFv4.3" > Equ_Parentv2.vcf
cat $equCab3_vcfContigs >> Equ_Parentv2.vcf
echo "#CHROM POS ID REF ALT QUAL FILTER INFO" | tr ' ' '\t' >> Equ_Parentv2.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$8FS$9FS$1FS$10FS$11FS"."FS"."FS".";next}{if(a[$2])print a[$2];}' $map670 <(cat sel_auto.listv2 "$work_dir"/x_chr/PAR60maf2_selv2.txt "$work_dir"/x_chr/xSp.sel_150.listv2) >> Equ_Parentv2.vcf
tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel.txt | awk 'BEGIN{FS=OFS="\t"}{print $8,$9}' | tr 'TCGA' 'AGCT' > temp_oppStrand
paste <(tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel.txt) temp_oppStrand | awk 'BEGIN{FS=OFS="\t"}{if($20=="norm")print "eMSYv3",$4,$7,$8,$9,".",".",".";else print "eMSYv3",$4,$7,$21,$22,".",".",".";}' >> Equ_Parentv2.vcf
bcftools sort -o Equ_Parentv2_sorted.vcf Equ_Parentv2.vcf
bcftools norm -c ws -f $equCab3_ref Equ_Parentv2_sorted.vcf 1> Equ_Parentv2.check.vcf 2> Equ_Parentv2.bcftools.log ## There are 5 problematic Markers

# Identify the problematic Markers
grep -v "^#" Equ_Parentv2.check.vcf > temp_check
grep -v "^#" Equ_Parentv2_sorted.vcf > temp_orig
diff temp_check temp_orig
diff temp_check temp_orig | awk '{if($4)print $4}' | sort | uniq | grep -Fwf - "$work_dir"/y_chr/Y_SNP-List_sel.txt | tr '[' '\t' | tr ']' '\t' > dif1
cat dif1 | cut -f15 | tr 'TCGA' 'AGCT' | rev > dif2
paste dif1 dif2 | awk 'BEGIN{FS=OFS="\t"}{if($22=="norm")print $17;else print $23}' | while read marker;do
grep $marker ../eMSY/eMSYv3_nontrimmed.fasta | sed "s/$marker/\n/" | head -n1 | tr -d '\n' | wc -c;done > dif3
paste dif1 dif3 | cut -f7,15-17,21,22,23
#fYR, fTI, fWM (strand swap)
#fWQ also suffers from strand swap. Moreover, it should map to eMSYv3:5042449 instead of eMSYv3:5042499
marker=$(echo "AAGTGATGACATTGTGTTTTGCTCTAATTCCAACAGCCTTAGATACACTT" | tr 'TCGA' 'AGCT' | rev)
grep $marker ../eMSY/eMSYv3_nontrimmed.fasta | sed "s/$marker/\n/" | head -n1 | tr -d '\n' | wc -c
#rAL also suffers from strand swap. Moreover, it should be labeled "reco" instead of "norm"


## correct the position of fWQ and orintation of rAL (bcftools will swap the alleles of the 5 SNPs)
cat "$work_dir"/y_chr/Y_SNP-List_sel.txt | awk 'BEGIN{FS=OFS="\t"}{if($7=="fWQ"){$4="5042449";$19="eMSYv3:5042449"};if($7=="rAL")$20="reco";print;}' > "$work_dir"/y_chr/Y_SNP-List_sel_cor.txt
echo "##fileformat=VCFv4.3" > Equ_Parentv2cor.vcf
cat $equCab3_vcfContigs >> Equ_Parentv2cor.vcf
echo "#CHROM POS ID REF ALT QUAL FILTER INFO" | tr ' ' '\t' >> Equ_Parentv2cor.vcf
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$8FS$9FS$1FS$10FS$11FS"."FS"."FS".";next}{if(a[$2])print a[$2];}' $map670 <(cat sel_auto.listv2 "$work_dir"/x_chr/PAR60maf2_selv2.txt "$work_dir"/x_chr/xSp.sel_150.listv2) >> Equ_Parentv2cor.vcf
tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel_cor.txt | awk 'BEGIN{FS=OFS="\t"}{print $8,$9}' | tr 'TCGA' 'AGCT' > temp_oppStrand_cor
paste <(tail -n+2 "$work_dir"/y_chr/Y_SNP-List_sel_cor.txt) temp_oppStrand_cor | awk 'BEGIN{FS=OFS="\t"}{if($20=="norm")print "eMSYv3",$4,$7,$8,$9,".",".",".";else print "eMSYv3",$4,$7,$21,$22,".",".",".";}' >> Equ_Parentv2cor.vcf
bcftools sort -o Equ_Parentv2cor_sorted.vcf Equ_Parentv2cor.vcf
bcftools norm -c ws -f $equCab3_ref Equ_Parentv2cor_sorted.vcf 1> Equ_Parentv2cor.check.vcf 2> Equ_Parentv2cor.bcftools.log


rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/Equ_Parentv2cor.check.vcf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected

#########
## check for distribution (using final list from aggressively pruned data)
grep -v "^#" Equ_Parentv2cor.check.vcf | awk -F"\t" '{print $1}' | sort | uniq -c | sort -k2,2n > Equ_Parentv2cor.check.chr.dist

rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/Equ_Parentv2cor.check.chr.dist remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected

Rscript rainfall.R ## I ran locally  ## Outputs: Rainfall_density.pdf & Rainfall-PAR.pdf

# special distribution check for PAR amrkers
bgzip -c Equ_Parentv2cor.check.vcf > Equ_Parentv2cor.check.vcf.gz
bcftools index Equ_Parentv2cor.check.vcf.gz
bcftools view --with-header --regions X:700000-2071050 -o Equ_Parentv2cor.check_PAR.vcf Equ_Parentv2cor.check.vcf.gz
echo "X 700000 PAR-start" | tr ' ' '\t' > PAR.bed
grep -v "^#" Equ_Parentv2cor.check_PAR.vcf >> PAR.bed
echo "X 2071050 PAR-end" | tr ' ' '\t' >> PAR.bed
echo "CHR,POS,ID,Dist" > PAR_dist.csv
cat PAR.bed | awk '{if(!x){x=$2;print $1,$2,$3,0;}else{print $1,$2,$3,$2-x;x=$2}}' | tr ' ' ',' >> PAR_dist.csv
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/PAR_dist.csv remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected


#######
## Assess similarities (PCA) using selected autosomal and X markers
## prep data set by extracting selected SNPs
cat Equ_Parentv2cor.check.vcf | grep -v "^#" | cut -f3 > sel_Ids
plink --bfile ../newIds/allStudies.nonAmb.dedup --chr-set 31 no-xy --allow-extra-chr --extract sel_Ids --make-bed --out allStudies.nonAmb.dedup.sel
## 1461 variants and 8465 samples pass filters and QC.

## remove extra-chr & make indvidual id unique
plink --bfile allStudies.nonAmb.dedup.sel --chr-set 31 no-xy --allow-extra-chr 0 --allow-no-sex --recode --out allStudies.nonAmb.dedup.sel.forEigensoft
cp allStudies.nonAmb.dedup.sel.forEigensoft.ped allStudies.nonAmb.dedup.sel.forEigensoft.ped_temp
cat allStudies.nonAmb.dedup.sel.forEigensoft.ped_temp | awk '{br=$2;gsub("_.*","",br);$2=$1"."$2;$1=br;$6=1;print $0}' > allStudies.nonAmb.dedup.sel.forEigensoft.ped
rm allStudies.nonAmb.dedup.sel.forEigensoft.ped_temp

## PCA by eigensoft
conda install -c bioconda eigensoft
cp $scripts/spca_parfile_sel .
smartpca -p spca_parfile_sel > smartpca_sel.log
#ploteig -i out.evec -c 1:2 -x 1vs2.xtxt

##visualize
cp $scripts/plot_sel.R .
cat allStudies.nonAmb.dedup.sel.forEigensoft.ped | \
    awk '{if(!a[$1]){n+=1;a[$1]=n;} \
    if(a[$1]<=18)x=a[$1];else x=(a[$1]-1)%18+1; \
    c=int((a[$1]-1)/18);if(c==0)col="black";else if(c==1)col="red";else if(c==2)col="green";else col="blue";
    print $1,$2,x,col}' > st_sample_sel.map  ## indifferent from st_sample.map
cat st_sample_sel.map | awk '{print $1,$3,$4}' | sort | uniq > st_sel.dat ## indifferent from st.dat

Rscript plot_sel.R "PC1" "PC2"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/PC1.PC2.out_sel.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
Rscript plot_sel.R "PC1" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/PC1.PC3.out_sel.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
Rscript plot_sel.R "PC1" "PC4"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/PC1.PC4.out_sel.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
Rscript plot_sel.R "PC1" "PC5"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/PC1.PC5.out_sel.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
Rscript plot_sel.R "PC2" "PC3"
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/PC2.PC3.out_sel.evec.pdf remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
###############
## prep GBS forms
echo "Marker ID,Chromosome,Start (Begin Location),End Location,Reference Allele,Variant Allele,Strand,Marker Type,Priority,Sequence" > GBS.csv

# Y: eMSYv3.0
#cp GBS.csv GBS_Y.csv
#tail -n+2 ../y_chr/Y_SNP-List_sel_cor.txt | awk 'BEGIN{FS="\t";OFS=","}{print $7,"eMSYv3",$4,$4,$8,$9,$20,"SNP","1",$15}' | sed 's/,norm,/,Plus,/' | sed 's/,reco,/,Minus,/' >> GBS_Y.csv
cp GBS.csv GBS_Y.csv
awk 'BEGIN{FS="\t";OFS=","}FNR==NR{a[$7]=$20",SNP,1,"$15;next}{if(a[$3])print $3,$1,$2,$2,$4,$5,a[$3]}' ../y_chr/Y_SNP-List_sel_cor.txt <(grep "^eMSYv3" Equ_Parentv2cor.check.vcf) | sed 's/,norm,/,Plus,/' | sed 's/,reco,/,Minus,/' >> GBS_Y.csv


echo "ID in_ref gen_ref in_prefix gen_prefix in_suffix gen_suffix" | tr ' ' '\t' > Y_flank_errors.tab
echo "ID in_ref gen_ref in_prefix gen_prefix in_suffix gen_suffix" | tr ' ' '\t' > GBS_Y_QC.tab
echo "ID Sequence" | tr ' ' '\t' > GBS_Y_seq.tab
seq=$(head -n2 $HOME/Horse_parentage_SNPs/eMSY/eMSYv3_nontrimmed.fasta | tail -n1 | tr 'acgt' 'ACGT')
tail -n+2 GBS_Y.csv | cut -d, -f1,2,3,5,7,10 | sed 's/\[/,\[/' | sed 's/\]/\],/' | tr ',' '\t' | while read id chr pos ref strand pref alleles suff;do
  tpos=$((pos-1))
  prefLen=${#pref}
  suffLen=${#suff}
  start=$((tpos-prefLen))
  fstart=$((tpos-100))
  tref=${seq:$tpos:1}
  if [ "$strand" == "Plus" ];then
    tpref=${seq:$start:$prefLen}
    tsuff=${seq:$pos:$suffLen}
    #tref=${seq:$tpos:1}
    fpref=${seq:$fstart:100}
    fsuff=${seq:$pos:100}
  else
    tpref=$(echo ${seq:$pos:$((suffLen-1))} | rev | tr 'ACGT' 'TGCA')
    tsuff=$(echo ${seq:$((start-1)):$((prefLen+1))} | rev | tr 'ACGT' 'TGCA')
    #tref=$(echo ${seq:$tpos:1} | tr 'ACGT' 'TGCA')
    fpref=$(echo ${seq:$pos:100} | rev | tr 'ACGT' 'TGCA')
    fsuff=$(echo ${seq:$fstart:100} | rev | tr 'ACGT' 'TGCA')
  fi
  if [[ "$ref" != "$tref" || "$pref" != "$tpref" || "$suff" != "$tsuff" ]];then echo $id $ref $tref $pref $tpref $suff $tsuff | tr ' ' '\t' >> Y_flank_errors.tab;fi
  echo $id $ref $tref $pref $tpref $suff $tsuff | tr ' ' '\t' >> GBS_Y_QC.tab
  echo $id $fpref$alleles$fsuff | tr ' ' '\t' >> GBS_Y_seq.tab
done

paste GBS_Y.csv GBS_Y_seq.tab | tr '\t' ',' | cut -d, -f-9,12 > GBS_Y_form.csv
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/GBS_Y_form.csv remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/Y_flank_errors.tab remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected


# Equcab3
cp GBS.csv GBS_nonY.csv
awk 'BEGIN{FS="\t";OFS=","}FNR==NR{a[$1]=$5",SNP,1,"$6;next}{if(a[$3])print $3,$1,$2,$2,$4,$5,a[$3]}' ../SNParrays/670k/Axiom_MNEc670.na35.r3.a2.annot.simple.csv <(grep -v "^#" Equ_Parentv2cor.check.vcf) | sed 's/,+,/,Plus,/' | sed 's/,-,/,Minus,/' >> GBS_nonY.csv

## Flank errors of illumina probes and find the 200bp flanking sequence for the GBS form
tchr=""
echo "ID in_ref gen_ref in_prefix gen_prefix in_suffix gen_suffix" | tr ' ' '\t' > nonY_flank_errors.tab
echo "ID in_ref gen_ref in_prefix gen_prefix in_suffix gen_suffix" | tr ' ' '\t' > GBS_nonY_QC.tab
echo "ID Sequence" | tr ' ' '\t' > GBS_nonY_seq.tab
tail -n+2 GBS_nonY.csv | cut -d, -f1,2,3,5,7,10 | sed 's/\[/,\[/' | sed 's/\]/\],/' | tr ',' '\t' | while read id chr pos ref strand pref alleles suff;do
  if [ "$chr" != "$tchr" ];then
    tchr=$chr;seq=$(grep -A1 "^>$chr$" $equCab3_ref_unwrap | tail -n1);fi
  tpos=$((pos-1))
  prefLen=${#pref}
  suffLen=${#suff}
  start=$((tpos-prefLen))
  tpref=$(echo ${seq:$start:$prefLen} | tr 'acgt' 'ACGT')
  tsuff=$(echo ${seq:$pos:$suffLen} | tr 'acgt' 'ACGT')
  tref=$(echo ${seq:$tpos:1} | tr 'acgt' 'ACGT')
  if [[ "$ref" != "$tref" || "$pref" != "$tpref" || "$suff" != "$tsuff" ]];then echo $id $ref $tref $pref $tpref $suff $tsuff | tr ' ' '\t' >> nonY_flank_errors.tab;fi
  echo $id $ref $tref $pref $tpref $suff $tsuff | tr ' ' '\t' >> GBS_nonY_QC.tab
  fstart=$((tpos-100))
  fpref=$(echo ${seq:$fstart:100} | tr 'acgt' 'ACGT')
  fsuff=$(echo ${seq:$pos:100} | tr 'acgt' 'ACGT')
  echo $id $fpref$alleles$fsuff | tr ' ' '\t' >> GBS_nonY_seq.tab
done

paste GBS_nonY.csv GBS_nonY_seq.tab | tr '\t' ',' | cut -d, -f-9,12 > GBS_nonY_form.csv
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/GBS_nonY_form.csv remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/nonY_flank_errors.tab remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected


cat GBS_Y_form.csv > GBS_form.csv
tail -n+2 GBS_nonY_form.csv >> GBS_form.csv
rclone -v --copy-links copy $HOME/Horse_parentage_SNPs/final2/GBS_form.csv remote_UCDavis_GoogleDr:Horse_parentage_share/agressive_pruning/selected
#######
## start a new screen
cd ~/Horse_parentage_SNPs
work_dir=$(pwd)
scripts=$work_dir/scripts
equCab3_vcfContigs="$work_dir"/equCab3/vcf_contigs.txt
map670="$work_dir/SNParrays/670k/map_670_complete.tab"
equCab3_ref="$work_dir"/equCab3/equCab3_genome.fa

####
## assessment of 80k markers
comm -23 <(grep -v "^#" Equ_Parentv2cor.check.vcf | grep -v eMSYv3 | cut -f3 | sort) <(cat ../GGPIlluminaEquineV4_Ecab2_Warmblood_KWPN.bim | cut -f2 | sort) > 80k_missing
#AX-103508102
#AX-104419981
cat 80k_missing | grep -Fwf - sel_auto.listv2
cat 80k_missing | grep -Fwf - "$work_dir"/x_chr/PAR60maf2_selv2.txt ## 2
cat 80k_missing | grep -Fwf - "$work_dir"/x_chr/xSp.sel_150.listv2
## This chip has all the Y markers but it misses 2 of PAR markers
