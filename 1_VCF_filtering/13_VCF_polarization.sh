#!/bin/bash
#SBATCH -p shared,edwards # partition to submit to 
#SBATCH -c 2      # cores requested 
#SBATCH -t 1-00:00  # time requested in days-hours:minutes
#SBATCH -o 0_fastga.ut # name of output file
#SBATCH -e 0_fastga.err  # name of error file
#SBATCH --mem=80000 #Mb
#SBATCH --mail-type=END # email me when my job is done running

WD1=/n/netscratch/edwards_lab/Lab/gizquierdoaranega/CALONECTRIS/CALPUF
cd $WD1

module load python/3.10.9-fasrc01
conda activate fastga

## 1. Align with FastGA

FastGA -pafS GCA_987198025.1_bPufPuf1.hap1.1_genomic.fna GCA_964196065.2_bCalBor7.hap2.2_genomic.fna > $WD1/CALPUF.paf

### 2. Create VCF with paftools

conda activate myenv

sort -k6,6 -k8,8n CALPUF.paf | paftools.js call -l 1000 -L 5000 -s Pufpuf -f GCA_964196065.2_bCalBor7.hap2.2_genomic.fna - | bgzip > PUFPUF.vcf.gz
tabix -p vcf PUFPUF.vcf.gz

#### 3. Filter your target VCF to biallelic SNPs (polarizing indels or multiallelic sites is notoriously error-prone):

conda activate bcftools

ALLVCF=Calonectris.dv.15x.combined_phasing.LEURANDOM.mac2.vcf.gz
PUFVCF=CHR1.vcf.gz

#Next, merge your target VCF with the outgroup VCF so you can evaluate the outgroup's genotype at every target site:

bcftools merge $ALLVCF $PUFVCF -Oz -o Calonectris_wPuf.vcf.gz
tabix Calonectris_wPuf.vcf.gz

#bcftools view -r CHR1 Calonectris_wPuf.vcf.gz -Oz -o CHR1.wPuf.vcf.gz
#tabix CHR1.wPuf.vcf.gz

### Use script to repolarize

conda activate myenv

PYSCRIPT=/n/home11/gizquierdoaranega/scripts/CALONECTRIS/polarizeVCFbyOutgroup.allfields.py

#python3 $PYSCRIPT -vcf CHR1.wPuf.vcf.gz -out CHR1.wPuf.polarized.vcf -ind Pufpuf -keep

### 4. Finally, remove Pufpuf and keep only polymorphic positions

conda activate vcftools

vcftools --vcf CHR1.wPuf.polarized.vcf --remove-indv Pufpuf --mac 1 --recode --recode-INFO-all --out CHR1.wPuf.polarized

conda activate bcftools

bgzip CHR1.wPuf.polarized.recode.vcf
mv CHR1.wPuf.polarized.recode.vcf.gz CHR1.wPuf.polarized.vcf.gz
tabix CHR1.wPuf.polarized.vcf.gz
