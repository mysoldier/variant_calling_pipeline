#!/usr/bin/env bash
#SBATCH --job-name=snp_calling_pipeline
#SBATCH --mail-type=END,ABORT
#SBATCH --mail-user=tpaisie@ufl.edu
#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --output=out_snp_calling_pipeline_test_log_%j
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=5gb
#SBATCH --time=90:00:00

#**********************************************#

#Author: Taylor Paisie
#Must have bam files made
#Uses bam files that were made with the galaxy pipeline
#Will give you a SNP alignment with "?" as gaps

#**********************************************#

module load intel/2016.0.109 
module load openmpi/1.10.2
module load gcc/5.2.0
module load parallel
module load gatk
module load vcflib
module load bcftools
module load samtools
module load htslib


REF=/ufrc/salemi/tpaisie/javiana/refseq/NC_020307.fa
#REF=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa

export _JAVA_OPTIONS="-Xms1g -Xmx10g"

#VCFFILE=$1

# step 2 - extract SNPs and Indels from each vcf file
# extracts SNPS
gatk SelectVariants -R $REF -V ${VCFFILE} --select-type-to-include SNP -O snps_${VCFFILE}

gatk SelectVariants -R $REF -V ${VCFFILE} --select-type-to-include INDEL -O indels_${VCFFILE}

# step 3 - filter the SNP and Indel files
# filters the SNP vcf
gatk VariantFiltration -R $REF -V snps_${VCFFILE} --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filter-name "basic_snp_filter" -O filtered_snps_${VCFFILE}

gatk VariantFiltration -R $REF -V indels_${VCFFILE} --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O filtered_indels_${VCFFILE}


# step 4 - BQSR #1
parallel 'gatk BaseRecalibrator -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}.bam --known-sites filtered_snps_javiana_uncal_041418.vcf --known-sites filtered_indels_javiana_uncal_041418.vcf -O {}.table' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)


# step 6 - Analyze the BQSR reports from the base recalibration steps (steps 4 & 5)
parallel 'gatk ApplyBQSR -R /ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa -I {}.bam --bqsr-recal-file {}.table -O recal_{}.bam' ::: $(ls *.bam | rev | cut -c 5- | rev | uniq)



LIST=$1
#calling variants on recalibrated bam files (will have only one vcf file as the output)
freebayes -L ${LIST}.txt -v ${LIST}.vcf -f $REF 0.001 -p 1 -i -X -n 0 -E 3 --min-repeat-size 5 -m 1 -q 20 -R 0 -Y 0 -e 1000 -F 0.5 -C 2 -3 0 -G 1 -! 0


# step 9 - filter vcf file from freebayes for SNPs only
vcffilter -f "TYPE = snp" ${LIST}.vcf > snps_${LIST}.vcf


# step 10 - compressing and indexing the snp only vcf file & variant normalization of the snp only vcf file & decompressing the vcf file

#compress vcf
bgzip snps_${LIST}.vcf

#index vcf
tabix -p vcf snps_${LIST}.vcf.gz

# normalizing the variant vcf file
bcftools norm -f $REF -o norm_${LIST}.vcf.gz snps_${LIST}.vcf.gz

#decompresses the output variant normalized file
bgzip -d norm_${LIST}.vcf.gz


# step 11 - filter normalized vcf file and create SNP alignment fasta file
bash vcflib_pipeline_HPC2.sh norm_${LIST}.vcf











