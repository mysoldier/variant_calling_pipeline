#!/usr/bin/env bash
#SBATCH --job-name=snp_calling_pipeline
#SBATCH --mail-type=END,ABORT
#SBATCH --mail-user=tpaisie@ufl.edu
#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --output=out_snp_calling_pipeline_test_log_%j
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=30gb
#SBATCH --time=96:00:00

#**********************************************#

#Author: Taylor Paisie
#Must have bam files made
#Uses bam files that were made with the galaxy pipeline
#Will give you a SNP alignment with "?" as gaps

#**********************************************#




module load gatk
module load samtools
module load freebayes
module load htslib
module load bcftools
module load vcflib
module load python


REF=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa

export _JAVA_OPTIONS="-Xms1g -Xmx10g"

# index bam files if they are not indexed
for i in *.bam
	do
		samtools index ${i}
	done


# step 1 - call variants from each sample using Haplotype Caller
for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
	do
		GenomeAnalysisTK -T HaplotypeCaller -R $REF -I ${i}.bam -o ${i}.vcf -ploidy 1
	done



# step 2 - extract SNPs and Indels from each vcf file
# extracts SNPS
for i in $(ls *.vcf | rev | cut -c 5- | rev | uniq)
	do 
		GenomeAnalysisTK -T SelectVariants -R $REF -V ${i}.vcf -selectType SNP -o ${i}_snps.vcf 
	done

# extracts Indels
for i in $(ls *.vcf | rev | cut -c 5- | rev | uniq)
	do
		GenomeAnalysisTK -T SelectVariants -R $REF -V ${i}.vcf -selectType INDEL -o ${i}_indels.vcf
	done



# step 3 - filter the SNP and Indel files
# filters the SNP vcf
for i in $(ls *_snps.vcf | rev | cut -c 10- | rev | uniq)
	do
		GenomeAnalysisTK -T VariantFiltration -R $REF -V ${i}_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o ${i}_snps_filtered.vcf
	done

# filters the Indel vcf
for i in $(ls *_indels.vcf | rev | cut -c 12- | rev | uniq)
	do
		GenomeAnalysisTK -T VariantFiltration -R $REF -V ${i}_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o ${i}_indels_filtered.vcf
	done



# step 4 - BQSR #1
for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
	do
		GenomeAnalysisTK -T BaseRecalibrator -R $REF -I ${i}.bam -knownSites ${i}_snps_filtered.vcf -knownSites ${i}_indels_filtered.vcf -o ${i}.table
	done


# step 5 - BQSR #2
for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
	do
		GenomeAnalysisTK -T BaseRecalibrator -R $REF -I ${i}.bam -knownSites ${i}_snps_filtered.vcf -knownSites ${i}_indels_filtered.vcf -BQSR ${i}.table -o post_${i}.table
	done



# step 6 - Analyze the BQSR reports from the base recalibration steps (steps 4 & 5)
for i in $(ls *.table | rev | cut -c 7- | rev | uniq)
	do
		GenomeAnalysisTK -T AnalyzeCovariates -R $REF -before ${i}.table -after post_${i}.table -l DEBUG -csv post_${i}.csv -plots post_${i}.pdf
		Rscript BQSR.R post_${i}.csv ${i}.table post_${i}.pdf 
	done


# make directory for original bam files
mkdir original_bams

# step 7 - applying the base recalibration scores to the bam files
for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
	do
		GenomeAnalysisTK -T PrintReads -R $REF -I ${i}.bam -BQSR ${i}.table -o ${i}_recal.bam
		mv ${i}.bam original_bams/
	done


# step 8 - rename bam files (if you want) and call variants with freebayes
#renaming recalibrated bam files
for i in $(ls *_recal.bam | rev | cut -c 11- | rev | uniq)
	do
		mv ${i}_recal.bam ${i}.bam
		mv ${i}_recal.bai ${i}.bai
	done


LIST=all_O1_inaba
#calling variants on recalibrated bam files (will have only one vcf file as the output)
freebayes -L ${LIST}.txt -v ${LIST}.vcf -f $REF -T 0.001 -p 1 -i -X -n 0 -E 3 --min-repeat-size 5 -m 1 -q 20 -R 0 -Y 0 -e 1000 -F 0.5 -C 2 -3 0 -G 1 -! 0


# step 9 - filter vcf file from freebayes for SNPs only
vcffilter -f "TYPE = snp" vch1786_allO1_fb_051717.vcf > vch1786_allO1_snps_051817.vcf


# step 10 - compressing and indexing the snp only vcf file & variant normalization of the snp only vcf file & decompressing the vcf file

#compress vcf
bgzip vch1786_allO1_snps_051817.vcf

#index vcf
tabix -p vcf vch1786_allO1_snps_051817.vcf.gz

# normalizing the variant vcf file
bcftools norm -f $REF -o vch1786_allO1_norm_051817.vcf.gz vch1786_allO1_snps_051817.vcf.gz

#decompresses the output variant normalized file
bgzip -d vch1786_allO1_norm_051817.vcf.gz



# step 11 - filter normalized vcf file and create SNP alignment fasta file
bash vcflib_pipeline_HPC2.sh vch1786_allO1_norm_051817.vcf











