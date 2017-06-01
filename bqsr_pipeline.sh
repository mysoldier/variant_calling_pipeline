#!/usr/bin/env bash
#SBATCH --job-name=bqsr_pipeline
#SBATCH --mail-type=END,ABORT
#SBATCH --mail-user=tpaisie@ufl.edu
#SBATCH --account=epi
#SBATCH --qos=epi
#SBATCH --output=out_bqsr_pipeline_test_log_%j
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=30gb
#SBATCH --time=96:00:00

#**********************************************#

#Author: Taylor Paisie
#Must have bam files made
#Uses bam files that were made with the galaxy pipeline
#only does the BQSR on bam files
#gives final, recalibrated bam files

#**********************************************#




module load gatk
module load samtools
module load R


REF=/ufrc/salemi/tpaisie/cholera/ref_seq/v_cholerae_o1_2010el_1786.fa

export _JAVA_OPTIONS="-Xms1g -Xmx10g"

# index bam files if they are not indexed
#for i in *.bam
	#do
		#samtools index ${i}
	#done


# step 1 - call variants from each sample using Haplotype Caller
#for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
	#do
		#GenomeAnalysisTK -T HaplotypeCaller -R $REF -I ${i}.bam -o ${i}.vcf -ploidy 1
	#done



# step 2 - extract SNPs and Indels from each vcf file
# extracts SNPS
#for i in $(ls *.vcf | rev | cut -c 5- | rev | uniq)
	#do 
		#GenomeAnalysisTK -T SelectVariants -R $REF -V ${i}.vcf -selectType SNP -o snps_${i}.vcf 
		#GenomeAnalysisTK -T SelectVariants -R $REF -V ${i}.vcf -selectType INDEL -o indels_${i}.vcf
	#done



# step 3 - filter the SNP and Indel files
# filters the SNP vcf
#for i in $(ls snps_*.vcf | rev | cut -c 5- | rev | uniq)
	#do
		
		#GenomeAnalysisTK -T VariantFiltration -R $REF -V ${i}.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o ${i}_filtered.vcf
	#done

# filters the Indel vcf
#for i in $(ls indels_*.vcf | rev | cut -c 5- | rev | uniq)
	#do
		#GenomeAnalysisTK -T VariantFiltration -R $REF -V ${i}.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o ${i}_filtered.vcf
	#done



# step 4 - BQSR #1
for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
	do
		GenomeAnalysisTK -T BaseRecalibrator -R $REF -I ${i}.bam -knownSites snps_${i}_filtered.vcf -knownSites indels_${i}_filtered.vcf -o ${i}.table
	done


# step 5 - BQSR #2
for i in $(ls *.bam | rev | cut -c 5- | rev | uniq)
	do
		GenomeAnalysisTK -T BaseRecalibrator -R $REF -I ${i}.bam -knownSites snps_${i}_filtered.vcf -knownSites indels_${i}_filtered.vcf -BQSR ${i}.table -o post_${i}.table
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

# gives final, recalibrated bam file
# perform freebayes on all samples after this pipeline
