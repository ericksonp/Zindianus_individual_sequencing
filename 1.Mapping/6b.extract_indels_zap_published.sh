#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --time=8:00:00
#SBATCH --partition=basic
#SBATCH --array=1-5
#SBATCH --dependency=afterany:39204



#Load local software


#Name of pipeline
PIPELINE=FilterSNPs

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling

#Reference genome
REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz

#Intervals to analyze
intervals=/scratch/perickso/private/ind_seq/haplotype_calling/intervals.txt


#Parameters
#Java
JAVAMEM=45G
CPU=1

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################

#Interval to analyze
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

echo ${i} "is being processed" $(date)

###########################################################################
###########################################################################
# Filter VCFs to retain just SNPs
###########################################################################
###########################################################################

	vcftools \
		--gzvcf $WORKING_FOLDER/${i}.updated.genotyped.raw.vcf.gz \
		--recode \
		--recode-INFO-all \
		--keep-only-indels  \
		--out ${i}.updated.genotypedINDELs.raw


#filter indels

gatk --java-options "-Xmx${JAVAMEM} -Xms${JAVAMEM}" VariantFiltration \
      -R $REFERENCE \
      -V $WORKING_FOLDER/${i}.updated.genotypedINDELs.raw.recode.vcf \
			-O $WORKING_FOLDER/${i}.updated.INDEL.GATKfilter.vcf \
			--filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || InbreedingCoeff < -0.8 || ReadPosRankSum < -20.0" \
			--filter-name "GATKfilter"


vcftools --vcf $WORKING_FOLDER/${i}.updated.INDEL.GATKfilter.vcf  \
			--remove-filtered-all \
			--recode \
			--out $WORKING_FOLDER/${i}.updated.INDEL.GATKfilter.filterPASS

#make table of positions +/- 20 bp from indels

grep -v "#" $WORKING_FOLDER/${i}.updated.INDEL.GATKfilter.filterPASS.recode.vcf  | cut -f1-2    | awk '{print $1"\t"$2-20"\t"$2+20}' | awk -vOFS='\t' '$2<0 {$2=0} 1' > $WORKING_FOLDER/${i}.updated.genotypedINDELs.bed

bedtools merge -i $WORKING_FOLDER/${i}.updated.genotypedINDELs.bed > $WORKING_FOLDER/${i}.updated.genotypedINDELs.merged.bed

bedtools intersect -v \
	-header \
	-a $WORKING_FOLDER/${i}.updated.hardfilterSNP.comeault2020.vcf.gz \
	-b $WORKING_FOLDER/${i}.updated.genotypedINDELs.merged.bed > \
	$WORKING_FOLDER/${i}.updated.hardfilterSNP.comeault2020.noINDEL.vcf
