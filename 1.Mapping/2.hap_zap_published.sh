#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20G
#SBATCH --time=72:00:00
#SBATCH --partition=basic
#SBATCH --array=20
#SBATCH --dependency=afterany:39048





#BEFORE the pipeline: make dictionary file of reference genome
PICARD=/usr/local/sw/picard/picard.jar

#java -jar $PICARD  CreateSequenceDictionary R=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta.gz

#Name of pipeline

PIPELINE=Haplocaller

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq

# User defined inputs -- this represents the name of the samples
#merge all the sample names into one files
#cat /project/berglandlab/priscilla/zaps/usftp21.novogene.com/raw_data/file_list.txt /scratch/pae3g/zaps/seqs/zap_to_download.txt | grep -v ZP > /scratch/pae3g/zaps/dovetail_mapping/samples.txt
BIO_SAMPLES=/scratch/perickso/private/raw_data/published/SRR_Acc_List.txt

#Where the bam files are located
BAMS_FOLDER=$WORKING_FOLDER/joint_bams

#Reference genome
#needs to be unzipped for gatk!!
REFERENCE=/scratch/perickso/private/ref/PO1791_Zaprionus_indianus.RepeatMasked.fasta

#Sample suffixes and post-fixes. What tags are expected across all samples?
# Understanding of this comes from the previous pipeline
SUFFIX="joint.srt.rmdp"

#Parameters

#Java
JAVAMEM=18G

#Read Information
Group_library="Zaprionus2020_Dovetail_repeatMasked"
Library_Platform="illumina"
Group_platform="Novogene"

#HaploCaller -- heterozygocity prior
HET=0.005

###########################################################################
###########################################################################
# Determine sample to process, "i"
###########################################################################
###########################################################################

i=`sed -n ${SLURM_ARRAY_TASK_ID}p $BIO_SAMPLES`

###########################################################################
###########################################################################
# Begin Pipeline
###########################################################################
###########################################################################
#This part of the pipeline will generate log files to record warnings and completion status

# Welcome message
echo "your unique run id is" $unique_run_id

if [[ -e "${PIPELINE}.completion.log" ]]
then
	echo "Warning log exist"
	echo "lets move on"
	date
else
	echo "Log doesnt exist. lets fix that"
	touch $WORKING_FOLDER/${PIPELINE}.completion.log
	date
fi

# Move to working directory
cd $WORKING_FOLDER

###########################################################################
###########################################################################
# Generate Folders and files
###########################################################################
###########################################################################
# this part of the script will check and generate, if necessary, all of the output folders used in the script

echo "have you checked if the folders where already built with mkdir?"

if [[ -d "RGSM_final_bams" ]]
then
	echo "Working RGSM_final_bams folder exist"
	echo "lets move on"
	date
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/RGSM_final_bams
	date
fi


if [[ -d "haplotype_calling" ]]
then
	echo "Working haplotype_calling folder exist"
	echo "lets move on"
	date
else
	echo "folder doesnt exist. lets fix that"
	mkdir $WORKING_FOLDER/haplotype_calling
	date
fi

###########################################################################
###########################################################################
# Forcing a uniform read group to the joint bam file
###########################################################################
###########################################################################


if [[ -e "$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam" ]]
then
	echo "bam file exists"
	echo "lets move on"
	date
else
	echo "bam doesn't exist"
	java -jar $PICARD AddOrReplaceReadGroups \
		I=$BAMS_FOLDER/${i}.$SUFFIX.bam \
		O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
		RGLB=$Group_library \
		RGPL=$Library_Platform \
		RGPU=$Group_platform \
		RGSM=${i}
	date
fi


###########################################################################
###########################################################################
# Index Bam files
###########################################################################
###########################################################################
if [[ -e "$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bai" ]]
then
	echo "bai file exists"
	echo "lets move on"
	date
else
	echo "bai doesn't exist"
	java -jar $PICARD BuildBamIndex \
      I=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam
	date
fi
###########################################################################
###########################################################################
# Haplotype Calling
###########################################################################
###########################################################################
# Call haplotypes with GATK

/usr/local/sw/gatk-4.2.0.0/gatk --java-options "-Xmx${JAVAMEM}" HaplotypeCaller \
	-R $REFERENCE \
	-I $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
	-O $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf \
	--heterozygosity $HET \
	-ERC GVCF

###########################################################################
###########################################################################
# Compress and index with Tabix
###########################################################################
###########################################################################

bgzip $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf
tabix $WORKING_FOLDER/haplotype_calling/${i}.raw.g.vcf.gz

echo ${i} "completed" $(date) >> $WORKING_FOLDER/${PIPELINE}.completion.log

echo "done" $(date)
