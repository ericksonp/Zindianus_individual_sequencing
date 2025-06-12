#! /bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=erickson


#Name of pipeline
PIPELINE=admix_autosomes_noinv

#Working folder is core folder where this pipeline is being run.
WORKING_FOLDER=/scratch/perickso/private/ind_seq/popgen

#Final VCF file
IN_GZVCF=/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.vcf.gz

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
# Filter final VCF
###########################################################################
###########################################################################
if [[ -e "${PIPELINE}.plink.bim" ]]
then
	echo "Warning plink exist"
	echo "lets move on"
else
	echo "plink doesnt exist. lets fix that"
	#run admixture on unrelated females only for scaffold 3
	  plink \
	  --vcf $IN_GZVCF \
	  --recode \
		--double-id \
	  --make-bed \
	  --allow-extra-chr \
	  --out $PIPELINE.plink
fi




###########################################################################
###########################################################################
# LD Prune PLINK file
###########################################################################
###########################################################################

plink  --file $PIPELINE.plink --indep-pairwise 1000 50 0.2 --out $PIPELINE.prune --allow-extra-chr

plink --file $PIPELINE.plink --extract $PIPELINE.prune.prune.in --make-bed --out  $PIPELINE.pruneddata --allow-extra-chr

#remove "scaffold" from bim file

sed -i -e 's/Scaffold_//g' ${PIPELINE}.pruneddata.bim
admixture=/scratch/perickso/private/sw/dist/admixture_linux-1.3.0/admixture

#run admixture

IN_PLINK=$PIPELINE.pruneddata.bed

for K in  2 3 4 5 6  ; do \
$admixture $IN_PLINK $K -j10 ; done

#do CV analysis
for K in  2 3 4 5 6  ; do \
$admixture --cv=10 $IN_PLINK $K -j10  | tee log.$PIPELINE..${K}.out; done
