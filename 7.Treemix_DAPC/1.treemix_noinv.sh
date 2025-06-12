
#workflow from https://speciationgenomics.github.io/Treemix/

cd /scratch/perickso/private/ind_seq/popgen/treemix

#vcf file prep=will use these files elsewhere!
#bcftools annotate --set-id '%CHROM\_%POS\' /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gz -o  /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.vcf.gz -Oz
vcftools \
--gzvcf /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf.gz \
--chr Scaffold_1 \
--chr Scaffold_2 \
--chr Scaffold_4 \
--chr Scaffold_5 \
--recode \
--stdout | gzip > /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.vcf.gz

vcftools \
--gzvcf  /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.vcf.gz \
--max-missing 1 \
--recode \
--stdout | gzip >  /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.vcf.gz


#see 0.treemix_sample_prep.R for list of samples to extract


#try again with only VA and FL samples
bcftools view \
-S /scratch/perickso/private/ind_seq/popgen/treemix/treemix_samps_CM_MIA.txt \
/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.vcf.gz \
-o  /scratch/perickso/private/ind_seq/popgen/treemix/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.vcf.gz -Oz

cd /scratch/perickso/private/ind_seq/popgen/treemix

#make plink files for LD pruning
plink \
--vcf  /scratch/perickso/private/ind_seq/popgen/treemix/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.vcf.gz \
--recode \
--make-bed \
--double-id \
--allow-extra-chr \
--out  zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink

#LD prune
plink  \
--file zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink \
--indep-pairwise 100 10 .2 \
--out zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune \
--allow-extra-chr

#extract LD-pruned SNPs
plink \
--file zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink \
--extract zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.prune.in \
--recode vcf-iid \
--out  zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.pruneddata \
--allow-extra-chr

#make treemix files
./vcf2treemix.sh zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.pruneddata.vcf /scratch/perickso/private/ind_seq/popgen/treemix/treemix_cluster_CM_MIA.txt

#run treemix
conda activate treemix

for i in {0..5}
do
 treemix -i zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.pruneddata.treemix.frq.gz \
 -m $i \
 -o zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.autosomesonly.noN.CM_MIA.plink.prune.pruneddata.treemix.frq.gz$i \
 -root MIA_2019 \
 -bootstrap \
 -k 500 \
 -noss > treemix_noinv_${i}_log &
done
