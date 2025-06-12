#post phasing file prep
cd /scratch/perickso/private/ind_seq/popgen


  for f in *.phased.vcf; do
     bgzip $f
     tabix ${f}.gz
  done
   ls $PWD/*.phased.vcf.gz > Phased_vcfs_to_merge.txt
  bcftools merge -l Phased_vcfs_to_merge.txt > zaprionus.individual.2023.all.phased.vcf
  bgzip zaprionus.individual.2023.all.phased.vcf
  tabix zaprionus.individual.2023.all.phased.vcf.gz ;
done



nohup vcf-info-annotator \
-o /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.vcf \
-d "ancestral allele" \
-f Character \
/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.vcf.gz \
/scratch/perickso/private/ind_seq/popgen/africa_ancestral_alleles.txt \
AA &




bgzip /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.vcf
tabix /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.vcf.gz

for chr in {1..5}; do
 bcftools annotate -r Scaffold_$chr -x ^INFO/AA,^FORMAT/GT /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.vcf.gz > /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.$chr.vcf
done


#doesnt seem to be needed but keeping this code just in case:
# for chr in {1..5}; do
# vcftools --vcf /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.vcf \
#   --chr Scaffold_${chr} \
#   --recode-INFO-all \
#   --out /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.2023.all.phased.ancestral.${chr}
# done
