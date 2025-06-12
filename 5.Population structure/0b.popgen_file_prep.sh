#popgen file preparation

#add snpids to vcf
bcftools annotate --set-id '%CHROM\_%POS\' /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gz -o  /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.vcf.gz -Oz
tabix /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.vcf.gz
#filter for mac of 3 for population genetic studies
vcftools --gzvcf /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.vcf.gz --mac 3 --recode --stdout | bgzip >  /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.vcf.gz
tabix /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.vcf.gz
#make a gds from this file for use in R

bedtools intersect -v \
	-header \
	-a /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.vcf.gz \
	-b /scratch/perickso/private/ind_seq/popgen/LDdecay_inversions.bed > \
	/scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf

bgzip /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf
tabix /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.withID.mac3.removeinv.vcf.gz
