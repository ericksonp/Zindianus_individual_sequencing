WORKING_FOLDER=/scratch/perickso/private/ind_seq/haplotype_calling


cp $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter_dropsingletons.recode.vcf.gz $WORKING_FOLDER/zaprionus.individual.nosingleton.2023.vcf.gz

cp $WORKING_FOLDER/zaprionus.updated_dovetailRepeatMasked_allSNPs_noINDEL_filterPASS_norepeat_dropind_secondfilter.recode.vcf.gz $WORKING_FOLDER/zaprionus.individual.2023.vcf.gz


mv zaprionus.individual* ../popgen
