
cd /scratch/perickso/private/ind_seq/popgen

#use this version below!
/usr/local/sw/plink-2.0-alpha/plink2 \
--vcf /scratch/perickso/private/ind_seq/popgen/zaprionus.individual.nosingleton.2023.vcf.gz \
--king-cutoff 0.0625 \
--allow-extra-chr \
--const-fid \
--out kingunrelated0.0625

#make a file that has identical columns for family id and individual id
paste  /scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id  /scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id >  /scratch/perickso/private/ind_seq/popgen/kingunrelated0.0625.king.cutoff.in.id.use
