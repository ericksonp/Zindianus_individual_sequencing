#extract CVs from admixture logs

for i in {2..8}; do
grep CV log.zaprionus.individual.nosingleton.2023.${i}.out
done > zaprionus.individual.nosingleton.2023.CVs.txt


#K=5 is best
