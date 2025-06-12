#This script will run through a script of the populations of Z. ind we have and choose the populations we would like to analyze
#This script chooses female flies from Miami and V
#If populations change change /scratch/perickso_shared/alexandra/baypass/data/zap_full_info_updated_v2.csv to the updated .csv
#If you want to change the destination files change the write values
#All 0 slurm scripts can be ignored if the files are already present in their folders (ex. we have already come up with the FL population or added ids to the VCF, etc)

srun \
--pty \
-t 6:00:00 \
--mem=100G \
--partition erickson \
--ntasks-per-node=24 \
bash

require(ape) ; require(data.table)
library(dplyr)

choosing_areas <-data.frame(fread("/scratch/perickso_shared/alexandra/baypass/data/zap_full_info_updated_v2.csv"))

choosing_areas <- choosing_areas %>% filter(assigned_sex == "F")

#check for other databases but this filtered to just VA and MIA

MIA = choosing_areas %>% filter(Location == "MIA")
MIA_names = MIA$sample.id
MIA_names <- na.omit(MIA_names)
write(MIA_names,"/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/MIA_sample_ids.txt")

VA = choosing_areas %>% filter(Location != "MIA") #this will get only VA samples because this row in the file is NA for all samples outside of VA and MIA
VA_names = VA$sample.id
VA_names <- na.omit(VA_names)
write(VA_names,"/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/VA_sample_ids.txt")