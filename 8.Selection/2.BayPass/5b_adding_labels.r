library(data.table)
library(ggplot2)
#library(slider)
library(foreach)
library(stringr)
#library(dpylr)
library(foreach)
library(zoo)

#
#data <- fread(file = '/scratch/perickso_shared/alexandra/BayPassAlexEdited/zaprionus_individual_2023_3a_outlier_SNPs_q999.txt')
data <- fread(file = '/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zap_ind_2023_core_summary_pi_xtx.out')
location <- fread(file = '/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/zap_ind_2023_plink_out_pruned.prune.in')
#data <- fread(file = "/Users/Home/Documents/Erickson-Lab/zaprionus_individual_2023_3a_outlier_SNPs_q999.xtx")
firstRow = colnames(location)
locations = rbind(firstRow, location, use.names = FALSE)
data$locations = locations

chr<- c()
i = 1
while(i <= nrow(data)){
 chr <-append(chr, as.integer(substr(data$locations[i], 10, 10)))
  i = i+1
  if(i%%1000 == 0)
    print(i)
} 

data$Scaffold <- chr
data$num <- rownames(data)
fwrite(data, '/scratch/perickso_shared/alexandra/BayPassAlexEdited/allVAFemalesvsFLfemales/updated_data_for_manhattan.txt')