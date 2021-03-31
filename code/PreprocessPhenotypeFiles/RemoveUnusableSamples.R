###
### Several individuals need to be excluded from the dataset for various reasons
### This will generate an updated phenotype file with the unusable individuals removed, 
###   and also subset to only those individuals with exome data
### It will generate a second file that is needed for Zach's script to filter the 
###   exome vcf and get AC info
###   Format: 
###     s	project	sample	is_forbidden
###   
###
###

library(readr)
library(dplyr)
library(tidyr)

source("../code/config.R")


###
### List of cohorts and individuals that need to be removed
###
cohorts_to_remove <- c("SEARCH", "SIGMA_extremes", "TODAY")
irb_to_remove <- read.table(forbidden_52k_samples_fp, header=T, sep = "\t", stringsAsFactors=F)$Individual.ID

###
### List of exome samples
###
exome_samples <- read.table(exome_52k_samples_fp, header=T, sep = "\t", stringsAsFactors=F)

###
### 52k phenotype file
###
pheno <- read.table(get_phenotype_path("52k", flitered=FALSE), header=T, sep = ",", stringsAsFactors=F)



print(paste("Number of individuals found in exome data before removing: ",length(pheno$combined_id)))
#Number of individuals found in exome data before removing:  51667
print(paste("Number of individuals found in phenotype data before removing: ",length(exome_samples$s)))
#Number of individuals found in phenotype data before removing:  54540
print(paste("Number of individuals found in exome and phenotype data before removing: ",length(intersect(pheno$combined_id,exome_samples$s))))
#Number of individuals found in exome and phenotype data before removing:  51082



###
### remove individuals that can't be used
###
pheno <- subset(pheno, !(COHORT %in% cohorts_to_remove))
pheno <- subset(pheno, !(id_55k %in% irb_to_remove))


print(paste("Number of individuals found in exome data after removing: ",length(pheno$combined_id)))
#Number of individuals found in exome data after removing:  47499
print(paste("Number of individuals found in exome and phenotype data after removing: ",length(intersect(pheno$combined_id,exome_samples$s))))
#Number of individuals found in exome and phenotype data after removing:  46963


###
### Filter phenotype file to only those that also have exome data
###
pheno <- subset(pheno, (combined_id %in% exome_samples$s))


###
### Write out new file
###
write.table(pheno,get_phenotype_path("52k"),quote = F,sep = "\t",row.names=F)


###
### Make file for Zach's hail script to filter exome data
###
exome_samples$temp = exome_samples$s
exome_samples %>% separate(temp, c("project", "sample"), "::") -> exome_samples
exome_samples$is_forbidden = !(exome_samples$s %in% pheno$combined_id)
colnames(exome_samples) <- c("s", "project", "sample", "is_forbidden")

###
### Write out file
###
write.table(exome_samples,"../data/samples_to_keep_new_phenotype_file_V5.txt",quote = F,sep = "\t",row.names=F)

write_google_tsv(exome_samples, "gs://gnomad-zach/data/gnomad52k/samples_to_keep_new_phenotype_file_V5.txt")
