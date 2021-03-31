library(readr)
library(dplyr)
library(lubridate)
library(tidyr) 
library(stringr)

source("../code/config.R")

###
### Import the tab-del submitters file that was pre-filtered and manually annotated for 
##  matching the phenotype of interest
###

# Note this was done for both UKBB and 52k
#submissions_post_phenotype_annotation_fp = "../output/GetClinVarPathogenic/ukbb_regeneron/submissions_2019-04-10_15.44.29_phenotype_annotation.csv"
#data_source = "ukbb_regeneron"
#genome_build = "GRCh38"


submissions_post_phenotype_annotation_fp <- "../output/GetClinVarPathogenic/52k/submissions_2019-04-10_14.58.53_annotate_phenotypes.csv"
data_source <- "52k"
genome_build <- "GRCh37"


submissions_post_phenotype_annotation <- read_csv(submissions_post_phenotype_annotation_fp)



curated_variants <- read_google_tsv(get_filtered_variant_path(data_source), col_types = 'ccciiccii')

curated_variants %>%
  separate(AC, c("Del1", "AC"), ",") %>%
  separate(AC, c("AC", "Del2"), "]") %>%
  select("variant","n_called","n_not_called","AC","AF","n_het") -> curated_variants



###
### Refilter the annotated submission file for Path and Likely path
###
submissions_post_phenotype_annotation %>% 
  filter(str_detect(submitter_collection_method, "clinical testing")) %>%
  mutate(submitter_date_last_evaluated = mdy(submitter_date_last_evaluated)) %>% #Reformat the date and time info where it exists
  group_by(VariationID) %>% #Group by VariationID
  mutate(flag_HC_lab_unspecified_phenotype = any((submitter %in% high_conf_lab) &
                                                   (clarification_unspecified_reported_phenotype == "No Information") &
                                                   (submitter_date_last_evaluated >= mdy("Jan 01, 2017")),na.rm = TRUE)) %>% 
  ungroup %>%
  filter(Correct_phenotype == TRUE) %>% #filter out submissions that are the incorrect phenotype
  group_by(VariationID) %>% #Group by VariationID
  filter(any(tolower(submitter_clinical_significance) %in% clinvar_path_levels)) %>% #Filter any Variant where there is not at least one Likely Path or Path
  arrange(desc(submitter_date_last_evaluated),.by_group = TRUE) %>%  #sort by date last evaluated
  ungroup -> sorted_submission_report




###
### Add info about if the submition is from a high confidence lab after Jan 2017
###
sorted_submission_report %>% 
  mutate(high_confidence_lab_2017 = ((submitter %in% high_conf_lab) & 
                                       #(USE_HC) &
                                       (submitter_date_last_evaluated >= mdy("Jan 01, 2017")))) -> sorted_submission_report



###
### Add info about if high confidence lab submitted Path or Likely path after Jan 2017
###
sorted_submission_report %>%  
  group_by(VariationID) %>% 
  mutate(high_confidence_lab_2017_any = any(high_confidence_lab_2017 == TRUE)) %>%
  mutate(high_confidence_lab_2017_clinical_significance = paste(submitter_clinical_significance[high_confidence_lab_2017 == TRUE], collapse=", ")) %>%
  ungroup -> sorted_submission_report
  



###
### Filter to minimal info about the remaining variants
###
sorted_submission_report %>%
  select(VariationID, 
         AlleleID,
         Type,
         Name,
         GeneSymbol,
         LastEvaluated,
         PhenotypeList, 
         NumberSubmitters,
         chr, 
         pos, 
         REF, 
         ALT, 
         Traits,
         last_evaluated,
         n_called,
         n_not_called,
         AC,
         n_het,
         flag_HC_lab_unspecified_phenotype,
         high_confidence_lab_2017_any,
         high_confidence_lab_2017_clinical_significance) %>%
  mutate(clinvar_link = paste0("https://www.ncbi.nlm.nih.gov/clinvar/variation/", VariationID)) %>%
  mutate(gnomad_link = paste0("http://gnomad.broadinstitute.org/variant/", paste(chr, pos, REF, ALT, sep = '-'))) %>%
  mutate(variant_str = paste(chr, pos, REF, ALT, sep=':')) %>%
  mutate(Variant = paste(chr, pos, REF, ALT, sep='-')) %>%
  filter(variant_str %in% curated_variants$variant) %>% 
  distinct -> exported_variant_data


###
### Write pathogenic variant info out to file
###
exported_variant_data %>%
  write_excel_csv(paste0(get_data_path(data_source, "clinvar_pathogenic_pheno"), format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), ".csv"))


###
### Write pathogenic submission info out to file
###
sorted_submission_report %>%
  write_excel_csv(paste0(get_data_path(data_source, "clinvar_submissions_pheno"), format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), ".csv"))



