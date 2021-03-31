library(readr)
library(dplyr)
library(lubridate)
library(tidyr) 
library(stringr)

source("../code/config.R")

# Note: Used for both 52K and UKBB which use different genome builds
genome_build <- "GRCh37"
#genome_build = "GRCh38"

print(paste("Using genome build: ", genome_build))
  
###
### Import the list of genes we are interested in
###
genes_list <- read_tsv(AD_gene_list_fp)
genes <- genes_list$Gene_name


### 
### Currently one line per trait per gene, collapse to one line 
### per gene comma seperate trait 
###
genes_list %>% 
  group_by(Gene_name) %>% 
  summarise(Traits = paste(Trait, collapse=", ")) -> genes_list


###
### Import the tab-del ClinVar data and filter to only clinical submissions
###
clinvar <- read_tsv(clinvar_path)
submissions <- read_tsv(submission_path,
                        comment = "",
                        skip = 16,
                        col_types = 'icccccccccccc',
                        col_names=c("VariationID",
                                   "submitter_clinical_significance",
                                   "submitter_date_last_evaluated",
                                   "submitter_description",
                                   "submitted_phenotype_info",
                                   "reported_phenotype_info",
                                   "submitter_review_status",
                                   "submitter_collection_method",
                                   "submitter_origin_counts",
                                   "submitter",
                                   "submitterSCV",
                                   "submitted_gene_symbol",
                                   "ExplanationOfInterpretation"))




###
### Merge clinvar tab-delimited data with the Position information from the vcf and filter
###
clinvar %>%
  filter(GeneSymbol %in% genes) %>% #Only keep variants in the genes of interest
  filter(ClinSigSimple == 1) %>% #Filter to only keep variants with at least one report of Likely Path or Path using the ClinSigSimple column from the tab-delimited data
  filter(Assembly == genome_build) -> test_clinvar #Filter to only keep variant info on genome Assembly GRCh37


###
### Join the filtered Likely Path and path variant information with all submissions 
###
test_clinvar %>% 
  inner_join(submissions,by="VariationID") %>% #Join variant information with submissions using VariationID
  left_join(genes_list,by = c("GeneSymbol" = "Gene_name")) %>% #Add our trait info
  mutate(last_evaluated = mdy(LastEvaluated),
         submitter_date_last_evaluated = mdy(submitter_date_last_evaluated)) %>% #Reformat the date and time info where it exists
  group_by(VariationID) %>% #Group by VariationID
  filter(any(tolower(submitter_clinical_significance) %in% clinvar_path_levels)) %>% #Filter any Variant where there is not at least one Likely Path or Path
  arrange(desc(submitter_date_last_evaluated),.by_group = TRUE) %>%  #sort by date last evaluated
  ungroup -> sorted_submission_report


###
### Keep all submissions with at least 1 path/likely path for the variant
### 
sorted_submission_report %>%
  group_by(VariationID) %>% #Group by VariationID
  filter(any(tolower(submitter_clinical_significance) %in% clinvar_path_levels)) %>% #Filter any Variant where there is not at least one Likely Path or Path
  ungroup -> sorted_submission_report_all


###
### Keep clinical submissions with at least 1 clinical path/likely path for the variant
### 
sorted_submission_report %>%
  filter(str_detect(submitter_collection_method, "clinical testing")) %>%
  #filter(!(reported_phenotype_info == "CN517202:not provided" |reported_phenotype_info == "CN169374:not specified")) %>%
  group_by(VariationID) %>% #Group by VariationID
  filter(any(tolower(submitter_clinical_significance) %in% clinvar_path_levels)) %>% #Filter any Variant where there is not at least one Likely Path or Path
  ungroup -> sorted_submission_report_clinic_only


###
### Keep non clinical submissions with at least 1 non clinical path/likely path for the variant
### 
sorted_submission_report %>%
  filter(!str_detect(submitter_collection_method, "clinical testing")) %>%
  #filter(!(reported_phenotype_info == "CN517202:not provided" |reported_phenotype_info == "CN169374:not specified")) %>%
  group_by(VariationID) %>% #Group by VariationID
  filter(any(tolower(submitter_clinical_significance) %in% clinvar_path_levels)) %>% #Filter any Variant where there is not at least one Likely Path or Path
  ungroup -> sorted_submission_report_no_clinic


###
### Summarize the numbers of submissions for all, clinical, and non clinical by gene
### 
sum_submissions_by_gene_clinic_only <- data.frame(table(sorted_submission_report_clinic_only$GeneSymbol))
sum_submissions_by_gene_no_clinic <- data.frame(table(sorted_submission_report_no_clinic$GeneSymbol))
sum_submissions_by_gene_all <- data.frame(table(sorted_submission_report_all$GeneSymbol))

sum_submissions_by_gene_all %>%
  full_join(sum_submissions_by_gene_clinic_only,by = "Var1") %>%
  full_join(sum_submissions_by_gene_no_clinic,by = "Var1") -> sum_submissions_by_gene

colnames(sum_submissions_by_gene) <- c("Gene",
                                       "sum_submissions_by_gene_all",
                                       "sum_submissions_by_gene_clinic_only",
                                       "sum_submissions_by_gene_no_clinic")


###
### Filter to minimal info about the remaining variants
### Want to get the total number of variants per gene that 
### have a clinical path/likely path submission
###

sorted_submission_report_all %>%
  select(VariationID, GeneSymbol, Traits) %>%
  distinct -> exported_variant_data_all


sorted_submission_report_clinic_only %>%
  select(VariationID, GeneSymbol, Traits) %>%
  distinct -> exported_variant_data_clinic_only


sorted_submission_report_no_clinic %>%
  select(VariationID, GeneSymbol,Traits) %>%
  distinct -> exported_variant_data_no_clinic


sum_variants_by_gene_all <- data.frame(table(exported_variant_data_all$GeneSymbol))
sum_variants_by_gene_clinic_only <- data.frame(table(exported_variant_data_clinic_only$GeneSymbol))
sum_variants_by_gene_no_clinic <- data.frame(table(exported_variant_data_no_clinic$GeneSymbol))

sum_variants_by_gene_all %>%
  full_join(sum_variants_by_gene_clinic_only,by = "Var1") %>%
  full_join(sum_variants_by_gene_no_clinic,by = "Var1") -> sum_variants_by_gene

colnames(sum_variants_by_gene) <- c("Gene", "sum_variants_by_gene_all",
                                    "sum_variants_by_gene_clinic_only",
                                    "sum_variants_by_gene_no_clinic")


###
### Keep all submissions that are path/likely path
###
sorted_submission_report %>%
  filter(tolower(submitter_clinical_significance) %in% clinvar_path_levels)  -> sorted_submission_report_all

###
### Keep all clinical testing submissions that are path/likely path
###
sorted_submission_report %>%
  filter(str_detect(submitter_collection_method, "clinical testing")) %>%
  filter(tolower(submitter_clinical_significance) %in% clinvar_path_levels)-> sorted_submission_report_clinic_only


###
### Get a list of all the phenotypes for the clinical path/likely path submissions
###
sorted_submission_report_clinic_only %>% 
  group_by(GeneSymbol) %>% 
  summarise(SubmittedTraits = paste(unique(reported_phenotype_info), collapse=", ")) %>% 
  #summarise(Traits = paste(unique(reported_phenotype_info), collapse=", ")) %>% 
  distinct -> sorted_submission_report_clinic_only_phenotypes

colnames(sorted_submission_report_clinic_only_phenotypes) <- c("Gene", "Clinical_P_LP_phenotypes")


###
### Keep all non clinical testing submissions that are path/likely path
###
sorted_submission_report %>%
  filter(!str_detect(submitter_collection_method, "clinical testing")) %>%
  filter(tolower(submitter_clinical_significance) %in% clinvar_path_levels)  -> sorted_submission_report_no_clinic

###
### Summarize the numbers of path/likely path submissions for all, clinical, and non clinical by gene
### 
sum_submissions_P_LP_by_gene_clinic_only <- data.frame(table(sorted_submission_report_clinic_only$GeneSymbol))
sum_submissions_P_LP_by_gene_no_clinic <- data.frame(table(sorted_submission_report_no_clinic$GeneSymbol))
sum_submissions_P_LP_by_gene_all <- data.frame(table(sorted_submission_report_all$GeneSymbol))

sum_submissions_P_LP_by_gene_all %>%
  full_join(sum_submissions_P_LP_by_gene_clinic_only,by = "Var1") %>%
  full_join(sum_submissions_P_LP_by_gene_no_clinic,by = "Var1") -> sum_submissions_P_LP_by_gene

colnames(sum_submissions_P_LP_by_gene) <- c("Gene",
                                            "sum_submissions_P_LP_by_gene_all",
                                            "sum_submissions_P_LP_by_gene_clinic_only",
                                            "sum_submissions_P_LP_by_gene_no_clinic")


sum_submissions_by_gene %>%
  full_join(sum_variants_by_gene,by = "Gene") %>%
  full_join(sum_submissions_P_LP_by_gene,by = "Gene") %>%
  full_join(sorted_submission_report_clinic_only_phenotypes,by = "Gene")-> sum_submissions_and_variants_by_gene


sum_submissions_and_variants_by_gene %>%
  write_excel_csv(paste0(get_sum_clinvar_sub_by_gene_prefix(genome_build), format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), ".csv"))




