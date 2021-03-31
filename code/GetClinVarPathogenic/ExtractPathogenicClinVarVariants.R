library(readr)
library(dplyr)
library(lubridate)
library(tidyr) 
library(stringr)

source("../code/config.R")

# Note: This was run for both UKBB and 52K
#genome_build = "GRCh38"
#data_source = "ukbb_regeneron"


genome_build <- "GRCh37"
data_source <- "52k"


###
### Download a list of variants that appear in the curated variant list (filtered to variants in regions of interest)
###
curated_variants <- read_google_tsv(get_filtered_variant_path(data_source), col_types = 'ccciiccii')

curated_variants %>%
  separate(AC, c("Del1", "AC"), ",") %>%
  separate(AC, c("AC", "Del2"), "]") %>%
  select("variant","n_called","n_not_called","AC","AF","n_het") -> curated_variants


###
### Import the list of genes we are interested in
###
genes_list <- read_tsv(AD_gene_list_fp)
genes <- genes_list$Gene_name

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
### Identify High confidence labs
###
###
submissions %>%
  filter(str_detect(submitter_collection_method, "clinical testing")) %>% #First filter to only submissions with clinical testing
  mutate(submitter_date_last_evaluated = mdy(submitter_date_last_evaluated)) %>% #Reformat the date and time info where it exists
  group_by(submitter) %>% #Group by submitter
  mutate(Freq = n()) %>% #Count up number of submissions per submitter
  mutate(last_eval_date = max(submitter_date_last_evaluated, na.rm=TRUE)) %>% #determine the most recent submission date of each submitter
  select(submitter, Freq, last_eval_date) %>% #Select only the columns: submitter, total number of submissions, and most recent evaluation date
  distinct %>% #Reduce to a distinct set of rows since freq and last eval will be on all submissions
  arrange(desc(Freq)) %>% #sort by number of submissions, most resent first
  filter(Freq > 15000 & last_eval_date >= mdy("Jan 01, 2017")) -> submitters_HC #Filter to only labs with >15,000 submissions and updated after Jan 01, 2017


submitters_HC %>%
  rename(`Number of submissions` = Freq) %>%
  rename(`Last evaluation date` = last_eval_date) %>%
  write_tsv(HC_lab_list) -> submitters_HC


saveRDS(submitters_HC, "../output/GetClinVarPathogenic/submitters_HC.rds")


##Filter out Counsyl from the high confidence labs because of conflicting phenotype info
submitters_HC %>% 
  filter(submitter != "Counsyl")->submitters_HC

###
### Import the ClinVar vcf
###
vcf <- read_tsv(get_vcf_path(genome_build),
                comment = "#",
                col_types = 'ciiccccc',
                col_names=c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"))


###
### Extract position info from the clinvar vcf because the info in the tab-delimited file 
### might not be correctly aligned for comparing to our vcf positions 
### Split the vcf INFO field to get the ALLELEID 
### Then only keep chromosome, position, VariationID, ref allele, alt allele, and AlleleID
###
vcf %>%  
  separate(INFO, c("Del1", "INFO"), "ALLELEID=") %>%  
  separate(INFO, c("AlleleID", "del2"), ";") %>% 
  mutate(CHROM = paste0("chr", CHROM)) %>%
  select(CHROM,POS,ID,REF,ALT,AlleleID) -> map_id_position

colnames(map_id_position) <- c("chr", "pos", "VariationID", "REF", "ALT", "AlleleID")

###
### Merge clinvar tab-delimited data with the Position information from the vcf and filter
###
clinvar %>%
  right_join(map_id_position,by="VariationID") %>% #Join clinvar tab-delimited data with position from vcf
  filter(GeneSymbol %in% genes) %>% #Only keep variants in the genes of interest
  filter(ClinSigSimple == 1) %>% #Filter to only keep variants with at least one report of Likely Path or Path using the ClinSigSimple column from the tab-delimited data
  filter(Assembly == genome_build) %>% #Filter to only keep variant info on genome Assembly GRCh37
  mutate(variant_str = paste(chr, pos, REF, ALT, sep=':')) %>% #Add another column with position and allele information that can be used to join
  filter(variant_str %in% curated_variants$variant) -> test_clinvar #Filter to include only variants that are actually in our dataset


###
### Join the filtered Likely Path and path variant information with all submissions 
###
test_clinvar %>% 
  inner_join(submissions,by="VariationID") %>% #Join variant information with submissions using VariationID
  left_join(genes_list,by = c("GeneSymbol" = "Gene_name")) %>% 
  mutate(last_evaluated = mdy(LastEvaluated),
         submitter_date_last_evaluated = mdy(submitter_date_last_evaluated)) %>% #Reformat the date and time info where it exists
  group_by(VariationID) %>% #Group by VariationID
  filter(any(tolower(submitter_clinical_significance) %in% clinvar_path_levels)) %>% #Filter any Variant where there is not at least one Likely Path or Path
  arrange(desc(submitter_date_last_evaluated),.by_group = TRUE) %>%  #sort by date last evaluated
  ungroup -> sorted_submission_report



sorted_submission_report %>%
  filter(str_detect(submitter_collection_method, "clinical testing") | grepl("Lipodystrophy", Traits)  | grepl("Neonatal Diabetes", Traits)) %>%
  #filter(!(reported_phenotype_info == "CN517202:not provided" |reported_phenotype_info == "CN169374:not specified")) %>% # CANT USE THIS FILTER, EXCLUDES ALL OF GENEDX!!!!!
  group_by(VariationID) %>% #Group by VariationID
  filter(any(tolower(submitter_clinical_significance) %in% clinvar_path_levels)) %>% #Filter any Variant where there is not at least one Likely Path or Path
  ungroup -> sorted_submission_report

###
### Filter to minimal info about the remaining variants
###
sorted_submission_report %>%
  select(VariationID, GeneSymbol, chr, pos, REF, ALT, HGNC_ID, Traits) %>%
  mutate(variant_str = paste(chr, pos, REF, ALT, sep=':')) %>%
  filter(variant_str %in% curated_variants$variant) %>% 
  distinct -> exported_variant_data


###
### Join with Allele count information
###
exported_variant_data %>% 
  left_join(curated_variants,by = c("variant_str" = "variant")) %>% 
  filter(AC > 0) -> bla


###
### Write out to file
###
bla %>%
  write_excel_csv(paste0(get_data_path(data_source, "clinvar_pathogenic"), format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), ".csv"))


###
### Join with Allele count information
###
sorted_submission_report %>% 
  left_join(curated_variants,by = c("variant_str" = "variant"))%>% 
  filter(AC > 0)-> bla


###
### Write out to file
###
bla %>%
  write_excel_csv(paste0(get_data_path(data_source, "clinvar_submissions"), format(Sys.time(), "%Y-%m-%d_%H.%M.%S"), ".csv"))



