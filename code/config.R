
### function to write to copy to a google cloud folder
write_google_tsv <- function(.x, remote_path) {
    local_path <- tempfile()
    try({
        write_tsv(.x, local_path)
        system(paste0("/Users/jgoodric/google-cloud-sdk/bin/gsutil cp ", local_path, " ", remote_path))
    })
    unlink(local_path)
}

read_google_tsv <- function(remote_path, ...) {
    local_path <- paste0(tempfile(), ".", tools::file_ext(remote_path))
    try({
        system(paste0("/Users/jgoodric/google-cloud-sdk/bin/gsutil cp ", remote_path, " ", local_path))
    })
    .x <- read_tsv(local_path, ...)
    unlink(local_path)
    .x
}


carrier_categories <- c("Carrier Pathogenic", "Carrier Likely Pathogenic", "Carrier VUS", "Carrier Likely Benign", "Carrier Benign", 
    "Carrier LOF", "Carrier Likely_LOF", "Carrier LOF_Unknown", "Carrier Likely_Not_LOF", "Carrier Not_LOF", "Non Carrier")

high_conf_lab <- c("Invitae", "GeneDx", "Laboratory for Molecular Medicine,Partners HealthCare Personalized Medicine", "Genetic Services Laboratory, University of Chicago", 
    "EGL Genetic Diagnostics,Eurofins Clinical Diagnostics", "PreventionGenetics,PreventionGenetics", "Ambry Genetics")

clinvar_levels <- c("pathogenic", "likely pathogenic", "vus", "likely benign", "benign")
clinvar_path_levels <- c("pathogenic", "likely pathogenic")
lof_levels <- c("lof", "likely_lof", "likely_not_lof", "not_lof")
lof_path_levels <- c("lof", "likely_lof")
path_levels <- c("pathogenic", "likely pathogenic", "lof", "likely_lof")
non_path_levels <- c("vus", "likely benign", "benign", "likely_not_lof", "not_lof")
non_path_no_vus_levels <- c("likely benign", "benign", "likely_not_lof", "not_lof")

all_levels <- c("pathogenic", "likely pathogenic", "lof", "likely_lof", "vus", "likely benign", "benign", "likely_not_lof", "not_lof")


get_filtered_variant_path <- function(data_source) {
    if (data_source == "52k") {
        return("gs://gnomad-zach/data/gnomad52k/breaking_julia_post_GQ_DP_AB_filter/filtered_variants.gz")
    }
    if (data_source == "ukbb_regeneron") {
        return("gs://gnomad-julia/UKBB/ukb_regeneron_plink.filtered_variants.gz")
    }
}

clinvar_path <- "../data/ClinVar/variant_summary_2019-04.txt.gz"
submission_path <- "../data/ClinVar/submission_summary_2019-04.txt.gz"

get_vcf_path <- function(reference) {
    vcf_path <- paste0("../data/ClinVar/clinvar_20190408_", reference, ".vcf.gz")
}
HC_lab_list <- "../output/GetClinVarPathogenic/list_of_HC_labs.txt"


clinvar_output_path <- "../output/GetClinVarPathogenic"
review_curation_path <- "../output/RadiantCuration"
carrier_output_path <- "../output/GetCarrierLists"
lof_curation_path <- "../output/LoFCuration"
gnomad_AF_path <- "../data/gnomAD_AF"

get_sum_clinvar_sub_by_gene_prefix <- function(genome_build) {
    return(paste0(clinvar_output_path, "/sum_all_submissions_and_variants_by_gene_", genome_build, "_"))
}

get_data_path <- function(data_source, data_type) {
    if (data_type == "clinvar_pathogenic") {
        return(paste(clinvar_output_path, data_source, "reported_pathogenics_", sep = "/"))
    }
    if (data_type == "clinvar_submissions") {
        return(paste(clinvar_output_path, data_source, "submissions_", sep = "/"))
    }
    if (data_type == "clinvar_pheno_annotation") {
        return(paste(clinvar_output_path, data_source, "VariantCuration_annotate_phenotype.csv", sep = "/"))
    }
    if (data_type == "clinvar_pathogenic_pheno") {
        return(paste(clinvar_output_path, data_source, "reported_pathogenics_after_phenotype_curation_", sep = "/"))
    }
    if (data_type == "clinvar_submissions_pheno") {
        return(paste(clinvar_output_path, data_source, "submissions_after_phenotype_curation_", sep = "/"))
    }
    if (data_type == "lof_carrier") {
        return(paste(carrier_output_path, data_source, "lof_carrier.txt", sep = "/"))
    }
    if (data_type == "clinvar_carrier") {
        return(paste(carrier_output_path, data_source, "clinvar_carrier.txt", sep = "/"))
    }
    if (data_type == "lof_curated") {
        return(paste(lof_curation_path, data_source, "lof_curated.txt", sep = "/"))
    }
    if (data_type == "clinvar_curated") {
        return(paste(clinvar_output_path, data_source, "clinvar_curated.txt", sep = "/"))
    }
    if (data_type == "review_carrier") {
        return(paste(carrier_output_path, data_source, "radiant_carrier.txt", sep = "/"))
    }
    if (data_type == "review_curated") {
        return(paste(review_curation_path, data_source, "radiant_curated.txt", sep = "/"))
    }
    if (data_type == "gnomad_AF") {
        return(paste(gnomad_AF_path, paste0(data_source, "_gnomad_AF.txt"), sep = "/"))
    }

}



get_phenotype_path <- function(data_source, filtered = TRUE) {
    if (data_source == "52k") {
        if (filtered) {
            return("../output/UpdatedPhenotypeFiles/FINAL_pheno_file_sequenced_onlyV5_Jul20_individuals_removed.txt")
        } else {
            return("../data/Phenotypes/FINAL_pheno_file_sequenced_onlyV5_Jul20.csv")
        }
    } else {
        if (filtered) {
            return("../data/Phenotypes/UKBiobank_genoQC_reportedANDgeneticEUR_N455146_diabetesphenotypes_complete_updateMay2019_add_lipids_pheno_interest_extra_phenotypes_unrelated_exomes.txt")
        } else {
            return("../data/Phenotypes/UKBiobank_genoQC_reportedANDgeneticEUR_N455146_diabetesphenotypes_complete_updateMay2019_add_lipids_pheno_interest_extra_phenotypes.txt")
        }
    }
}


get_related_exclude_path <- function(data_source, exclude = TRUE) {
    if (data_source == "52k") {
        if (exclude) {
            return("../data/t2d.basic.related.exclude")
        } else {
            return("../data/t2d.all.related.dat.txt")
        }
    }
}



forbidden_52k_samples_fp <- "../data/forbidden_samples.txt"
exome_52k_samples_fp <- "../data/sample_ids.txt"


AD_gene_list_fp <- "../data/AD_trait_genelist.txt"
map_gene_pheno_fp <- "../data/map_gene_to_phenotype.txt"

ancestry_order <- c("East Asian", "South Asian", "European", "African American", "Hispanic", "Other")
disease_order <- c("High LDL", "Low LDL", "High HDL", "High triglycerides", "Monogenic obesity", "MODY", "MODY extended", "Lipodystrophy", 
    "Neonatal diabetes")
threshold_for_plot <- data.frame(Trait = c("High HDL", "High LDL", "Low LDL", "High triglycerides", "Obesity"), threshold = c(70, 
    190, 80, 200, 30), stringsAsFactors = F)

variant_levels <- c("NC", "C", "Non Carrier", "Carrier", "ClinVar", "LOF", "Review", "ClinVar and LOF", "LOF and Review", "VUS Novel", 
    "VUS <1/10000", "VUS <1/1000", "VUS >=1/1000")

colors <- c(Review = "#5B9555", ClinVar = "#377eb8", LOF = "#d6604d", `Non Carrier` = "#666666", Carrier = "#d6604d", NC = "#666666", 
    C = "#d6604d", `ClinVar and LOF` = "#807dba", `LOF and Review` = "#e88c30", `VUS Novel` = "#2b8cbe", `VUS <1/10000` = "#80b1d3", 
    `VUS <1/1000` = "#fdb462", `VUS >=1/1000` = "#b3de69")

ancestry_colors <- c(Hispanic = "#ED1E24", European = "#6AA5CD", `African American` = "#941494", `South Asian` = "#FF9912", `East Asian` = "#108C44")

my_paired_palette <- c("#ecb4ac", "#d6604d", "#9dc3e2", "#377eb8", "#99C098", "#5B9555")
write_google_tsv <- function(.x, remote_path) {
    local_path <- tempfile()
    try({
        write_tsv(.x, local_path)
        system(paste0("/Users/jgoodric/google-cloud-sdk/bin/gsutil cp ", local_path, " ", remote_path))
    })
    unlink(local_path)
}

read_google_tsv <- function(remote_path, ...) {
    local_path <- paste0(tempfile(), ".", tools::file_ext(remote_path))
    try({
        system(paste0("/Users/jgoodric/google-cloud-sdk/bin/gsutil cp ", remote_path, " ", local_path))
    })
    .x <- read_tsv(local_path, ...)
    unlink(local_path)
    .x
}


get_PRS_path <- function(phenotype) {
    PRS_fp <- switch(phenotype, ldl_mgdl = "../data/PRS/LDL_PRS_ukb_v3_imp_info_0.3_inf_ENGAGE_hm3_h2_0.1347_score.txt.profile", 
        hdl_mgdl = "../data/PRS/HDL_PRS_ukb_v3_imp_info_0.3_inf_ENGAGE_hm3_h2_0.1572_score.txt.profile", tg_mgdl = "../data/PRS/TG_PRS_ukb_v3_imp_info_0.3_inf_ENGAGE_hm3_h2_0.1572_score.txt.profile", 
        tg_mgdl_ln = "../data/PRS/TG_PRS_ukb_v3_imp_info_0.3_inf_ENGAGE_hm3_h2_0.1572_score.txt.profile", bmi = "../data/PRS/penetrance_ukb_imp_chr_Khera.et.al_GPS_BMI_Cell_2019_mod_fix_import.tsv", 
        bmi_ln = "../data/PRS/penetrance_ukb_imp_chr_Khera.et.al_GPS_BMI_Cell_2019_mod_fix_import.tsv", t2d = "../data/PRS/penetrance_ukb_imp_chr_Type2Diabetes_PRS_LDpred_rho0.01_v3_fix_import.tsv")
    return(PRS_fp)
}
