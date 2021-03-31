library(scales)
library(reshape2)
library(plyr)
library(dplyr)

source("../code/config.R")
load_data <- function(data_source, age_cutoff, keep_ascertained=FALSE){
  if(data_source == "52k"){
    pheno <- read.table(get_phenotype_path(data_source), header=T, sep = "\t", stringsAsFactors=F)
    ###
    ### Redefine the T2D diabetes category  
    ### Fast glu > 100; 2hr > 140; HbA1c > 5.7% == pre-diabetes AND diabetes
    ###
    pheno$t2d_with_pre = pheno$t2d_2_plus_or_2hg_11_1
    pheno$t2d = pheno$t2d_2_plus_or_2hg_11_1
    pheno$t2d_with_pre[which(((pheno$fast_glu_mmoll*18)>=100) | ((pheno$h_glu*18)>=140) | (pheno$hba1c>=5.7))] = 1
    pheno$t2d[which(pheno$hba1c>=6.5)] = 1
    
    
    ###
    ### Clean up the race3 field of the phenotype file
    ###
    pheno$Ancestry = pheno$race3
    pheno %>%
      mutate(Ancestry=replace(Ancestry, Ancestry == "asian","Other")) %>%
      mutate(Ancestry=replace(Ancestry, Ancestry == "other","Other")) %>%
      mutate(Ancestry=replace(Ancestry, Ancestry == "ea","East Asian")) %>%
      mutate(Ancestry=replace(Ancestry, Ancestry == "sa","South Asian")) %>%
      mutate(Ancestry=replace(Ancestry, Ancestry == "aa","African American")) %>%
      mutate(Ancestry=replace(Ancestry, Ancestry == "eu","European")) %>%
      mutate(Ancestry=replace(Ancestry, Ancestry == "hs","Hispanic")) -> pheno

    epacts_unrelated <- read.table("../data/Phenotypes/t2d.epacts.covar.unrelated.ped.txt", header=T, sep = "\t", stringsAsFactors=F, comment.char = "", check.names = FALSE)
    
    #Filter to samples in unrelated ped file and filter out samples with phenotypes related to our traits of interest
    pheno <- pheno %>% filter(combined_id %in% epacts_unrelated$IID)
    if(!keep_ascertained){
      pheno <- pheno %>%
        filter(!(esp_phenotype %in% c("BMI_High","BMI_High_Nondiab","BMI_High_diab","BMI_Low","LDL_High","LDL_Low")))
    }
    
    #pheno = pheno %>% filter(combined_id %in% epacts_unrelated$IID)
    pheno$fast_glu = pheno$fast_glu_mmoll
    
    
    ###
    ### Load allele frequency information
    ###
    gnomad_AF <- read.delim(get_data_path(data_source, "gnomad_AF"), header=T, sep="\t", stringsAsFactors=F)
    gnomad_AF <- gnomad_AF %>% select(locus, alleles,
                                      gnomad_exomes_AF = gnomad_exomes_popmax_AF,
                                      gnomad_exomes_AC = gnomad_exomes_52K_adj_AC)
    gnomad_AF$gnomad_exomes_AF[which(is.na(gnomad_AF$gnomad_exomes_AF))] = 0
    
  }else{
    pheno <- read.table(get_phenotype_path(data_source, filtered=FALSE), header=T, sep = "\t", stringsAsFactors=F)
    exome_sample_ids <- read.table("../data/ukbb_exome_sample_ids.txt", header=FALSE, stringsAsFactors=F)
    pheno <- pheno %>% filter(FLOREZ_FID %in% exome_sample_ids[, 1])
    pheno <- pheno %>% filter(unrelated == 1)
    pheno <- pheno %>% mutate(sex=case_when(SEX == "Male" ~ 1, SEX=="Female" ~ 2),
                              age=age/12,
                              cov_GENO_ARRAY=case_when(cov_GENO_ARRAY == "UKBB" ~ 1, cov_GENO_ARRAY=="UKBL" ~ 2),
                              hba1c = hba1c.30750.NGSP.max)
    pheno["combined_id"] <- pheno["FLOREZ_IID"]
    
    ###
    ### Load allele frequency information
    ###
    
    gnomad_AF <- read.delim(get_data_path(data_source, "gnomad_AF"), header=T, sep="\t", stringsAsFactors=F)
    gnomad_AF <- gnomad_AF %>% select(locus, alleles,
                                      gnomad_exomes_AF = gnomad_exomes_popmax_AF,
                                      gnomad_exomes_AC = gnomad_exomes_AC)
    gnomad_AF$gnomad_exomes_AF[which(is.na(gnomad_AF$gnomad_exomes_AF))] = 0
    
  }
  pheno <- pheno %>% filter(age >= age_cutoff)
  pheno <- pheno %>% mutate(ldl_mgdl_adj=ifelse(is.na(lipidmeds) | !lipidmeds, ldl_mgdl, ldl_mgdl/0.7), tg_mgdl_adj=ifelse(is.na(lipidmeds) | !lipidmeds, tg_mgdl, tg_mgdl/0.85))
  
  gnomad_AF <- gnomad_AF %>%
    separate(locus, c("chr", "pos"), ":") %>%
    separate(alleles, c("ref", "alt"), ",") %>%
    separate(alt, c("alt", "Del1"), "]") %>%
    separate(ref, c("Del2", "ref"), "\\[") %>%
    mutate(Variant = paste(chr,pos,ref,alt,sep = "-"))%>%
    select(-Del1,-Del2)
  
  ###
  ### LOF carrier list
  ###
  lof <- read.delim(get_data_path(data_source, "lof_carrier"), header=T, sep="\t", stringsAsFactors=F)
  
  ###
  ### Annotation of LOF variants
  ###
  lof_to_keep <- read.delim(get_data_path(data_source, "lof_curated"), header=T, sep="\t", stringsAsFactors=F)

  ###
  ### ClinVar carrier list
  ###
  clinvar <- read.delim(get_data_path(data_source, "clinvar_carrier"), header=T, sep="\t", stringsAsFactors=F)
  
  ###
  ### Annotation of ClinVar variants
  ###
  clinvar_to_keep <- read.delim(get_data_path(data_source, "clinvar_curated"), header=T, sep="\t", stringsAsFactors=F)

  ###
  ### Review carrier list
  ###
  review <- read.delim(get_data_path(data_source, "review_carrier"), header=T, sep="\t", stringsAsFactors=F)
  
  ###
  ### Annotation of review variants
  ###
  review_to_keep <- read.delim(get_data_path(data_source, "review_curated"), header=T, sep="\t", stringsAsFactors=F)
  
  ###
  ### Load in mapping from gene to relevant disease for LoF
  ###
  AD_gene_list <- read.delim(AD_gene_list_fp, header=T, sep = "\t", stringsAsFactors=F)

  ###
  ### Load in list of genes with known LoF mechanism
  ###
  LoF_Genes <- subset(AD_gene_list, LoF == 1)$Gene
  
  ###
  ### Load in list of phenotype to plot for each gene
  ###
  phenos_to_plot_by_gene <- read.table(map_gene_pheno_fp, header=T, sep = "\t", stringsAsFactors=F)
  diseases_of_interest <- unique(phenos_to_plot_by_gene$Disease)
  
  ###
  ### Reduce LOF carrier list to only variants of interest
  ###
  lof_reduce_to_keep <- lof[which(lof$Variant %in% lof_to_keep$Variant),]
  
  ###
  ### Reduce ClinVar carrier list to only variants of interest
  ###
  clinvar_reduce_to_keep <- clinvar[which(clinvar$Variant %in% clinvar_to_keep$Variant),]
  #lof_split = split(lof_reduce_to_keep, lof_reduce_to_keep$V4)
  
  ###
  ### Reduce ClinVar carrier list to only variants of interest
  ###
  review_reduce_to_keep <- review[which(review$Variant %in% review_to_keep$Variant),]

  ###
  ### Merge lof carrier list with the Gene name
  ###
  lof_reduce_to_keep_w_gene <- merge(lof_reduce_to_keep, lof_to_keep,
                                     by.x="Variant",
                                     by.y="Variant")
  ###
  ### Merge ClinVar carrier list with the Gene name
  ###
  clinvar_reduce_to_keep_w_gene <- merge(clinvar_reduce_to_keep, clinvar_to_keep,
                                         by.x="Variant",
                                         by.y="Variant")
  ###
  ### Merge ClinVar carrier list with the Gene name
  ###
  review_reduce_to_keep_w_gene <- merge(review_reduce_to_keep, review_to_keep,
                                        by.x="Variant",
                                        by.y="Variant")

  lof_reduce_to_keep_w_gene_and_pheno <- merge(lof_reduce_to_keep_w_gene,
                                               AD_gene_list %>% filter(LoF == 1),
                                               by.x="Gene",
                                               by.y="Gene_name")
  
  review_reduce_to_keep_w_gene_and_pheno <- merge(review_reduce_to_keep_w_gene,
                                                  AD_gene_list,
                                                  by.x="Gene",
                                                  by.y="Gene_name")
  
  clinvar_reduce_to_keep_w_gene_and_pheno <- clinvar_reduce_to_keep_w_gene %>% separate_rows(Traits, sep=", ")
  clinvar_reduce_to_keep_w_gene_and_pheno$Trait = clinvar_reduce_to_keep_w_gene_and_pheno$Traits
  
  if(data_source == "52k"){
    cols_to_select <- c("SubjectID", "Gene", "Level", "Variant", "ClinVarOrLOF", "Trait", "Notes")
  }else{
    cols_to_select <- c("SubjectID", "Gene", "Level", "Variant", "ClinVarOrLOF", "Trait", "Variant_b37", "Notes")
  }
  
  ###
  ### Reduce the clinvar carrier list to only the columns we need 
  ###
  clinvar_reduce_to_keep_w_gene_and_pheno %>% 
    select(c(cols_to_select,"VariationID","AlleleID","LastEvaluated","high_confidence_lab_2017_clinical_significance")) -> carrier_list_clinvar_clean
  
  ###
  ### Reduce the lof carrier list to only the columns we need 
  ###
  lof_reduce_to_keep_w_gene_and_pheno %>% 
    select(cols_to_select) %>%
    left_join(carrier_list_clinvar_clean %>% 
                select("Variant","Trait","VariationID","AlleleID","LastEvaluated","high_confidence_lab_2017_clinical_significance") %>%
                unique(), by=c("Variant"="Variant","Trait"="Trait")) -> carrier_list_lof_clean
  
  ###
  ### Reduce the review carrier list to only the columns we need 
  ###
  review_reduce_to_keep_w_gene_and_pheno %>% 
    select(cols_to_select) %>%
    left_join(carrier_list_clinvar_clean %>% 
                select("Variant","Trait","VariationID","AlleleID","LastEvaluated","high_confidence_lab_2017_clinical_significance") %>%
                unique(), by=c("Variant"="Variant","Trait"="Trait"))-> carrier_list_review_clean
  
  ###
  ### combine the clinvar carrier list and the lof carrier list
  ###
  full_carrier_list <- rbind(carrier_list_clinvar_clean, carrier_list_lof_clean, carrier_list_review_clean)
  if(data_source == "52k"){
    full_carrier_list$Variant_b37 = full_carrier_list$Variant
  }
  full_carrier_list <- full_carrier_list %>%
    filter(SubjectID %in% pheno$combined_id) %>%
    left_join(gnomad_AF, by=c("Variant" = "Variant"))

  full_carrier_list %>% 
    select("SubjectID","ClinVarOrLOF","Variant","Gene","Trait","Level", "gnomad_exomes_AF", "gnomad_exomes_AC") %>%
    unique() %>%
    group_by(SubjectID,Variant,Trait,Gene) %>% 
    mutate(ClinVarOrLOF = paste(sort(unique(ClinVarOrLOF)),collapse=", ")) %>%
    ungroup() %>% 
    mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="ClinVar, LOF", "ClinVar and LOF")) %>%
    mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="LOF, Review", "LOF and Review")) %>%
    unique() -> all_carriers_collapse
  
  full_carrier_list %>%
    select("SubjectID","ClinVarOrLOF","Variant","Variant_b37","Gene","Trait","Level","gnomad_exomes_AF","gnomad_exomes_AC") %>%
    unique() %>%
    group_by(ClinVarOrLOF,Variant,Gene,Trait,Level,gnomad_exomes_AF,gnomad_exomes_AC) %>%
    mutate(Carrier_count = dplyr::n()) %>%
    select("ClinVarOrLOF","Variant","Variant_b37","Gene","Trait","Level","Carrier_count","gnomad_exomes_AF","gnomad_exomes_AC") %>%
    unique() %>%
    ungroup()-> full_variant_list
  
  
  write.table(full_variant_list,paste("../output/PostCurationAnalysis", data_source,"full_variant_list.txt",sep = "/"),quote=F,sep="\t",col.names=T,row.names=F)
  
  ###
  ### List of all pathogenic carriers, remove trait info and annotate if in both clinvar and lof
  ###
  full_carrier_list %>% 
    filter(tolower(Level) %in% path_levels) %>%
    select("SubjectID","ClinVarOrLOF","Variant","Gene","Trait","gnomad_exomes_AF","gnomad_exomes_AC") %>%
    unique() %>%
    group_by(SubjectID,Variant,Trait,Gene) %>% 
    mutate(ClinVarOrLOF = paste(sort(unique(ClinVarOrLOF)),collapse=", ")) %>%
    ungroup() %>% 
    mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="ClinVar, LOF", "ClinVar and LOF")) %>%
    mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="LOF, Review", "LOF and Review")) %>%
    unique() %>%
    ungroup() -> pathogenic_carriers_collapse
  
  ###
  ### List of all carriers, remove trait info and annotate if in both clinvar and lof
  ###
  full_carrier_list %>% 
    select("SubjectID","ClinVarOrLOF","Variant","Gene","Trait","gnomad_exomes_AF","gnomad_exomes_AC") %>%
    unique() %>%
    group_by(SubjectID,Variant,Trait,Gene) %>% 
    mutate(ClinVarOrLOF = paste(sort(unique(ClinVarOrLOF)),collapse=", ")) %>%
    ungroup() %>% 
    mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="ClinVar, LOF", "ClinVar and LOF")) %>%
    mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="LOF, Review", "LOF and Review")) %>%
    unique() %>%
    ungroup() -> all_carriers_collapse
  
  
  
  ###
  ### Add the phenotype info to the carrier list
  ### BEWARE: This will have multiple entries per individual if they are 
  ### carriers of multiple variants or a variant is both clinvar and lof
  ### 
  pheno_annotate_carrier <- merge(pheno, full_carrier_list, by.x="combined_id", by.y="SubjectID")
  
  ###
  ### Keep only carriers of pathogenic variants
  ### 
  pheno_annotate_carrier_more_severe <- subset(pheno_annotate_carrier, tolower(Level) %in% path_levels)
  
  
  write.table(pheno_annotate_carrier,paste("../output/PostCurationAnalysis", data_source,"phenotype_file_all_carriers_with_gene_and_trait.txt",sep = "/"),quote=F,sep="\t",col.names=T,row.names=F)
  write.table(pheno_annotate_carrier_more_severe,paste("../output/PostCurationAnalysis", data_source,"phenotype_file_all_carriers_with_gene_and_trait_most_confident.txt",sep = "/"),quote=F,sep="\t",col.names=T,row.names=F)

  ###
  ### Load the samples that have to be excluded based on PCA ancestry exclusion done by Jason
  ###
  exclude_files <- grep(".mds.exclude", dir("../data/"), value=TRUE)
  
  samples_exclude <- NULL
  for (i in exclude_files[-1]){
    to_exclude <- read_tsv(paste0("../data/", i), col_names = "combined_id")
    samples_exclude <- rbind(samples_exclude, to_exclude)
  }

  return(list("pheno" = pheno,
         "AD_gene_list" = AD_gene_list,
         "LoF_Genes" = LoF_Genes,
         "phenos_to_plot_by_gene" = phenos_to_plot_by_gene,
         "diseases_of_interest" = diseases_of_interest,
         "full_variant_list" = full_variant_list,
         "full_carrier_list" = full_carrier_list,
         "pheno_annotate_carrier" = pheno_annotate_carrier,
         "pheno_annotate_carrier_more_severe" = pheno_annotate_carrier_more_severe,
         "pathogenic_carriers_collapse" = pathogenic_carriers_collapse,
         "all_carriers_collapse" = all_carriers_collapse))
}

load_both_datasets <- function(age_cutoff, keep_ascertained=FALSE){
  my_52k_data <- load_data("52k", age_cutoff, keep_ascertained)
  my_ukbb_data <- load_data("ukbb_regeneron", age_cutoff, keep_ascertained)
  
  #Same in both datasets
  AD_gene_list <- my_52k_data$AD_gene_list
  LoF_Genes <- my_52k_data$LoF_Genes
  phenos_to_plot_by_gene <- my_52k_data$phenos_to_plot_by_gene
  diseases_of_interest <- my_52k_data$diseases_of_interest
  
  #Only AMP-T2D-GENES
  my_52k_data$pheno$data_source = "AMP-T2D-GENES"
  my_52k_data$full_variant_list$data_source = "AMP-T2D-GENES"
  my_52k_data$full_carrier_list$data_source = "AMP-T2D-GENES"
  my_52k_data$pheno_annotate_carrier$data_source = "AMP-T2D-GENES"
  my_52k_data$pheno_annotate_carrier_more_severe$data_source = "AMP-T2D-GENES"
  my_52k_data$pathogenic_carriers_collapse$data_source = "AMP-T2D-GENES"
  my_52k_data$all_carriers_collapse$data_source = "AMP-T2D-GENES"
  my_52k_pheno <- my_52k_data$pheno
  my_52k_full_variant_list <- my_52k_data$full_variant_list
  my_52k_full_carrier_list <- my_52k_data$full_carrier_list
  my_52k_pheno_annotate_carrier <- my_52k_data$pheno_annotate_carrier
  my_52k_pheno_annotate_carrier_more_severe <- my_52k_data$pheno_annotate_carrier_more_severe
  my_52k_pathogenic_carriers_collapse <- my_52k_data$pathogenic_carriers_collapse
  my_52k_all_carriers_collapse <- my_52k_data$all_carriers_collapse

  #Only UKBB
  my_ukbb_pheno <- my_ukbb_data$pheno
  my_ukbb_full_variant_list <- my_ukbb_data$full_variant_list
  my_ukbb_full_carrier_list <- my_ukbb_data$full_carrier_list
  my_ukbb_pheno_annotate_carrier <- my_ukbb_data$pheno_annotate_carrier
  my_ukbb_pheno_annotate_carrier_more_severe <- my_ukbb_data$pheno_annotate_carrier_more_severe
  my_ukbb_pathogenic_carriers_collapse <- my_ukbb_data$pathogenic_carriers_collapse
  my_ukbb_all_carriers_collapse <- my_ukbb_data$all_carriers_collapse
  my_ukbb_pheno$combined_id = as.character(my_ukbb_pheno$combined_id)
  my_ukbb_pheno$data_source = "UKBB"
  my_ukbb_full_variant_list$data_source = "UKBB"
  my_ukbb_full_carrier_list$data_source = "UKBB"
  my_ukbb_full_carrier_list$SubjectID = as.character(my_ukbb_full_carrier_list$SubjectID)
  my_ukbb_pheno_annotate_carrier$data_source = "UKBB"
  my_ukbb_pheno_annotate_carrier$combined_id = as.character(my_ukbb_pheno_annotate_carrier$combined_id)
  my_ukbb_pheno_annotate_carrier_more_severe$data_source = "UKBB"
  my_ukbb_pheno_annotate_carrier_more_severe$combined_id = as.character(my_ukbb_pheno_annotate_carrier_more_severe$combined_id)
  my_ukbb_pathogenic_carriers_collapse$data_source = "UKBB"
  my_ukbb_pathogenic_carriers_collapse$SubjectID = as.character(my_ukbb_pathogenic_carriers_collapse$SubjectID)
  my_ukbb_all_carriers_collapse$data_source = "UKBB"
  my_ukbb_all_carriers_collapse$SubjectID = as.character(my_ukbb_all_carriers_collapse$SubjectID)
  
  
  #Merge datasets
  joint_data <- NULL
  joint_data$pheno = full_join(my_52k_pheno, my_ukbb_pheno)
  joint_data$full_variant_list = full_join(my_52k_full_variant_list, my_ukbb_full_variant_list)
  joint_data$full_carrier_list = full_join(my_52k_full_carrier_list, my_ukbb_full_carrier_list)
  joint_data$pheno_annotate_carrier = full_join(my_52k_pheno_annotate_carrier, my_ukbb_pheno_annotate_carrier)
  joint_data$pheno_annotate_carrier_more_severe = full_join(my_52k_pheno_annotate_carrier_more_severe, my_ukbb_pheno_annotate_carrier_more_severe)
  joint_data$pathogenic_carriers_collapse = full_join(my_52k_pathogenic_carriers_collapse, my_ukbb_pathogenic_carriers_collapse)
  joint_data$all_carriers_collapse = full_join(my_52k_all_carriers_collapse, my_ukbb_all_carriers_collapse)
  return(list("AD_gene_list" = AD_gene_list,
              "LoF_Genes" = LoF_Genes,
              "phenos_to_plot_by_gene" = phenos_to_plot_by_gene,
              "diseases_of_interest" =  diseases_of_interest,
              "data_52k" = my_52k_data,
              "data_ukbb" = my_ukbb_data,
              "joint_data" = joint_data))
}

annotate_phenotype_with_carriers <- function(pheno,
                                             gene_list,
                                             carriers_list,
                                             trait,
                                             phenotypes_to_keep,
                                             levels_label_as_carrier,
                                             collapse_carrier_multiple = FALSE,
                                             remove_vus = FALSE,
                                             stratify_vus = FALSE,
                                             by_gene = FALSE){
  ###
  ### Function to provide a table containing phenotype and carrier info 
  ###   pheno - full phenotype file to be annotated with carrier info
  ###   gene_list - list of genes to include in carrier annotation
  ###   carrier_list - list of all carriers with gene, trait, variant, and level info
  ###   trait - trait to subset carrier list by
  ###   phenotypes_to_keep - any phenotype to keep in the table, should include the phenotype of interest
  ###   levels_label_as_carrier - what carrier levels should be labeled as a carrier
  ###   collapse_carrier_multiple - individuals can have multiple variants making them a carrier
  ###        or a variant can be both clinvar and LOF, do you want that info collapsed so all you know 
  ###        is they are a carrier not the specific variant or the specific level
  ###   remove_vus - should VUS carriers be removed, default is false
  ###   stratify_vus - should VUS carriers be an additional category, default is false,
  ###        will only be used if remove_vus is false
  ###
  
  if(remove_vus & stratify_vus){
    stop("Can't remove VUS and stratify by VUS")
  }
  
  #filter carrier list to only include gene and trait of interest
  carriers_list %>% 
    filter(Gene %in% gene_list & Trait == trait) %>%
    select(-data_source)-> carriers_list_gene_trait
  
  #If we are going to remove the VUS non carriers, get a list of them
  if(remove_vus | stratify_vus){
    carriers_list_gene_trait %>%
      filter(tolower(Level) == "vus") -> carriers_list_gene_trait_vus
  }
  
  #Now filter carrier list to only the levels we want to keep
  carriers_list_gene_trait %>%
    filter(tolower(Level) %in% levels_label_as_carrier) -> carriers_list_gene_trait_level 
  
  #Keep combined_id if not already included
  phenotypes_to_keep <- unique(c("combined_id", phenotypes_to_keep, "data_source"))
  
  #Only keep phenotypes requested, add a Carrier column, then join with other carrier info
  pheno %>%
    select(phenotypes_to_keep) %>% 
    mutate(Carrier = ifelse(combined_id %in% carriers_list_gene_trait_level$SubjectID,"Carrier","Non Carrier")) %>%
    left_join(carriers_list_gene_trait,by = c("combined_id" = "SubjectID")) -> pheno_carrier
  
  #If requested, remove the VUS non carriers
  if(remove_vus){
    pheno_carrier %>%
      mutate(Carrier = replace(Carrier,(combined_id %in% carriers_list_gene_trait_vus$SubjectID) & (Carrier == "Non Carrier"),NA)) -> pheno_carrier
  }
  
  if(stratify_vus){
    pheno_carrier %>%
      mutate(Carrier = replace(Carrier,(combined_id %in% carriers_list_gene_trait_vus$SubjectID) & (Carrier == "Non Carrier"),"VUS")) -> pheno_carrier
  }
  
  #If a variant is in both clinvar and lof then a carrier will be represented twice
  #Also, if an individual has multiple variants making them a carrier they will be represented twice
  #If requested, this will collapse those and then relabel the ClinVarOrLOF field if appropriate
  #If a person has multiple pathogenic variants then they are labeled a carrier and the retained gnomAD AF will be the lowest AF
  #If a person has a pathogenic variant and a VUS they are a carrier and the retained gnomAD AF will be that of the pathogenic variant
  #If a person has multiple vus variants they will be labeled VUS and the retained gnomAD AF will be the lowest AF
  if(collapse_carrier_multiple){
    pheno_carrier <- pheno_carrier %>%
      group_by(combined_id,Variant,data_source,.drop = FALSE) %>% 
      mutate(ClinVarOrLOF = paste(sort(unique(ClinVarOrLOF)),collapse=", ")) %>% 
      ungroup() %>% 
      mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="ClinVar, LOF", "ClinVar and LOF")) %>% 
      mutate(ClinVarOrLOF=replace(ClinVarOrLOF, ClinVarOrLOF=="LOF, Review", "LOF and Review")) %>% 
      group_by(combined_id,data_source,.drop = FALSE) %>% 
      mutate(Carrier_collapse = paste(sort(unique(Carrier)),collapse=",")) %>%
      filter((Carrier == "Carrier")  | (Carrier_collapse != "Carrier,VUS")) %>% 
      ungroup() %>% 
      group_by(combined_id,data_source,Carrier,.drop = FALSE) %>% 
      slice(ifelse(length(which.min(gnomad_exomes_AF)),which.min(gnomad_exomes_AF),1)) %>%
      ungroup() 
    if(by_gene){
      pheno_carrier <- pheno_carrier %>%
        mutate(Gene = ifelse((Carrier == "Non Carrier"),"Non Carrier",Gene)) %>%
        select(-Variant,-Level) %>%
        unique()
    }else{
      pheno_carrier <- pheno_carrier %>%
        select(-Variant,-Level,-Gene) %>%
        unique()
    }
  }
  
  if(stratify_vus){
    pheno_carrier <- pheno_carrier %>%
      mutate(AF_group = case_when(gnomad_exomes_AC == 0 ~ "Novel",
                                  is.na(gnomad_exomes_AC) ~ "Novel",
                                  gnomad_exomes_AF < (1/10000) ~ "<1/10000",
                                  gnomad_exomes_AF < (1/1000) ~ "<1/1000",
                                  gnomad_exomes_AF >= (1/1000) ~ ">=1/1000")) %>%
      mutate(Carrier_simple = Carrier,
             Carrier = ifelse(Carrier == "VUS", paste("VUS",AF_group,sep = " "), Carrier))
  }
  pheno_carrier$Trait = trait
  
  return(pheno_carrier)
}


# These pvalue files were created by combining the output from multiple EPACTS runs where 
# group files of variants were made for each combination displayed in Supplementary Tables 4 
# and 5 and other details are described in the paper Methods.
load_epacts_results <- function(by_gene=FALSE){
  if(by_gene){
    prefix <- "gene"
  }else{
    prefix <- "combined"
  }
  p_values_52k <- read_tsv(paste0("../output/EpactsResults/ALL_", prefix, "_52k_epacts_results_w_adj.txt"))
  p_values_52k$data_source = "AMP-T2D-GENES"
  
  p_values_ukbb <- read_tsv(paste0("../output/EpactsResults/ALL_", prefix, "_ukbb_regeneron_epacts_results_w_adj.txt"))
  p_values_ukbb$data_source = "UKBB"
  
  p_values <- rbind(p_values_52k, p_values_ukbb)
  return(p_values)
}

ukbb_filter_carrier_list_AB_DP_GQ <- function(carrier_list, carrier_list_qual_info){
  carrier_list_qual_info$Variant = paste0("chr", carrier_list_qual_info$V2)
  if(nrow(carrier_list) != nrow(carrier_list_qual_info)){print("WARNING: unequal rows")}
  carrier_list_join_qual <- carrier_list %>%
    left_join(carrier_list_qual_info, by = c("Variant"="Variant","SubjectID"="V1")) %>%
    filter(V7=="False")
    return(carrier_list_join_qual)
}

carrier_list <- read.delim("../output/GetCarrierLists/ukbb_regeneron/clinvar_carriers_before_AB_GQ_DP_filter.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
carrier_list_qual_info <- read.delim("../output/GetCarrierLists/ukbb_regeneron/ukbb_regeneron_clinvar_carriers_tabix_results_clean.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
carrier_list_join_qual <- ukbb_filter_carrier_list_AB_DP_GQ(carrier_list, carrier_list_qual_info)
carrier_list_join_qual %>% write_delim("../output/GetCarrierLists/ukbb_regeneron/clinvar_carrier.txt",delim="\t")

carrier_list <- read.delim("../output/GetCarrierLists/ukbb_regeneron/lof_carriers_before_AB_GQ_DP_filter.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
carrier_list_qual_info <- read.delim("../output/GetCarrierLists/ukbb_regeneron/ukbb_regeneron_lof_carriers_tabix_results_clean.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
carrier_list_join_qual <- ukbb_filter_carrier_list_AB_DP_GQ(carrier_list, carrier_list_qual_info)
carrier_list_join_qual %>% write_delim("../output/GetCarrierLists/ukbb_regeneron/lof_carrier.txt",delim="\t")

carrier_list <- read.delim("../output/GetCarrierLists/ukbb_regeneron/radiant_carriers_before_AB_GQ_DP_filter.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
carrier_list_qual_info <- read.delim("../output/GetCarrierLists/ukbb_regeneron/ukbb_regeneron_radiant_carriers_tabix_results_clean.txt", sep = "\t", header=FALSE, stringsAsFactors = FALSE)
carrier_list_join_qual <- ukbb_filter_carrier_list_AB_DP_GQ(carrier_list, carrier_list_qual_info)
carrier_list_join_qual %>% write_delim("../output/GetCarrierLists/ukbb_regeneron/radiant_carrier.txt",delim="\t")
