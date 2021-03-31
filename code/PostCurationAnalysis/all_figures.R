library(ggplot2)
library(Hmisc)
library(scales)
library(ggbeeswarm)
library(ggExtra)
library(gridExtra)
library(ggforce)
library(cowplot)
library(stringr)
library(janitor)
library(tidyr)
library(plyr)
library(reshape2)
library(ggpubr)
library(readr)
library(dplyr)

age_cutoff <- 40


###
### Load all data
###

source("../code/PostCurationAnalysis/Load_dataBoth.R")
source("../code/PostCurationAnalysis/plotting.R")
all_data <- load_both_datasets(age_cutoff)
joint_data <- all_data$joint_data
AD_gene_list <- all_data$AD_gene_list
p_values <- load_epacts_results()
p_values <- p_values %>% filter(subset == "All")
p_values_gene <- load_epacts_results(by_gene=TRUE)
p_values_gene <- p_values_gene %>% filter(subset == "All")

all_penetrance <- penetrance_table(joint_data)
all_penetrance_non_path_no_vus <- penetrance_table(joint_data, levels_to_plot=non_path_no_vus_levels)

all_penetrance_gene <- penetrance_table(joint_data, by_gene=TRUE)
all_penetrance_gene_non_path_no_vus <- penetrance_table(joint_data, by_gene=TRUE, levels_to_plot=non_path_no_vus_levels)

joint_data_adj <- joint_data
joint_data_adj$pheno = joint_data_adj$pheno %>% mutate(ldl_mgdl = ldl_mgdl_adj, ldl_mgdl = ldl_mgdl_adj, tg_mgdl = tg_mgdl_adj)
all_penetrance_adj <- penetrance_table(joint_data_adj)
all_penetrance_gene_adj <- penetrance_table(joint_data_adj, by_gene=TRUE)


###
### Figure 1. Curation of ClinVar and pLoF variants across the monogenic conditions.
###

p <- plot_num_variants_by_condition_all_stacked(joint_data)

pdf('../output/PostCurationAnalysis/Figure1.pdf',height=5,width=12)
print(p)
dev.off()


###
### Figure 2. Carriers of rare clinically significant monogenic variants for lipid conditions and monogenic diabetes have more extreme effect size estimates than individuals with the top 1% of global extended polygenic scores (gePS). Lipids adjusted for medication use
###

p <- PRS_plots(joint_data_adj, age_cutoff = 60)
pdf("../output/PostCurationAnalysis/Figure2.pdf", width = 10, height=14)
print(p$p1)
dev.off()


percentile_table <- p$percentile_table %>%
select(Condition, percentile,N,phenotype_mean = phenotype_value,CI=ci) %>%
mutate(percentile = ifelse(percentile=="C","Carrier",percentile))

# Figure2 in table form
write.table(percentile_table,"../output/PostCurationAnalysis/Figure2_table.txt",sep = "\t",quote=F,row.names = F,col.names = T)


###
### Figure 3. Phenotype distributions and penetrance estimates of clinically significant variant carriers.
###

all_penetrance_main <- all_penetrance_adj %>%
filter((condition_simple %in% c("High LDL",
"Low LDL",
"High HDL",
"High TG",
"Monogenic obesity",
"MODY")) & (!lipidmeds)) %>%
mutate(order = condition)

all_penetrance_gene_main <- all_penetrance_gene_adj %>%
filter(Gene %in% c("GCK","HNF1A","APOB","PCSK9","LDLR","MC4R","APOA5","LPL","CETP")) %>%
filter((condition_simple %in% c("High LDL",
"Low LDL",
"High TG",
"MODY")) & (!lipidmeds)) %>%
mutate(order = paste(condition,Gene,sep = " ")) %>%
mutate(condition = Gene) %>%
select(-Gene)

all_penetrance_main <- rbind(all_penetrance_main, all_penetrance_gene_main)

all_penetrance_main$order = factor(all_penetrance_main$order,
levels = c("High LDL (>190 mg/dl)",
"High LDL (>190 mg/dl) APOB",
"High LDL (>190 mg/dl) LDLR",
"Low LDL (<80 mg/dl)",
"Low LDL (<80 mg/dl) APOB",
"Low LDL (<80 mg/dl) PCSK9",
"High HDL (>70 mg/dl)",
"High triglycerides (>200 mg/dl)",
"High triglycerides (>200 mg/dl) APOA5",
"High triglycerides (>200 mg/dl) LPL",
"Monogenic obesity (BMI >30)",
"MODY (T2D)",
"MODY (T2D) GCK",
"MODY (T2D) HNF1A",
"MODY (T2D and pre-diabetes)",
"MODY (T2D and pre-diabetes) GCK",
"MODY (T2D and pre-diabetes) HNF1A"))

p_value_adj <- p_values
p_value_adj <- p_value_adj %>% filter(!((grepl("ldl_", Trait) | grepl("tg_", Trait)) & !grepl("_adj", Trait))) %>% mutate(Trait = str_replace(Trait, '_adj', ''))
p1 <- plot_carrier_phenotypes_by_condition(joint_data_adj, AD_gene_list, p_value_adj %>% filter(type=='pathogenic'))
all_penetrance_main <- all_penetrance_main %>% filter(!is.na(order))
p2 <- penetrance_plot(all_penetrance_main, labels = rev(c("High LDL (>190 mg/dl)",
                                                          "APOB",
                                                          "LDLR",
                                                          "Low LDL (<80 mg/dl)",
                                                          "APOB",
                                                          "PCSK9",
                                                          "High HDL - CETP (>70 mg/dl)",
                                                          "High triglycerides (>200 mg/dl)",
                                                          "APOA5",
                                                          "LPL",
                                                          "Monogenic obesity - MC4R (BMI >30)",
                                                          "MODY (T2D)",
                                                          "GCK",
                                                          "HNF1A",
                                                          "MODY (T2D and pre-diabetes)",
                                                          "GCK",
                                                          "HNF1A")))
p2 <- p2 +
theme(axis.text.y=element_text(face=rev(c("plain","italic","italic","plain","italic","italic","plain","plain","italic","italic","plain","plain","italic","italic","plain","italic","italic"))))

pdf("../output/PostCurationAnalysis/Figure3.pdf",height=15,width=12)
print(ggarrange(p1, p2, heights = c(1,1.5), nrow=2,labels = c("A","B")))
dev.off()


###
### Figure 4. Ascertainment bias significantly impacts expressivity of clinically significant variants for LDL cholesterol conditions.
###

all_data_ascertained <- load_data("52k", age_cutoff, keep_ascertained=TRUE)
ascertainment_data <- ascertainment_plots(all_data_ascertained)

pdf("../output/PostCurationAnalysis/Figure4.pdf", width=13, height=8)
ascertainment_data$ascertainment_plot
dev.off()


###
### Figure 5. The combination of clinically significant monogenic variants and corresponding polygenic scores significantly improves prediction for high HDL cholesterol and high triglyceride conditions.
###

p <- PRS_plots(joint_data_adj, age_cutoff = 60)

pdf("../output/PostCurationAnalysis/Figure5.pdf", width = 8, height=9)
print(p$p2)
dev.off()


###
### Supplementary Figure 1. Distribution of clinically significant variants across ancestries.
###

p <- plot_carrier_ancestry_by_condition_proportion(all_data$data_52k)

pdf("../output/PostCurationAnalysis/Supplementary_Figure1.pdf", width = 10, height=4.8)
print(p)
dev.off()



###
### Supplementary Figure 2. Carriers of clinically significant variants in MODY genes show a younger age of diabetes diagnosis compared to the rest of the AMP-T2D-GENES cohorts and UK Biobank population.
###

p <- age_of_diagnosis_survival_plot(joint_data, AD_gene_list)

pdf("../output/PostCurationAnalysis/Supplementary_Figure2.pdf",height=11,width=8.5)
print(p)
dev.off()




###
### Getting data for Supplementary Table 1.
###

participant_characteristic_table <- participant_characteristics(joint_data$pheno)
write.table(participant_characteristic_table$all, "../output/PostCurationAnalysis/Supplementary_Table1_data.txt", sep = "\t",quote = FALSE, row.names = FALSE)


###
### Getting data for Supplementary Table 3.
###

all_trait_carrier_count_table <- carrier_count_table_by_condition(joint_data$pheno, joint_data$full_carrier_list %>% filter(tolower(Level) %in% path_levels))

all_trait_carrier_count_table <- all_trait_carrier_count_table %>%
filter(data_source == "AMP-T2D-GENES") %>%
full_join(all_trait_carrier_count_table %>% filter(data_source == "UKBB"), by=c("Trait"="Trait")) %>%
ungroup()

all_trait_carrier_count_table <- all_trait_carrier_count_table %>%
select("Condition" = Trait,
"AMP-T2D-GENES path variants" = number_variants.x,
"AMP-T2D-GENES carriers" = number_carriers.x,
"AMP-T2D-GENES percent carriers" = percent_carriers.x,
"UKBB path variants" = number_variants.y,
"UKBB carriers" = number_carriers.y,
"UKBB percent carriers" = percent_carriers.y)

write_delim(all_trait_carrier_count_table, "../output/PostCurationAnalysis/Supplementary_Table3_data.txt", delim = "\t")



###
### Getting data for Supplementary Table 4.
###

all_penetrance_format <- formatted_penetrance(all_penetrance, p_values)
pheno_mean_ci <- carrier_phenotypes_mean_ci_pvalues(joint_data, AD_gene_list, conditions_cont_phenotypes, p_values)
all_penetrance_format_pheno <- all_penetrance_format %>%
left_join(pheno_mean_ci, by=c("Condition"="Condition"), suffix = c("_binary", "_cont")) %>%
select("Condition","Restricted no lipid medications",
"total_Non Carrier_AMP-T2D-GENES","Mean (95% CI) Non Carrier_AMP-T2D-GENES","count_Non Carrier_AMP-T2D-GENES","prop_Non Carrier_AMP-T2D-GENES",
"total_Carrier_AMP-T2D-GENES","Mean (95% CI) Carrier_AMP-T2D-GENES","count_Carrier_AMP-T2D-GENES","prop_Carrier_AMP-T2D-GENES",
"Beta (se)_AMP-T2D-GENES_cont","P value_AMP-T2D-GENES_cont","Beta (se)_AMP-T2D-GENES_binary","P value_AMP-T2D-GENES_binary","OR_AMP-T2D-GENES",
"total_Non Carrier_UKBB","Mean (95% CI) Non Carrier_UKBB","count_Non Carrier_UKBB","prop_Non Carrier_UKBB",
"total_Carrier_UKBB","Mean (95% CI) Carrier_UKBB","count_Carrier_UKBB","prop_Carrier_UKBB",
"Beta (se)_UKBB_cont","P value_UKBB_cont","Beta (se)_UKBB_binary","P value_UKBB_binary","OR_UKBB")
write.table(all_penetrance_format_pheno,"../output/PostCurationAnalysis/Supplementary_Table4_data_combined.txt",sep = "\t",quote=F,row.names = F,col.names = T)

all_penetrance_format <- formatted_penetrance(all_penetrance_gene, p_values_gene, by_gene=TRUE)
pheno_mean_ci_gene <- carrier_phenotypes_mean_ci_pvalues(joint_data, AD_gene_list, conditions_cont_phenotypes, p_values_gene, by_gene = TRUE)
all_penetrance_format_pheno_gene <- all_penetrance_format %>%
left_join(pheno_mean_ci_gene, by=c("Condition"="Condition","Gene"="Gene"), suffix = c("_binary", "_cont")) %>%
select("Condition","Gene","total_AMP-T2D-GENES", "Mean (95% CI)_AMP-T2D-GENES", "prop_AMP-T2D-GENES","Beta (se)_AMP-T2D-GENES_cont","P value_AMP-T2D-GENES_cont","OR_AMP-T2D-GENES",
"P value_AMP-T2D-GENES_binary", "total_UKBB","Mean (95% CI)_UKBB", "prop_UKBB","Beta (se)_UKBB_cont","P value_UKBB_cont","OR_UKBB", "P value_UKBB_binary")
write.table(all_penetrance_format_pheno_gene,"../output/PostCurationAnalysis/Supplementary_Table4_data_gene.txt",sep = "\t",quote=F,row.names = F,col.names = T)





###
### Getting data for Supplementary Table 5.
###

all_penetrance_format_non_path_no_vus <- formatted_penetrance(all_penetrance_non_path_no_vus, p_values %>% filter(type=='non_pathogenic_no_vus'), carrier_type = "non_pathogenic_no_vus")
pheno_mean_ci_non_path_no_vus <- carrier_phenotypes_mean_ci_pvalues(joint_data, AD_gene_list, conditions_cont_phenotypes, p_values %>% filter(type=='non_pathogenic_no_vus'), carrier_type='non_pathogenic_no_vus', levels_to_use=non_path_no_vus_levels)
all_penetrance_format_pheno_non_path_no_vus <- all_penetrance_format_non_path_no_vus %>%
  left_join(pheno_mean_ci_non_path_no_vus, by=c("Condition"="Condition"), suffix = c("_binary", "_cont")) %>%
  select("Condition","Restricted no lipid medications",
         "total_Non Carrier_AMP-T2D-GENES","Mean (95% CI) Non Carrier_AMP-T2D-GENES","count_Non Carrier_AMP-T2D-GENES","prop_Non Carrier_AMP-T2D-GENES",
         "total_Carrier_AMP-T2D-GENES","Mean (95% CI) Carrier_AMP-T2D-GENES","count_Carrier_AMP-T2D-GENES","prop_Carrier_AMP-T2D-GENES",
         "Beta (se)_AMP-T2D-GENES_cont","P value_AMP-T2D-GENES_cont","Beta (se)_AMP-T2D-GENES_binary","P value_AMP-T2D-GENES_binary","OR_AMP-T2D-GENES",
         "total_Non Carrier_UKBB","Mean (95% CI) Non Carrier_UKBB","count_Non Carrier_UKBB","prop_Non Carrier_UKBB",
         "total_Carrier_UKBB","Mean (95% CI) Carrier_UKBB","count_Carrier_UKBB","prop_Carrier_UKBB",
         "Beta (se)_UKBB_cont","P value_UKBB_cont","Beta (se)_UKBB_binary","P value_UKBB_binary","OR_UKBB")
write.table(all_penetrance_format_pheno_non_path_no_vus,"../output/PostCurationAnalysis/Supplementary_Table5_data_combined.txt",sep = "\t",quote=F,row.names = F,col.names = T)

all_penetrance_format_non_path_no_vus <- formatted_penetrance(all_penetrance_gene_non_path_no_vus, p_values_gene %>% filter(type=='non_pathogenic_no_vus'), by_gene=TRUE, carrier_type = "non_pathogenic_no_vus")
pheno_mean_ci_gene_non_path_no_vus <- carrier_phenotypes_mean_ci_pvalues(joint_data, AD_gene_list, conditions_cont_phenotypes, p_values_gene %>% filter(type=='non_pathogenic_no_vus'), by_gene = TRUE, carrier_type='non_pathogenic_no_vus', levels_to_use=non_path_no_vus_levels)
all_penetrance_format_pheno_gene_non_path_no_vus <- all_penetrance_format_non_path_no_vus %>%
  left_join(pheno_mean_ci_gene_non_path_no_vus, by=c("Condition"="Condition","Gene"="Gene"), suffix = c("_binary", "_cont")) %>%
  select("Condition","Gene","total_AMP-T2D-GENES", "Mean (95% CI)_AMP-T2D-GENES", "prop_AMP-T2D-GENES","Beta (se)_AMP-T2D-GENES_cont","P value_AMP-T2D-GENES_cont","OR_AMP-T2D-GENES",
         "P value_AMP-T2D-GENES_binary", "total_UKBB","Mean (95% CI)_UKBB", "prop_UKBB","Beta (se)_UKBB_cont","P value_UKBB_cont","OR_UKBB", "P value_UKBB_binary")
write.table(all_penetrance_format_pheno_gene_non_path_no_vus,"../output/PostCurationAnalysis/Supplementary_Table5_data_genes.txt",sep = "\t",quote=F,row.names = F,col.names = T)


###
### Getting data for Supplementary Table 6. Comparison of Top 1% gePS with interquartile range and monogenic carriers.
###

p <- PRS_plots(joint_data_adj, age_cutoff = 60)
write.table(p$glm_interquartile_1_per_pvalues_all, "../output/PostCurationAnalysis/Supplementary_Table6_data.txt", sep= "\t", quote=F, row.names=F,col.names = T)



###
### Getting data for Supplementary Table 7. Mean serum LDL values based on ascertainment.
###

write.table(ascertainment_data$ascertainment_stats, "../output/PostCurationAnalysis/Supplementary_Table7_data.txt", sep= "\t", quote=F, row.names=F,col.names = T)


###
### Getting data for Supplementary Table 9. Impact of polygenic score on trait expressivity in monogenic carriers.
###
# Table 9A
p <- PRS_plots(joint_data_adj, age_cutoff = 60)
write.table(p$glm_results, "../output/PostCurationAnalysis/Supplementary_Table9A_data.txt", sep= "\t", quote=F, row.names=F,col.names = T)


# Table 9B
p <- PRS_plots(joint_data_adj)
write.table(p$glm_results_all, "../output/PostCurationAnalysis/Supplementary_Table9B_data.txt", sep= "\t", quote=F, row.names=F,col.names = T)



# Getting data for Supplementary Table 12. Variant curation assessments and carrier counts.
variant_table <- final_variant_table(joint_data$full_carrier_list)
write.table(variant_table,"../output/PostCurationAnalysis/Supplementary_Table12_data.txt",sep = "\t",quote=F,row.names = F,col.names = T)



###
### Summary table of Numbers of pathogenic variants and carriers across conditions and Percent with condition used in text
###
all_curated_variant_table <- variant_count_table(joint_data$full_carrier_list)
all_curated_variant_table_path <- variant_count_table(joint_data$full_carrier_list %>% filter(tolower(Level) %in% path_levels))
colnames(all_curated_variant_table) <- c("data_source1", "ClinVar and Reviews1", "LOF1", "Total Curated1")

all_curated_variant_table <- cbind(all_curated_variant_table, all_curated_variant_table_path)

all_curated_variant_table <- all_curated_variant_table %>%
  select("Data set" = data_source, 
         "Number of curated ClinVar and Review variants" = `ClinVar and Reviews`,
         "Number of curated predicted LoF variants" = LOF,
         "Total number of curated variants" = `Total Curated`,
         "Number of P/LP ClinVar and Review variants" = `ClinVar and Reviews1`,
         "Number of likely LoF variants" = LOF1,
         "Total number of likely pathogenic variants" = `Total Curated1`) 

write_delim(all_curated_variant_table, "../output/PostCurationAnalysis/total_variant_counts.txt", delim = "\t")


