library(Hmisc)
library(ggplot2)
library(scales)
library(ggbeeswarm)
library(ggExtra)
library(gridExtra)
library(ggforce)
library(cowplot)
library(stringr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(plyr)
library(reshape2)
library(ggpubr)
library(gtable)
library(grid)
library(epiDisplay)
require(survival)
library(survminer)
library(dplyr)

source("../code/config.R")

plot_num_variants_by_condition_all_stacked <- function(data){
  pheno_annotate_carrier <- data$pheno_annotate_carrier
  variants <- unique(pheno_annotate_carrier[, c("ClinVarOrLOF", "Variant", "Gene", "Level", "Trait", "data_source")])
  variants$Trait[which(variants$Trait == "Obesity")] = "Monogenic obesity"
  variants <- variants %>%
    filter(!((Trait == "MODY Extended")&(Gene %in% c("GCK","HNF1A","HNF1B","PDX1","HNF4A")))) %>%
    mutate(Trait = replace(Trait, Trait == "MODY Extended", "MODY extended")) %>%
    mutate(Trait = replace(Trait, Trait == "Neonatal Diabetes", "Neonatal diabetes")) %>%
    filter(Trait %in% disease_order) %>%
    mutate(Trait = factor(Trait, rev(disease_order)),
           Excluded = ifelse(tolower(Level) %in% path_levels, "Pathogenic", "Excluded"),
           ClinVarOrLOF = ifelse(ClinVarOrLOF %in% c("Review","ClinVar"), "ClinVar and Reviews", ClinVarOrLOF)) %>%
    mutate(CombineClinVarExclude = factor(paste(Excluded, ClinVarOrLOF),
                                          levels = c("Excluded LOF",
                                                     "Pathogenic LOF",
                                                     "Excluded ClinVar and Reviews",
                                                     "Pathogenic ClinVar and Reviews"))) %>%
    group_by(Trait, CombineClinVarExclude, data_source) %>%
    dplyr::summarise(count = dplyr::n())
  

  
  p <- ggplot(data=variants, aes(x=Trait, y=count, fill=CombineClinVarExclude)) +
    geom_bar(stat="identity")+
    theme_classic() + theme(strip.background = element_blank(), 
                            strip.placement = "outside") + #,
    ylab("Number of variants") + xlab("Condition")+
    scale_fill_manual(values=my_paired_palette)+ 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.text.x = element_text(size = 14),
          legend.text=element_text(size=12),
          legend.title=element_blank()) +
    facet_wrap(~data_source,nrow=1) +
    coord_flip()
  
  return(p)
}


penetrance_table <- function(data,
                             include_mody=TRUE,
                             levels_to_plot=path_levels,
                             stratify_vus = FALSE,
                             by_gene = FALSE){
  pheno <- data$pheno
  full_carrier_list <- data$full_carrier_list
  
  all_threshold_results <- NULL
  
  if(include_mody){
    #MODY
    threshold_results <- penetrance_threshold(data, "MODY", t2d, "MODY (T2D)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
    
    threshold_results <- penetrance_threshold(data, "MODY", t2d_with_pre, "MODY (T2D and pre-diabetes)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
    
    #MODY Extended
    threshold_results <- penetrance_threshold(data, "MODY Extended", t2d, "MODY extended (T2D)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
    
    threshold_results <- penetrance_threshold(data, "MODY Extended", t2d_with_pre, "MODY extended (T2D and pre-diabetes)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
    
    #Lipodystrophy
    threshold_results <- penetrance_threshold(data, "Lipodystrophy", t2d, "Lipodystrophy (T2D)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
    
    threshold_results <- penetrance_threshold(data, "Lipodystrophy", t2d_with_pre, "Lipodystrophy (T2D and pre-diabetes)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
    
    #Neonatal diabetes
    threshold_results <- penetrance_threshold(data, "Neonatal Diabetes", t2d, "Neonatal diabetes (T2D)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
    
    threshold_results <- penetrance_threshold(data, "Neonatal Diabetes", t2d_with_pre, "Neonatal diabetes (T2D and pre-diabetes)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
    all_threshold_results <- rbind(all_threshold_results, threshold_results)
  }
  #Obesity 
  threshold_results <- penetrance_threshold(data, "Obesity", bmi, "Monogenic obesity (BMI >30)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  threshold_results$condition_simple = "Monogenic obesity"
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  ######HDL > 70
  threshold_results <- penetrance_threshold(data, "High HDL", hdl_mgdl, "High HDL (>70 mg/dl)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  ###### High LDL
  #clinical diagnosis of FH: low density lipoprotein cholesterol (LDL-C) >190 mg/dL
  threshold_results <- penetrance_threshold(data, "High LDL", ldl_mgdl, "High LDL (>190 mg/dl)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  threshold_results <- penetrance_threshold(data, "High LDL", ldl_mgdl, "High LDL (>190 mg/dl)", lipidmeds = TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  ###### High LDL ADJ
  #clinical diagnosis of FH: low density lipoprotein cholesterol (LDL-C) >190 mg/dL
  threshold_results <- penetrance_threshold(data, "High LDL", ldl_mgdl_adj, "High LDL ADJ (>190 mg/dl)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  threshold_results <- penetrance_threshold(data, "High LDL", ldl_mgdl_adj, "High LDL ADJ (>190 mg/dl)", lipidmeds = TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  ###### Low LDL
  threshold_results <- penetrance_threshold(data, "Low LDL", ldl_mgdl, "Low LDL (<80 mg/dl)", threshold_under=TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  threshold_results <- penetrance_threshold(data, "Low LDL", ldl_mgdl, "Low LDL (<80 mg/dl)", threshold_under=TRUE, lipidmeds = TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  ###### Low LDL ADJ
  threshold_results <- penetrance_threshold(data, "Low LDL", ldl_mgdl_adj, "Low LDL ADJ (<80 mg/dl)", threshold_under=TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  threshold_results <- penetrance_threshold(data, "Low LDL", ldl_mgdl_adj, "Low LDL ADJ (<80 mg/dl)", threshold_under=TRUE, lipidmeds = TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  
  ###### Triglycerides
  #rationale for 200 - it's considered "high triglycidemia"
  threshold_results <- penetrance_threshold(data, 'High triglycerides', tg_mgdl, "High triglycerides (>200 mg/dl)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  threshold_results$condition_simple = "High TG"
  all_threshold_results <- rbind(all_threshold_results, threshold_results)

  threshold_results <- penetrance_threshold(data, 'High triglycerides', tg_mgdl, "High triglycerides (>200 mg/dl)", lipidmeds = TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  threshold_results$condition_simple = "High TG"
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  
  
  ###### Triglycerides ADJ
  #rationale for 200 - it's considered "high triglycidemia"
  threshold_results <- penetrance_threshold(data, 'High triglycerides', tg_mgdl_adj, "High triglycerides ADJ (>200 mg/dl)", levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  threshold_results$condition_simple = "High TG"
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  threshold_results <- penetrance_threshold(data, 'High triglycerides', tg_mgdl_adj, "High triglycerides ADJ (>200 mg/dl)", lipidmeds = TRUE, levels_to_plot=levels_to_plot, stratify_vus = stratify_vus, by_gene = by_gene)
  threshold_results$condition_simple = "High TG"
  all_threshold_results <- rbind(all_threshold_results, threshold_results)
  
  
  if(by_gene){
    get_CI <- function(x){
      if((!is.na(x[5])) & (!is.na(x[7]) & as.numeric(x[7])>0)){
        CI <- binom.test(as.numeric(x[5]), as.numeric(x[7]), conf.level = 0.95)[[4]]
        return(list(CI[[1]],CI[[2]]))
      }
      else{
        return(list(NA,NA))
      }
    }
  }else{
    get_CI <- function(x){
      if((!is.na(x[4])) & (!is.na(x[6]) & as.numeric(x[6])>0)){
        CI <- binom.test(as.numeric(x[4]), as.numeric(x[6]), conf.level = 0.95)[[4]]
        return(list(CI[[1]],CI[[2]]))
      }
      else{
        return(list(NA,NA))
      }
    }
  }
  all_threshold_results <- all_threshold_results %>% ungroup()
  CI <- matrix(unlist(apply(all_threshold_results, 1, get_CI)), ncol=2, byrow=T)
  colnames(CI) <- c("low", "high")
  all_threshold_results <- cbind(all_threshold_results, CI)
  all_threshold_results$Carrier = factor(all_threshold_results$Carrier, levels = variant_levels)#c("Non Carrier", "Carrier"))
  all_threshold_results$condition = factor(all_threshold_results$condition, 
                                           levels = c("High LDL (>190 mg/dl)",
                                                      "High LDL (>190 mg/dl) no lipid meds",
                                                      "High LDL ADJ (>190 mg/dl)",
                                                      "High LDL ADJ (>190 mg/dl) no lipid meds",
                                                      "Low LDL (<80 mg/dl)",
                                                      "Low LDL (<80 mg/dl) no lipid meds",
                                                      "Low LDL ADJ (<80 mg/dl)",
                                                      "Low LDL ADJ (<80 mg/dl) no lipid meds",
                                                      "High HDL (>70 mg/dl)",
                                                      "High triglycerides (>200 mg/dl)",
                                                      "High triglycerides (>200 mg/dl) no lipid meds",
                                                      "High triglycerides ADJ (>200 mg/dl)",
                                                      "High triglycerides ADJ (>200 mg/dl) no lipid meds",
                                                      "Monogenic obesity (BMI >30)",
                                                      "MODY (T2D)",
                                                      "MODY (T2D and pre-diabetes)",
                                                      "MODY extended (T2D)",
                                                      "MODY extended (T2D and pre-diabetes)",
                                                      "Lipodystrophy (T2D)",
                                                      "Lipodystrophy (T2D and pre-diabetes)",
                                                      "Neonatal diabetes (T2D)",
                                                      "Neonatal diabetes (T2D and pre-diabetes)"))
  all_threshold_results$condition_simple = factor(all_threshold_results$condition_simple, 
                                                  levels = c("High LDL",
                                                             "Low LDL",
                                                             "High HDL",
                                                             "High TG",
                                                             "Monogenic obesity",
                                                             "MODY",
                                                             "MODY Extended",
                                                             "Lipodystrophy",
                                                             "Neonatal Diabetes"))
  
  
  return(all_threshold_results)
}


penetrance_plot <- function(all_threshold_results, labels = NULL){
  all_threshold_results <- all_threshold_results %>%
    mutate(order = factor(order, rev(levels(order))))
  p <- ggplot(all_threshold_results, aes(order, prop, color=Carrier, fill=Carrier)) +
    geom_point(size=5) + 
    theme_classic() + 
    ylab("Proportion of individuals with the condition") + 
    xlab("") + 
    labs(title="")+
    scale_color_manual(values = colors, limits=c("Non Carrier","Carrier")) +
    scale_fill_manual(values = colors, limits=c("Non Carrier","Carrier")) +
    coord_flip() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=13),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=13),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_blank())+
    ylim(0,1)+
    geom_errorbar(aes(ymin=low, ymax=high,colour=Carrier), width=.7) +
    facet_wrap(data_source~., nrow=1) + 
    theme(panel.spacing = unit(4, "lines"))
  if(!is.null(labels)){
    p <- p + scale_x_discrete(labels=labels)
  }
  return(p)
}


get_pvalues <- function(p_values, trait=NULL, phenotype=NULL, log_axis=NULL, no_lipid_meds=NULL, carrier_type="pathogenic"){
  phenotype_pval <- phenotype
  if(!is.null(trait) & !is.null(phenotype)){
    if(log_axis){phenotype_pval <- paste0(phenotype, "_ln") }
    p_values <- p_values %>%
      filter((type == carrier_type) & (Condition == gsub(" ", '_', trait)) & (Trait == phenotype_pval))
    if(no_lipid_meds){
      p_values <- p_values %>%
        filter(additional_covariates == "covart2d,covarlipidmeds")
    }else{
      p_values <- p_values %>%
        filter(additional_covariates == "covart2d")
    }
  }else{
    p_values <- p_values %>%
      filter((type == carrier_type) &(additional_covariates %in% c(NA,"","covart2d,covarlipidmeds","covart2d")))
  }
  p_values <- p_values %>%
    mutate(lipidmeds = ifelse(additional_covariates %in% c(NA,"","covart2d"),FALSE,TRUE),
           "P value" = formatC(PVALUE, digits = 3),
           "Beta (se)" = paste0(formatC(BETA, digits = 3), " (", formatC(SEBETA, digits = 3), ")"),
           "OR" = paste(formatC(exp(BETA), digits = 3)," (",formatC(exp(BETA-SEBETA), digits = 3),"-",formatC(exp(BETA+SEBETA), digits = 3),")"))
  
  return(p_values)
}

conditions_cont_phenotypes <- rbind.data.frame(
                   list("Obesity", "bmi", "Monogenic obesity (BMI >30)",FALSE,TRUE),
                   list("High HDL", "hdl_mgdl", "High HDL (>70 mg/dl)",FALSE,FALSE),
                   list("High LDL", "ldl_mgdl", "High LDL (>190 mg/dl)",FALSE,FALSE),
                   list("High LDL", "ldl_mgdl", "High LDL (>190 mg/dl) no lipid meds",TRUE,FALSE),
                   list("High LDL", "ldl_mgdl_adj", "High LDL ADJ (>190 mg/dl)",FALSE,FALSE),
                   list("High LDL", "ldl_mgdl_adj", "High LDL ADJ (>190 mg/dl) no lipid meds",TRUE,FALSE),
                   list("Low LDL", "ldl_mgdl", "Low LDL (<80 mg/dl)",FALSE,FALSE),
                   list("Low LDL", "ldl_mgdl", "Low LDL (<80 mg/dl) no lipid meds",TRUE,FALSE),
                   list("Low LDL", "ldl_mgdl_adj", "Low LDL ADJ (<80 mg/dl)",FALSE,FALSE),
                   list("Low LDL", "ldl_mgdl_adj", "Low LDL ADJ (<80 mg/dl) no lipid meds",TRUE,FALSE),
                   list("High triglycerides", "tg_mgdl", "High triglycerides (>200 mg/dl)",FALSE,TRUE),
                   list("High triglycerides", "tg_mgdl", "High triglycerides (>200 mg/dl) no lipid meds",TRUE,TRUE),
                   list("High triglycerides", "tg_mgdl_adj", "High triglycerides ADJ (>200 mg/dl)",FALSE,TRUE),
                   list("High triglycerides", "tg_mgdl_adj", "High triglycerides ADJ (>200 mg/dl) no lipid meds",TRUE,TRUE) ,stringsAsFactors = FALSE,
                   make.row.names=FALSE)
colnames(conditions_cont_phenotypes) <- c("condition", "phenotype", "label", "lipidmeds", "log_value")

carrier_phenotypes_mean_ci_pvalues <- function(data,
                                               AD_gene_list,
                                               conditions_cont_phenotypes,
                                               p_values,
                                               by_gene = FALSE,
                                               carrier_type="pathogenic",
                                               levels_to_use=path_levels,
                                               subsets=FALSE){
  
  p_values <- get_pvalues(p_values, carrier_type=carrier_type)
  if(subsets){
    if(by_gene){
      p_values <- p_values %>%
      select(Condition,Trait,Gene,data_source,lipidmeds,"subset","P value","Beta (se)",OR)
      join_by <- c("condition_simple"="Condition", "phenotype"="Trait", "data_source"="data_source", "lipidmeds"="lipidmeds", "subset"="subset", "Gene"="Gene")
    }else{
      p_values <- p_values %>%
        select(Condition,Trait,data_source,lipidmeds,"subset","P value","Beta (se)",OR)
      join_by <- c("condition_simple"="Condition", "phenotype"="Trait", "data_source"="data_source", "lipidmeds"="lipidmeds", "subset"="subset")
    }
  }else{
    if(by_gene){
      p_values <- p_values %>%
        select(Condition,Trait,Gene,data_source,lipidmeds,"P value","Beta (se)",OR)
      join_by <- c("condition_simple"="Condition", "phenotype"="Trait", "data_source"="data_source", "lipidmeds"="lipidmeds", "Gene"="Gene")
    }else{
      p_values <- p_values %>%
        select(Condition,Trait,data_source,lipidmeds,"P value","Beta (se)",OR)
      join_by <- c("condition_simple"="Condition", "phenotype"="Trait", "data_source"="data_source", "lipidmeds"="lipidmeds")
      
    }
  }
  pheno <- data$pheno
  full_carrier_list <- data$full_carrier_list
  all_means <- NULL
  for(row in seq_len(nrow(conditions_cont_phenotypes))){
    condition <- conditions_cont_phenotypes[row, 1]
    phenotype <- conditions_cont_phenotypes[row, 2]
    label <- conditions_cont_phenotypes[row, 3]
    lipidmeds <- conditions_cont_phenotypes[row, 4]
    if(subsets){
      subset <- conditions_cont_phenotypes[row, 6]
    }
    else{
      subset <- "All"
    }
    if(conditions_cont_phenotypes[row,5]){phenotype_pval <- paste0(phenotype, "_ln") }else{phenotype_pval <- phenotype}
    genes <- subset(AD_gene_list, Trait==condition)$Gene_name
    if(subsets){
      carrier_subset <- annotate_phenotype_with_carriers(pheno, genes, full_carrier_list, condition,
                                                         c("combined_id", phenotype, "age","lipidmeds","subset"),
                                                         levels_to_use,
                                                         collapse_carrier_multiple = TRUE,
                                                         by_gene=by_gene)
    }else{
      carrier_subset <- annotate_phenotype_with_carriers(pheno, genes, full_carrier_list, condition,
                                                         c("combined_id", phenotype, "age","lipidmeds"),
                                                         levels_to_use,
                                                         collapse_carrier_multiple = TRUE,
                                                         by_gene=by_gene)
    }
    if(lipidmeds){
      carrier_subset <- carrier_subset %>% filter(lipidmeds == 0)
    }
    carrier_subset <- carrier_subset %>% filter(Trait == condition & !is.na(Carrier))
    if(subsets){
      if(by_gene){
        df2 <- summarySE(carrier_subset, measurevar=phenotype, groupvars=c("Carrier","data_source","Gene","subset"), na.rm=TRUE)
        to_keep <- c("Carrier", "data_source", "Mean (95% CI)", "Gene", "subset")
      }else{
        df2 <- summarySE(carrier_subset, measurevar=phenotype, groupvars=c("Carrier","data_source","subset"), na.rm=TRUE)
        to_keep <- c("Carrier", "data_source", "Mean (95% CI)", "subset")
      }
    }else{
      if(by_gene){
        df2 <- summarySE(carrier_subset, measurevar=phenotype, groupvars=c("Carrier","data_source","Gene"), na.rm=TRUE)
        to_keep <- c("Carrier", "data_source", "Mean (95% CI)", "Gene")
      }else{
        df2 <- summarySE(carrier_subset, measurevar=phenotype, groupvars=c("Carrier","data_source"), na.rm=TRUE)
        to_keep <- c("Carrier", "data_source", "Mean (95% CI)")
      }
    }
    df2$phenotype_value = df2[,phenotype]
    means <- df2 %>%
            mutate("Mean (95% CI)"= paste0(formatC(phenotype_value, digits = 2, width = 0, format = "f"), "   (", formatC(phenotype_value - ci, digits = 2, width = 0, format = "f"), "-", formatC(phenotype_value + ci, digits = 2, width = 0, format = "f"), ")")) %>%
            select(to_keep) %>% 
      mutate(Condition = label, phenotype=phenotype_pval, condition_simple=gsub(" ", '_', condition), lipidmeds=lipidmeds)

    all_means <- rbind(all_means, means)
  }
  if(!by_gene){
    all_means <- all_means %>%
      pivot_wider(names_from = c("Carrier"),values_from=c("Mean (95% CI)"),names_prefix="Mean (95% CI) ") %>%  #,values_from=c("Mean (95% CI)","P value","Beta (se)")) %>%
      left_join(p_values,by=join_by) %>%
      select(-OR) %>%
      pivot_wider(names_from = c("data_source"),values_from=c("Mean (95% CI) Carrier","Mean (95% CI) Non Carrier",
                                                              "P value",
                                                              "Beta (se)"))#,names_prefix="Mean (95% CI)")
  }
  else{
    all_means <- all_means %>%
      left_join(p_values,by=join_by) %>% 
      select(-OR) %>%
      pivot_wider(names_from = c("data_source"),values_from=c("Mean (95% CI)",
                                                              "P value",
                                                              "Beta (se)"))
  }
  return(all_means)
}

plot_carrier_phenotypes_one_condition <- function(data,
                                                  AD_gene_list,
                                                  trait,
                                                  phenotype,
                                                  log_axis,
                                                  ylab,
                                                  title,
                                                  threshold=NULL,
                                                  color_col = "ClinVarOrLOF",
                                                  p_values = NULL,
                                                  above_or_below = "above",
                                                  no_lipid_meds = FALSE){
  if(!is.null(p_values)){
    p_values <- get_pvalues(p_values, trait, phenotype, log_axis, no_lipid_meds) %>% select(PVALUE, data_source)
  }
  pheno <- data$pheno
  full_carrier_list <- data$full_carrier_list
  genes <- subset(AD_gene_list, Trait==trait)$Gene_name
  carrier_subset <- annotate_phenotype_with_carriers(pheno, genes, full_carrier_list, trait,
                                                     c("combined_id", phenotype, "age","lipidmeds"),
                                                     path_levels,
                                                     collapse_carrier_multiple = TRUE)
  if(no_lipid_meds){
    carrier_subset <- carrier_subset %>% filter(lipidmeds == 0)
  }
  carrier_subset <- carrier_subset %>%
    mutate(ClinVarOrLOF = replace(ClinVarOrLOF, Carrier=="Non Carrier","Non Carrier")) %>%
    mutate(Carrier_abr = Carrier) %>%
    mutate(Carrier_abr = replace(Carrier_abr, Carrier_abr=="Non Carrier","NC")) %>%
    mutate(Carrier_abr = replace(Carrier_abr, Carrier_abr=="Carrier","C")) %>%
    mutate(Carrier_abr = factor(Carrier_abr, levels = c("NC","C"))) %>%
    filter(Trait == trait & !is.na(Carrier))
  carrier_subset <- carrier_subset %>% left_join(p_values, by=c("data_source"="data_source")) %>%
    mutate(data_source = paste(data_source, paste0("P=", formatC(PVALUE, digits = 3)), sep = "\n"))
    
  data_sources <- sort(unique(carrier_subset$data_source))
  carrier_subset$data_source = factor(carrier_subset$data_source, 
                                      levels = data_sources)
  p <- ggplot(carrier_subset, aes_string("Carrier_abr", phenotype)) +
    geom_quasirandom(aes_string(colour=color_col), size=1, groupOnX=TRUE) + 
    stat_summary(fun.data=mean_cl_normal, size=0.6,geom="errorbar")  +
    stat_summary(fun.data=mean_cl_normal, size=3, geom="point")  +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_blank()) + 
    ylab(ylab) + xlab("") + labs(title=title)+
    facet_wrap(~data_source, nrow=1) + 
    scale_colour_manual(values=colors,
                        limits=subset(variant_levels, variant_levels %in% unique(as.character(carrier_subset %>% pull(color_col))))) +
    guides(colour = guide_legend(override.aes = list(size=3)))

  if(log_axis & phenotype == "bmi"){
    p <- p + scale_y_log10(breaks = c(10, 20, 30, 40, 50, 60, 70, 80)) #+ a
  }else if(log_axis){
    p <- p + scale_y_log10(breaks = c(10, 100, 1000, 5000)) #+ a
  }
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  legend <- get_legend(p)
  p <- p + theme(legend.position="none")
  gb <- ggplot_build(p)
  g <- ggplot_gtable(gb)
  ranges <- gb$layout$panel_params
  if(!is.null(threshold)){
    data2npc <- function(x, range) scales::rescale(c(range, x), c(0,1))[-c(1,2)]
    if(log_axis){
      start <- c(data2npc(0, ranges[[1]][["x.range"]]), 
                 data2npc(log10(threshold), ranges[[1]][["y.range"]]))
      
      end <- c(data2npc(6, ranges[[2]][["x.range"]]), 
               data2npc((ranges[[2]][["y.range"]][[2]] - log10(threshold))+ranges[[2]][["y.range"]][[1]], ranges[[2]][["y.range"]])+0.05)
    }else{
      start <- c(data2npc(0, ranges[[1]][["x.range"]]), 
                 data2npc(threshold, ranges[[1]][["y.range"]]))
      
      end <- c(data2npc(6, ranges[[2]][["x.range"]]), 
               data2npc(ranges[[2]][["y.range"]][[2]] - threshold, ranges[[2]][["y.range"]]))
    }
  
    if(above_or_below == "above"){
      g <- gtable_add_grob(g, rectGrob(x=0, y=start[2], width = 2.1, height=end[2]-0.05,
                                     just=c("left","bottom"),gp = gpar( fill="blue", col="blue", alpha=0.12,lwd=0)), t=7, b=9, l=5)
    }else{
      g <- gtable_add_grob(g, rectGrob(x=0, y=0, width = 2.1, height=start[2],
                                       just=c("left","bottom"),gp = gpar( fill="blue", col="blue", alpha=0.12,lwd=0)), t=7, b=9, l=5)
      
    }
    g$layout$clip <- "off"
  }
  return(list(p=g,legend=legend))
}


plot_carrier_phenotypes_by_condition <- function(data, AD_gene_list, p_values=NULL, no_lipid_meds=FALSE){
  p_high_ldl <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                      "High LDL", "ldl_mgdl", FALSE, "LDL (mg/dl)", "High LDL",
                                                      threshold_for_plot$threshold[which(threshold_for_plot$Trait=="High LDL")],
                                                      color_col = "Carrier", p_values=p_values, no_lipid_meds=no_lipid_meds)
  p_low_ldl <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                     "Low LDL", "ldl_mgdl", FALSE, "LDL (mg/dl)", "Low LDL",
                                                     threshold_for_plot$threshold[which(threshold_for_plot$Trait=="Low LDL")],
                                                     color_col = "Carrier", p_values=p_values, above_or_below = "below"
                                                    , no_lipid_meds=no_lipid_meds)
  p_high_hdl <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                      "High HDL", "hdl_mgdl", FALSE, "HDL (mg/dl)", "High HDL",
                                                      threshold_for_plot$threshold[which(threshold_for_plot$Trait=="High HDL")],
                                                      color_col = "Carrier", p_values=p_values, no_lipid_meds=no_lipid_meds)
  p_high_tg <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                     "High triglycerides", "tg_mgdl", TRUE, "Triglycerides (mg/dl)", "High TG",
                                                     threshold_for_plot$threshold[which(threshold_for_plot$Trait=="High triglycerides")],
                                                     color_col = "Carrier", p_values=p_values, no_lipid_meds=no_lipid_meds)
  p_obesity <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                     "Obesity", "bmi", TRUE, expression("Body mass index (kg/m"^"2"*")"), "Obesity",
                                                     threshold_for_plot$threshold[which(threshold_for_plot$Trait=="Obesity")],
                                                     color_col = "Carrier", p_values=p_values, no_lipid_meds=no_lipid_meds)

  if(no_lipid_meds){
    return(grid.arrange(p_high_ldl$p,
                        p_low_ldl$p,
                        p_high_hdl$p,
                        p_high_tg$p, nrow=1))
  }else{
    return(grid.arrange(p_high_ldl$p,
                     p_low_ldl$p,
                     p_high_hdl$p,
                     p_high_tg$p,
                     p_obesity$p, nrow=1)) #,ncol=5,common.legend = TRUE,legend="right",widths = c(1,1,1,1,1)))

  }
}


plot_carrier_phenotypes_by_condition_supp <- function(data, AD_gene_list, p_values=p_values){
  p_high_ldl <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                      "High LDL", "ldl_mgdl", FALSE, "LDL (mg/dl)", "High LDL", p_values=p_values)
  p_low_ldl <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                     "Low LDL", "ldl_mgdl", FALSE, "LDL (mg/dl)", "Low LDL", p_values=p_values, above_or_below = "below")
  p_high_hdl <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                      "High HDL", "hdl_mgdl", FALSE, "HDL (mg/dl)", "High HDL", p_values=p_values)
  p_high_tg <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                     "High triglycerides", "tg_mgdl", TRUE, "Triglycerides (mg/dl)", "High TG", p_values=p_values)
  p_obesity <- plot_carrier_phenotypes_one_condition(data, AD_gene_list,
                                                     "Obesity", "bmi", TRUE, expression("Body mass index (kg/m"^"2"*")"), "Obesity", p_values=p_values)
  

  leg <- p_high_ldl$legend
  return(grid.arrange(p_high_ldl$p,
                   p_low_ldl$p,
                   p_high_hdl$p,
                   p_high_tg$p,
                   p_obesity$p,
                   leg,
                   ncol=2,nrow=3))
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE, binomial=FALSE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=TRUE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  length1 <- function (x, value) {
    length(which(x==value))
  }
  measurevar_param <- enquo(measurevar)
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  if(binomial){
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                         .fun = function(xx, col) {
                            c(N    = length2(xx[[col]], na.rm=na.rm),
                              X    = length1(xx[[col]],1),
                              mean = mean   (xx[[col]], na.rm=na.rm),
                              sd   = sd     (xx[[col]], na.rm=na.rm)
                            )
                          },
                         measurevar
    )
  }else{
    datac <- plyr::ddply(data, groupvars, .drop=.drop,
                         .fun = function(xx, col) {
                            c(N    = length2(xx[[col]], na.rm=na.rm),
                              mean = mean   (xx[[col]], na.rm=na.rm),
                              sd   = sd     (xx[[col]], na.rm=na.rm)
                            )
                          },
                         measurevar
    )
  }
  
  
  # Rename the "mean" column    

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  if(binomial){
    get_CI <- function(datac){
      if(!is.na(datac[[3]]) & as.numeric(datac[[3]])>0 & !is.na(datac[[2]]) & as.numeric(datac[[2]])>0){
        CI <- binom.test(as.numeric(datac[[3]]), as.numeric(datac[[2]]), conf.level = 0.95)[[4]]
        return(list(CI[[1]],CI[[2]]))
      }
      else{
        return(list(NA,NA))
      }
    }
    CI <- matrix(unlist(apply(datac, 1, get_CI)), ncol=2, byrow=T)
    colnames(CI) <- c("low", "high")
    datac <- cbind(datac, CI)
  }else{
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    datac$low = datac$mean - datac$ci
    datac$high = datac$mean + datac$ci
  }
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  
  return(datac)
}


get_or <- function(data, group){
  if(group == 1){
    return(tibble(group=1,OR=1,low=1,high=1))
  }
  data$prs = data$prs
  data <- data %>%
    mutate(prs_under_equal = case_when(percentile==1 ~ 0,
                                       percentile==group ~ 1))
  logit.full <- glm(phenotype_by_threshold~prs_under_equal+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=data, family=binomial)

  # Calculate OR for specific increment step of continuous variable
  return(tibble(group = group,
              OR=exp(coef(logit.full)['prs_under_equal']), 
              low = exp(confint.default(logit.full, level = 0.95)['prs_under_equal',1]),
              high = exp(confint.default(logit.full, level = 0.95)['prs_under_equal',2])))
}

PRS_plots_one_condition <- function(my_data,
                                    PRS_fp,
                                    phenotype,
                                    trait,
                                    threshold,
                                    prs_label,
                                    trait_label,
                                    trait_label2,
                                    p_values=NULL,
                                    restrict_lipid_meds = FALSE,
                                    threshold_under = FALSE,
                                    log_breaks = NULL,
                                    age_cutoff = NULL,
                                    by_gene = FALSE,
                                    genes_to_keep = NULL)
{

  pheno <- my_data$pheno
  full_carrier_list <- my_data$full_carrier_list

  pheno["phenotype_value"] <- pheno[phenotype]
  PRS <- read.table(PRS_fp, header=T)
  PRS$IID = as.character(PRS$IID)
  PRS <- pheno %>% left_join(PRS, by = c("combined_id"="IID"))
  PRS$prs = scale(PRS$SCORE)

  trait_title <- trait
  if(restrict_lipid_meds){
    PRS <- PRS %>% filter(lipidmeds == 0)
    trait_title <- paste0(trait_title, " restricted to\nno lipid medication")
  }
  if(threshold_under){
    PRS$prs = -1*PRS$prs
  }
  PRS %>% filter(!is.na(phenotype_value) & !is.na(prs))-> PRS
  if(is.null(threshold) & (trait == "MODY")){
    PRS %>% mutate(phenotype_by_threshold = phenotype_value==1) -> PRS
  }else{
    if(threshold_under){
      PRS %>% mutate(phenotype_by_threshold = phenotype_value<threshold) -> PRS
    }else{
      PRS %>% mutate(phenotype_by_threshold = phenotype_value>=threshold) -> PRS
    }
  }
  carrier_subset_phenotype_value <- annotate_phenotype_with_carriers(PRS,
                                                                     subset(AD_gene_list, Trait==trait)$Gene_name,
                                                                     full_carrier_list,
                                                                     trait,
                                                                     c("combined_id","phenotype_value","age","prs","phenotype_by_threshold","SEX","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"),
                                                                     path_levels,
                                                                     collapse_carrier_multiple = TRUE,
                                                                     by_gene = by_gene )
  if(!is.null(age_cutoff)){
    carrier_subset_phenotype_value <- carrier_subset_phenotype_value %>% filter((age >= age_cutoff) | (Carrier=="Carrier"))
  }

  carrier_subset_phenotype_value %>%
    mutate(quantile = ntile(prs, 4) ) %>% 
    mutate(decile = ntile(prs, 10) ) %>% 
    mutate(percentile = ntile(prs, 100) ) %>% 
    mutate(forty = ntile(prs, 40) ) %>% 
    filter(!is.na(phenotype_value) & !is.na(prs))-> carrier_subset_phenotype_value
  
  if(by_gene){
    carrier_subset_phenotype_value <- carrier_subset_phenotype_value %>% select("combined_id", "quantile", "decile", "percentile", "forty", "phenotype_value", "age", "prs", "phenotype_by_threshold", "SEX", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Carrier", "Gene")
  }else{
    carrier_subset_phenotype_value <- carrier_subset_phenotype_value %>% select("combined_id", "quantile", "decile", "percentile", "forty", "phenotype_value", "age", "prs", "phenotype_by_threshold", "SEX", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Carrier")
  }
  carrier_subset_phenotype_value <- carrier_subset_phenotype_value %>% mutate(Carrier_and_prs_quantile = paste(Carrier, quantile))

  top_1_per <- carrier_subset_phenotype_value %>% filter(percentile == 100 | Carrier == "Carrier")
  non_top_1_per <- carrier_subset_phenotype_value %>% filter(percentile < 100 | Carrier == "Carrier")
  interquartile <- carrier_subset_phenotype_value %>% mutate(grouping=ifelse((quantile == 2) | (quantile == 3), 'interquartile', 'outer')) %>% mutate(grouping=ifelse(percentile == 100, 'top_1', grouping)) %>% filter((quantile == 2) | (quantile == 3) | (percentile == 100)) %>% filter(Carrier != "Carrier")
  top_1_per_vs_99 <- carrier_subset_phenotype_value %>% filter(!is.na(percentile)) %>% mutate(grouping=ifelse(percentile < 100, 'top_1_per', 'non_top_1')) %>% filter(Carrier != "Carrier")
  glm_interquartile_1_per_pvalues <- NULL #data.frame(Condition=character(),
                                          #     Gene=character(), 
                                           #    `Top 1% of gePS vs Interquartile range (25-75%) estimate`=double(), 
                                            #   `Top 1% of gePS vs Interquartile range (25-75%) pval`=double(), 
                                             #  `Monogenic carriers vs Top 1% gePS estimate`=double(), 
                                              # `Monogenic carriers vs Top 1% gePS	pval`=double(), 
                                               #stringsAsFactors=FALSE) 
  if(trait == "MODY"){
    fit <- glm(phenotype_value~Carrier+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=top_1_per, family=binomial)
    top_1_per_p <- coef(summary(fit))[2, c(1, 4)]
    fit <- glm(phenotype_value~grouping+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=interquartile, family=binomial)
    interquartile_p <- coef(summary(fit))[2, c(1, 4)]
    glm_interquartile_1_per_pvalues <- glm_interquartile_1_per_pvalues %>% rbind(c("MODY", NA, interquartile_p[1], top_1_per_p[1], interquartile_p[2], top_1_per_p[2]))
    
    if(by_gene){
      fit <- glm(phenotype_value~Carrier+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=top_1_per%>% filter((Gene == "GCK") | (Gene == "Non Carrier")), family=binomial)
      top_1_per_p <- coef(summary(fit))[2, c(1, 4)]
      fit <- glm(phenotype_value~grouping+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=interquartile%>% filter((Gene == "GCK") | (Gene == "Non Carrier")), family=binomial)
      interquartile_p <- coef(summary(fit))[2, c(1, 4)]
      glm_interquartile_1_per_pvalues <- glm_interquartile_1_per_pvalues %>% rbind(c("MODY", "GCK", interquartile_p[1], top_1_per_p[1], interquartile_p[2], top_1_per_p[2]))
      
      fit <- glm(phenotype_value~Carrier+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=top_1_per%>% filter((Gene == "HNF1A") | (Gene == "Non Carrier")), family=binomial)
      top_1_per_p <- coef(summary(fit))[2, c(1, 4)]
      fit <- glm(phenotype_value~grouping+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=interquartile%>% filter((Gene == "HNF1A") | (Gene == "Non Carrier")), family=binomial)
      interquartile_p <- coef(summary(fit))[2, c(1, 4)]
      glm_interquartile_1_per_pvalues <- glm_interquartile_1_per_pvalues %>% rbind(c("MODY", "HNF1A", interquartile_p[1], top_1_per_p[1], interquartile_p[2], top_1_per_p[2]))
    }
  }else{
    fit <- glm(phenotype_value~Carrier+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=top_1_per)
    top_1_per_p <- coef(summary(fit))[2, c(1, 4)]
    fit <- glm(phenotype_value~grouping+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=interquartile)
    interquartile_p <- coef(summary(fit))[2, c(1, 4)]
    glm_interquartile_1_per_pvalues <- glm_interquartile_1_per_pvalues %>% rbind(c(trait, NA, interquartile_p[1], top_1_per_p[1], interquartile_p[2], top_1_per_p[2]))
  }
  if(trait=="MODY"){
    df2 <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                     groupvars=c("percentile"),
                      binomial=TRUE)
    df4 <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                     groupvars=c("forty"),
                     binomial=TRUE)
  }else{
    df2 <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                     groupvars=c("percentile"))
    df4 <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                     groupvars=c("forty"))
  }
  
  if(by_gene){
    carrier_subset_phenotype_value <- carrier_subset_phenotype_value %>%
      filter(Gene != "Non Carrier")
    if(!is.null(genes_to_keep)){
      carrier_subset_phenotype_value <- carrier_subset_phenotype_value %>%
        filter(Gene %in% genes_to_keep)
    }
    if(trait=="MODY"){
      one_mean <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                          groupvars=c("Gene"),
                          binomial=TRUE)
    }else{
      one_mean <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                            groupvars=c("Gene"))
    }

    one_mean <- one_mean %>%
      mutate(percentile = Gene, Carrier="Carrier") %>%
      select(-Gene)
  
    
  }else{
    if(trait=="MODY"){
      one_mean <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                            groupvars=c("Carrier"),
                            binomial=TRUE)
      head(one_mean)
    }else{
      one_mean <- summarySE(carrier_subset_phenotype_value, measurevar="phenotype_value", 
                       groupvars=c("Carrier"))
    }
    one_mean <- one_mean %>%
      mutate(percentile = "C")  %>% 
      filter(Carrier == "Carrier")
  }
  one_mean$color = "Carrier"
  df2$color = "gePS percentile"
  df2 <- rbind(df2, one_mean %>% select(-Carrier))
  df4$color = "gePS percentile"
  p_percentile <- ggplot(na.omit(df2 %>% filter(color != "Carrier")), aes(x=factor(percentile, levels = c(seq(1, 100))), y=phenotype_value, color=color)) +
    geom_point(size=2, position=position_dodge(.4))+
    xlab(paste0(prs_label, " percentile\n")) +
    ylab(trait_label) + 
    theme_classic() +
        theme(axis.text.y=element_text(size=12),
          axis.title=element_text(size=12),
          legend.text=element_text(size=12),
          axis.text.x =element_blank(),
          axis.ticks.x = element_blank(),
          legend.title=element_blank(),
          plot.margin=margin(t=25)) +
    geom_errorbar(aes(ymin=low, ymax=high), width=.1,
                  position=position_dodge(.9)) +
    ylim(min(min(df2$low,na.rm=TRUE),min(one_mean$low,na.rm=TRUE)),
         max(max(df2$high,na.rm=TRUE),max(one_mean$high,na.rm=TRUE))) +
  scale_colour_manual(values=c("gePS percentile"="black","Carrier"="#d6604d")) +
    theme(legend.position="none")
  
  p_percentile_carrier <- ggplot(na.omit(df2 %>% filter(color == "Carrier")), aes(x=percentile, y=phenotype_value, color=color)) +
    geom_point(size=2.2, position=position_dodge(.4))+
    xlab(paste0("Clinically significant\nvariant carrier")) +
    ylab("") + 
    theme_classic() +
    theme(axis.text.y=element_text(size=12),
          axis.title=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_blank(),
          plot.margin=margin(t=20)) +
    geom_errorbar(aes(ymin=low, ymax=high), width=.2,size=1.2,
                  position=position_dodge(.9)) +
    ylim(min(min(df2$low,na.rm=TRUE),min(one_mean$low,na.rm=TRUE)),
         max(max(df2$high,na.rm=TRUE),max(one_mean$high,na.rm=TRUE))) +
    scale_colour_manual(values=c("gePS percentile"="black","Carrier"="#d6604d")) +
    theme(legend.position="none")
  if(!by_gene){
    p_percentile_carrier <- p_percentile_carrier +
      theme(axis.text.x =element_blank(),
          axis.ticks.x = element_blank())
  }
  p_percentile_non_carrier <- ggplot(na.omit(df4 %>% filter(color != "Carrier")), aes(x=factor(forty, levels = c(seq(1, 40))), y=phenotype_value, color=color)) +
    geom_point(size=3, position=position_dodge(.4))+
    xlab(paste0(prs_label, " gePS grouping\n")) +
    ylab(trait_label) + 
    theme_classic() +
    theme(axis.text.y=element_text(size=12),
          axis.title=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_blank()) +
    geom_errorbar(aes(ymin=low, ymax=high), width=.1,
                  position=position_dodge(.9)) +
    scale_colour_manual(values=c("gePS percentile"="black","Carrier"="#d6604d")) +
    theme(legend.position="none")
  
  p_percentile_non_carrier_no_bars <- ggplot(na.omit(df4 %>% filter(color != "Carrier")), aes(x=factor(forty, levels = c(seq(1, 40))), y=phenotype_value, color=color)) +
    geom_point(size=3, position=position_dodge(.4))+
    xlab(paste0(prs_label, " gePS grouping\n")) +
    ylab(trait_label) + 
    theme_classic() +
    theme(axis.text.y=element_text(size=12),
          axis.title=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_blank()) +
    scale_colour_manual(values=c("gePS percentile"="black","Carrier"="#d6604d")) +
    theme(legend.position="none")
  
  
  p_percentile <- ggarrange(p_percentile, p_percentile_carrier, widths = c(1, 0.2), nrow=1, align="hv", legend="none")

  p_ecdf <- ggplot(na.omit(carrier_subset_phenotype_value %>% filter(quantile %in% c(1, 4))),
                   aes(phenotype_value,
                      group=as.factor(Carrier_and_prs_quantile),
                      color=as.factor(Carrier_and_prs_quantile))) + 
    stat_ecdf(size = 1.2) + ylab("Empirical CDF") + 
    xlab(trait_label2) +     
    scale_color_discrete(name = "", labels = c("Carrier gePS bottom quartile", 
                                               "Carrier gePS top quartile",
                                               "Non carrier gePS bottom quartile", 
                                               "Non carrier gePS top quartile"))
  if(!is.null(log_breaks)){
    p_ecdf <- p_ecdf + scale_x_log10(breaks = log_breaks) +
      annotation_logticks(sides="b")
  }
  p_ecdf <- p_ecdf + theme_classic() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_blank(),
          plot.title = element_text(size=13,hjust=0.5))+
    ggtitle(trait_title)


  df2$Condition = trait
  carrier_subset_phenotype_value_carrier_only <- carrier_subset_phenotype_value %>% filter(Carrier == "Carrier")

  glm_results <- summary(glm(phenotype_value~prs+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data= carrier_subset_phenotype_value_carrier_only))
  glm_results <- as.data.frame(glm_results$coefficients) %>% rownames_to_column('Variable')
  glm_results$Condition = trait
  glm_results$restrict_lipid_meds = restrict_lipid_meds

  if(trait == "MODY"){
    glm_results_all <- NULL
  }else{
    glm_results_all <- summary(glm(phenotype_value~factor(Carrier, levels=c("Non Carrier", "Carrier"))*prs+age+SEX+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data=carrier_subset_phenotype_value))
    glm_results_all <- as.data.frame(glm_results_all$coefficients) %>% rownames_to_column('Variable')
    glm_results_all$Condition = trait
    glm_results_all$restrict_lipid_meds = restrict_lipid_meds
  }
  
  return(list("p_percentile" = p_percentile,
              "p_percentile_non_carrier" = p_percentile_non_carrier,
              "p_percentile_non_carrier_no_bars" = p_percentile_non_carrier_no_bars,
              "p_ecdf" = p_ecdf, 
              "glm_results" = glm_results,
              "glm_results_all" = glm_results_all,
              "percentile_table" = df2,
              "glm_interquartile_1_per_pvalues" = glm_interquartile_1_per_pvalues))
}


PRS_plots <- function(data, age_cutoff=NULL){
  low_ldl_p <- PRS_plots_one_condition(data,
                                       get_PRS_path("ldl_mgdl"),
                                       "ldl_mgdl", "Low LDL",
                                       threshold_for_plot$threshold[which(threshold_for_plot$Trait=="Low LDL")],
                                       "UKBB LDL gePS", "Mean LDL (mg/dL)", "LDL (mg/dL)", threshold_under = TRUE,
                                       age_cutoff = age_cutoff)
  high_ldl_p <- PRS_plots_one_condition(data,
                                        get_PRS_path("ldl_mgdl"),
                                        "ldl_mgdl", "High LDL",
                                        threshold_for_plot$threshold[which(threshold_for_plot$Trait=="High LDL")],
                                        "UKBB LDL gePS", "Mean LDL (mg/dL)", "LDL (mg/dL)",
                                        age_cutoff = age_cutoff)
  low_ldl_lipid_meds_p <- PRS_plots_one_condition(data,
                                                  get_PRS_path("ldl_mgdl"),
                                                  "ldl_mgdl", "Low LDL",
                                                  threshold_for_plot$threshold[which(threshold_for_plot$Trait=="Low LDL")],
                                                  "UKBB LDL gePS", "Mean LDL (mg/dL)", "LDL (mg/dL)",
                                                  threshold_under = TRUE,
                                                  restrict_lipid_meds = TRUE,
                                                  age_cutoff = age_cutoff)
  high_ldl_lipid_meds_p <- PRS_plots_one_condition(data,
                                                   get_PRS_path("ldl_mgdl"),
                                                   "ldl_mgdl", "High LDL",
                                                   threshold_for_plot$threshold[which(threshold_for_plot$Trait=="High LDL")],
                                                   "UKBB LDL gePS", "Mean LDL (mg/dL)", "LDL (mg/dL)",
                                                   restrict_lipid_meds = TRUE,
                                                   age_cutoff = age_cutoff)
  high_hdl_p <- PRS_plots_one_condition(data,
                                        get_PRS_path("hdl_mgdl"),
                                        "hdl_mgdl", "High HDL",
                                        threshold_for_plot$threshold[which(threshold_for_plot$Trait=="High HDL")],
                                        "UKBB HDL gePS", "Mean HDL (mg/dL)", "HDL (mg/dL)",
                                        age_cutoff = age_cutoff)
  high_tg_p <- PRS_plots_one_condition(data,
                                       get_PRS_path("tg_mgdl"),
                                       "tg_mgdl", "High triglycerides",
                                       threshold_for_plot$threshold[which(threshold_for_plot$Trait=="High triglycerides")],
                                       "UKBB TG gePS", "Mean TG (mg/dL)", "TG (mg/dL)",
                                       age_cutoff = age_cutoff)
  obesity_p <- PRS_plots_one_condition(data,
                                       get_PRS_path("bmi"),
                                       "bmi", "Obesity",
                                       threshold_for_plot$threshold[which(threshold_for_plot$Trait=="Obesity")],
                                       "UKBB BMI gePS", expression("Mean BMI (kg/m"^"2"*")"),
                                       expression("BMI (kg/m"^"2"*")"),
                                       log_breaks = c(10,20,30,40,50,60,70),
                                       age_cutoff = age_cutoff)
  
  mody_p <- PRS_plots_one_condition(data,
                                    get_PRS_path("t2d"),
                                    "t2d", "MODY",
                                    threshold=NULL,
                                    "UKBB T2D gePS", "Proportion diabetes",
                                    "Proportion diabetes",
                                    age_cutoff = age_cutoff)
  
  legend <- get_legend(
    # create some space to the left of the legend
    obesity_p$p_ecdf + 
      guides(color = guide_legend(nrow = 4)) +
      theme(legend.position = "top")
  )
  prow1 <- plot_grid(low_ldl_p$p_percentile + theme(legend.position="none"),
                     #low_ldl_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,scale=0.93)
  prow2 <- plot_grid(high_ldl_p$p_percentile + theme(legend.position="none"),
                     #high_ldl_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,scale=0.93)
  prow3 <- plot_grid(high_hdl_p$p_percentile + theme(legend.position="none"),
                     #high_hdl_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,scale=0.93)
  prow4 <- plot_grid(high_tg_p$p_percentile + theme(legend.position="none"),
                     #high_tg_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,scale=0.93)
  prow5 <- plot_grid(obesity_p$p_percentile + theme(legend.position="none"),
                     #obesity_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,scale=0.93)
  prow6 <- plot_grid(mody_p$p_percentile + theme(legend.position="none"),
                     #mody_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,scale=0.93)
  
  p1 <- plot_grid(prow1, prow2, prow3, prow4, prow5, prow6, nrow = 6,
                  rel_heights = c(1,1,1,1,1,2.2),
                  labels = c("A - Low LDL", "B - High LDL", "C - High HDL", "D - High TG", "E - Obesity","F - Diabetes"),
                  label_x = 0.01, label_y = 0.9, hjust = -0.1, vjust = -0.1)
  
  p1_supp <- prow6

  
  prow1 <- plot_grid(low_ldl_p$p_ecdf + theme(legend.position="none"),
                     high_ldl_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,
                     labels = c("A", "B"),scale=0.93)
  prow2 <- plot_grid(#high_hdl_p$p_percentile + theme(legend.position="none"),
                     high_hdl_p$p_ecdf + theme(legend.position="none"),
                     high_tg_p$p_ecdf + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,
                     labels = c("C", "D"),scale=0.93)
  prow3 <- plot_grid(obesity_p$p_ecdf + theme(legend.position="none"),
                     legend,
                     align = 'vh',
                     hjust = -1,
                     nrow = 1,
                     labels = c("E",""),scale=0.93)
  
  p2 <- plot_grid(prow1, prow2, prow3, nrow = 3,
                  rel_heights = c(1,1,1),
                  label_x = 0, label_y = 1, hjust = -0.1, vjust = -0.1)
  
  
  
  percentile_table <- rbind(low_ldl_p$percentile_table,
                            high_ldl_p$percentile_table,
                            high_hdl_p$percentile_table,
                            high_tg_p$percentile_table,
                            obesity_p$percentile_table)
  
  glm_results <- rbind(low_ldl_p$glm_results,
                       low_ldl_lipid_meds_p$glm_results,
                       high_ldl_p$glm_results,
                       high_ldl_lipid_meds_p$glm_results,
                       high_hdl_p$glm_results,
                       high_tg_p$glm_results,
                       obesity_p$glm_results,
                       mody_p$glm_results)
  
  glm_results_all <- rbind(low_ldl_p$glm_results_all,
                           low_ldl_lipid_meds_p$glm_results_all,
                           high_ldl_p$glm_results_all,
                           high_ldl_lipid_meds_p$glm_results_all,
                           high_hdl_p$glm_results_all,
                           high_tg_p$glm_results_all,
                           obesity_p$glm_results_all)
  
  glm_interquartile_1_per_pvalues_all <- rbind(
    low_ldl_p$glm_interquartile_1_per_pvalues,
    low_ldl_lipid_meds_p$glm_interquartile_1_per_pvalues,
    high_ldl_p$glm_interquartile_1_per_pvalues,
    high_ldl_lipid_meds_p$glm_interquartile_1_per_pvalues,
    high_hdl_p$glm_interquartile_1_per_pvalues,
    high_tg_p$glm_interquartile_1_per_pvalues,
    obesity_p$glm_interquartile_1_per_pvalues,
    mody_p$glm_interquartile_1_per_pvalues
    )
  colnames(glm_interquartile_1_per_pvalues_all) <- c(
    "Condition", 
    "Gene",
    "Top 1% of gePS vs Interquartile range (25-75%) estimate", 
    "Top 1% of gePS vs Interquartile range (25-75%) pval",
    "Monogenic carriers vs Top 1% gePS estimate",
    "Monogenic carriers vs Top 1% gePS	pval")
  
  return(list(p1 = p1, p2 = p2,percentile_table=percentile_table,glm_results=glm_results,glm_results_all=glm_results_all, p_percentile_non_carrier_mody = mody_p$p_percentile_non_carrier,p_percentile_non_carrier_mody_no_bars = mody_p$p_percentile_non_carrier_no_bars,
              p_percentile_non_carrier_obesity = obesity_p$p_percentile_non_carrier,p_percentile_non_carrier_obesity_no_bars = obesity_p$p_percentile_non_carrier_no_bars, glm_interquartile_1_per_pvalues_all=glm_interquartile_1_per_pvalues_all))
}




plot_carrier_ancestry_by_condition_proportion <- function(data){
  pheno <- data$pheno
  pheno_annotate_carrier_more_severe <- data$pheno_annotate_carrier_more_severe
  
  pheno_annotate_carrier_more_severe$Trait[which(pheno_annotate_carrier_more_severe$Trait == "Obesity")] = "Monogenic obesity"
  pheno_annotate_carrier_more_severe$Ancestry = factor(pheno_annotate_carrier_more_severe$Ancestry, levels=ancestry_order)
  pheno_annotate_carrier_more_severe$Ancestry = factor(pheno_annotate_carrier_more_severe$Ancestry, levels=ancestry_order)
  pheno_annotate_carrier_more_severe <- pheno_annotate_carrier_more_severe %>%
    mutate(Trait = replace(Trait, Trait == "MODY Extended", "MODY extended")) %>%
    mutate(Trait = replace(Trait, Trait == "Neonatal Diabetes", "Neonatal diabetes"))
  pheno_annotate_carrier_more_severe$Trait = factor(pheno_annotate_carrier_more_severe$Trait, levels=disease_order)
  df_t <- pheno_annotate_carrier_more_severe %>%
    mutate(ClinVarOrLOF = ifelse(ClinVarOrLOF %in% c("Review","ClinVar"), "ClinVar and Reviews", ClinVarOrLOF)) %>%
    filter(!is.na(Trait)) %>%
    filter(!((Trait == "MODY Extended")&(Gene %in% c("GCK","HNF1A","HNF1B","PDX1","HNF4A")))) %>%
    select("combined_id","Variant","Trait","Ancestry","ClinVarOrLOF") %>%
    unique() %>%
    group_by(Trait, Ancestry, ClinVarOrLOF) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    group_by(Trait, ClinVarOrLOF) %>%
    mutate((Frequency <- n/sum(n))*100) %>%
    arrange(ClinVarOrLOF, Trait, Ancestry)

  colnames(df_t) <- c("Condition", "Ancestry", "ClinVarOrLOF", "n", "Frequency")

  p <- ggplot(data=df_t, aes(x=Condition, y=Frequency, fill=Ancestry)) +
    geom_col(position="stack") +
    scale_fill_manual(name="Ancestry",values = ancestry_colors)+
    theme_minimal() +
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1, size=12),
          strip.text.y = element_blank(),
          strip.text.x = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14, hjust = 0.5))  +
    ylab("Percent of LoF or ClinVar pathogenic carriers") + xlab("")
  p <- p + geom_text(aes(label=n), size=3, position=position_stack(vjust=0.5), colour="black") +
    facet_wrap(~ClinVarOrLOF, ncol=2, scales="free_x")
  
  pheno$Ancestry = factor(pheno$Ancestry, levels=ancestry_order)
  
  race3_table_full <- pheno %>%
    select("combined_id","Ancestry") %>%
    unique() %>%
    group_by(Ancestry) %>%
    dplyr::summarize(n = dplyr::n()) %>%
    mutate(Frequency=(n/sum(n))*100) %>%
    arrange(Ancestry) %>%
    filter(Ancestry != "Other")
  
  race3_table_full$Condition <- "            Full Cohort"
  p2 <- ggplot(data=race3_table_full, aes(x=Condition, y=Frequency, fill=Ancestry)) +
    geom_col(position="stack") +
    scale_fill_manual(name="Ancestry",values = ancestry_colors)+
    theme_minimal() +
    ylab("Percent of individuals in full cohort") + xlab("") + 
    ggtitle("") + 
    theme(strip.background = element_blank(), 
          strip.placement = "outside",
          axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1, size=12),
          strip.text.y = element_blank(),
          strip.text.x = element_text(size=14),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=14, hjust = 0.5),
          plot.title = element_text(size = 16)) 
  
  return(ggarrange(p, p2, ncol=2, common.legend = TRUE, legend="right", widths = c(1, 0.23)))
}



plot_carrier_phenotypes_one_condition_vus <- function(data,
                                                      AD_gene_list,
                                                      trait,
                                                      phenotype,
                                                      log_axis,
                                                      ylab,
                                                      title,
                                                      color_col = "Carrier"){
  pheno <- data$pheno
  full_carrier_list <- data$full_carrier_list
  genes <- subset(AD_gene_list, Trait==trait)$Gene_name
  carrier_subset <- annotate_phenotype_with_carriers(pheno, genes, full_carrier_list, trait,
                                                     c("combined_id", phenotype, "age"),
                                                     path_levels,
                                                     collapse_carrier_multiple = TRUE, stratify_vus = TRUE)
  carrier_subset <- carrier_subset %>%
    mutate(ClinVarOrLOF = replace(ClinVarOrLOF, Carrier=="Non Carrier","Non Carrier")) %>%
    filter(Trait == trait & !is.na(Carrier))
  
  carrier_subset <- carrier_subset %>%
    filter(Carrier %in% c("Non Carrier","Carrier","VUS Novel"))
  
  carrier_subset$Carrier = factor(carrier_subset$Carrier,levels = c("Non Carrier","Carrier","VUS Novel")) #,"VUS <1/10000","VUS <1/1000","VUS >=1/1000"))
  
  data_sources <- sort(unique(carrier_subset$data_source))
  carrier_subset$data_source = factor(carrier_subset$data_source, 
                                      levels = data_sources)
  p <- ggplot(carrier_subset, aes_string("Carrier", phenotype)) +
    geom_quasirandom(aes_string(colour=color_col), size=1, groupOnX=TRUE) + 
    stat_summary(fun.data=mean_cl_normal, size=0.8, geom="errorbar")  +
    stat_summary(fun.data=mean_cl_normal, size=3, geom="point")  +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=13),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_blank()) + 
    ylab(ylab) + xlab("") + labs(title=title)+
    facet_wrap(~data_source, nrow=1) + 
    scale_colour_manual(values=colors,
                        limits=subset(variant_levels, variant_levels %in% unique(as.character(carrier_subset %>% pull(color_col))))) +
    guides(colour = guide_legend(override.aes = list(size=3)))
  
  #a <- annotation_logticks(sides='l')
  #a$data <- data.frame(x=NA, data_source=data_sources[1])
  if(log_axis & phenotype == "bmi"){
    p <- p + scale_y_log10(breaks = c(10, 20, 30, 40, 50, 60, 70, 80)) #+ a
  }else if(log_axis){
    p <- p + scale_y_log10(breaks = c(10, 100, 1000, 5000)) #+ a
  }
  return(p)
}


variant_count_table <- function(full_carrier_list){
  full_carrier_list_clean <- full_carrier_list %>%
    mutate(ClinVarOrLOF = ifelse(ClinVarOrLOF %in% c("Review","ClinVar"), "ClinVar and Reviews", ClinVarOrLOF))
  
  variant_count_by_curation_type <- full_carrier_list_clean %>%
    select(data_source, ClinVarOrLOF,Variant_b37) %>% 
    unique() %>%
    group_by(data_source, ClinVarOrLOF) %>% 
    dplyr::summarise(number_variants = dplyr::n()) %>% 
    spread(ClinVarOrLOF, number_variants)
  
  variant_count_by_curation_type_total_curated <- full_carrier_list_clean %>%
    select(data_source, Variant_b37) %>% 
    unique() %>%
    group_by(data_source) %>% 
    dplyr::summarise(number_variants = dplyr::n())
  
  variant_count_by_curation_type_both_total_curated <- full_carrier_list_clean %>%
    select(Variant_b37) %>% 
    unique() %>%
    dplyr::summarise(number_variants = dplyr::n())
  
  variant_count_by_curation_type_total_curated <- rbind(variant_count_by_curation_type_total_curated,
                                                        c("Both",variant_count_by_curation_type_both_total_curated[[1]]))
  
  variant_count_by_curation_type_both <- full_carrier_list_clean %>%
    mutate(data_source = "Both") %>%
    select(data_source,ClinVarOrLOF,Variant_b37) %>% 
    unique() %>%
    group_by(data_source,ClinVarOrLOF) %>% 
    dplyr::summarise(number_variants = dplyr::n()) %>% 
    spread(ClinVarOrLOF, number_variants)
  
  variant_count_by_curation_all <- rbind(variant_count_by_curation_type, variant_count_by_curation_type_both)
  variant_count_by_curation_all[,"Total Curated"] <- variant_count_by_curation_type_total_curated$number_variants
  return(variant_count_by_curation_all)
}



carrier_count_table_by_condition <- function(pheno, full_carrier_list){
  total_n <- pheno %>%
    group_by(data_source) %>% 
    dplyr::summarise(count = dplyr::n()) %>% 
    spread(data_source, count)
  
  full_carrier_list_clean <- full_carrier_list %>%
    filter(!((Trait == "MODY Extended")&(Gene %in% c("GCK","HNF1A","HNF1B","PDX1","HNF4A")))) %>%
    mutate(ClinVarOrLOF = ifelse(ClinVarOrLOF %in% c("Review","ClinVar"), "ClinVar and Reviews", ClinVarOrLOF))
  
  variant_count_by_trait <- full_carrier_list_clean %>%
    select(data_source,Trait,Variant) %>% 
    unique() %>%
    group_by(data_source,Trait) %>% 
    dplyr::summarise(number_variants = dplyr::n())
  
  carrier_count_by_trait <- full_carrier_list_clean %>%
    select(SubjectID,data_source,Trait) %>% 
    unique() %>%
    group_by(data_source,Trait) %>% 
    dplyr::summarise(number_carriers = dplyr::n()) %>% 
    mutate(percent_carriers = (number_carriers/total_n[data_source][[1]])*100) 
  
  carrier_count_by_trait <- variant_count_by_trait %>%
    full_join(carrier_count_by_trait) %>%
    filter(Trait %in% c("High LDL",
                        "Low LDL",
                        "High HDL",
                        "High triglycerides",
                        "Obesity",
                        "MODY",
                        "MODY Extended",
                        "Lipodystrophy",
                        "Neonatal Diabetes"))
  
  carrier_count_by_trait$Trait = factor(carrier_count_by_trait$Trait, 
                                                  levels = c("High LDL",
                                                             "Low LDL",
                                                             "High HDL",
                                                             "High triglycerides",
                                                             "Obesity",
                                                             "MODY",
                                                             "MODY Extended",
                                                             "Lipodystrophy",
                                                             "Neonatal Diabetes"))
  carrier_count_by_trait <- carrier_count_by_trait %>% arrange(Trait)
  
  return(carrier_count_by_trait )
}


formatted_penetrance <- function(all_penetrance, p_values, by_gene=FALSE, carrier_type="pathogenic"){
  print("pval subset1")
  print(p_values)
  p_values_subset <- p_values %>%
    filter((type == carrier_type) & (additional_covariates %in% c(NA,"","covart2d,covarlipidmeds","covart2d"))) %>% 
    mutate(Condition = gsub("_", ' ', Condition)) %>%
    mutate(lipidmeds = ifelse(additional_covariates %in% c(NA,"","covart2d"),FALSE,TRUE),
           "P value" = formatC(PVALUE, digits = 3),
           "Beta (se)" = paste0(formatC(BETA, digits = 3), " (", formatC(SEBETA, digits = 3), ")"),
           "OR" = paste0(formatC(exp(BETA), digits = 3, width = 0), " (", formatC(exp(BETA - SEBETA), digits = 3, width = 0), "-", formatC(exp(BETA + SEBETA), digits = 3, width = 0), ")"),
           Condition = case_when(Condition == "High triglycerides" ~ "High TG",
                                 Condition == "High triglycerides ADJ" ~ "High TG ADJ",
                                 Condition == "Obesity" ~ "Monogenic obesity",
                                 TRUE ~ Condition))
  print("pval subset2")
  print(p_values_subset)
  if(by_gene){
    p_values_subset <- p_values_subset %>%
      select(Condition,Trait,Gene,data_source,lipidmeds,"P value","Beta (se)",OR)
    join_by <- c("condition_simple"="Condition", "phenotype"="Trait", "data_source"="data_source", "lipidmeds"="lipidmeds", "Gene"="Gene")
  }else{
    p_values_subset <- p_values_subset %>%
      select(Condition,Trait,data_source,lipidmeds,"P value","Beta (se)",OR)
    join_by <- c("condition_simple"="Condition", "phenotype"="Trait", "data_source"="data_source", "lipidmeds"="lipidmeds")
    
  }
  print("pval subset")
  print(p_values_subset)
  all_penetrance_format <- all_penetrance %>%
    mutate(condition_simple = as.character(condition_simple)) %>%
    mutate(phenotype = case_when(
      condition_simple == "Low LDL" & (phenotype == "ldl_mgdl_adj") ~ "low_ldl_adj",
      condition_simple == "High LDL" & (phenotype == "ldl_mgdl_adj") ~ "high_ldl_adj",
      condition_simple == "High TG" & (phenotype == "tg_mgdl_adj") ~ "high_tg_adj",
      condition_simple == "Low LDL" ~ "low_ldl",
      condition_simple == "High LDL" ~ "high_ldl",
      condition_simple == "High HDL" ~ "high_hdl",
      condition_simple == "High TG" ~ "high_tg",
      condition_simple == "Monogenic obesity" ~ "Obese",
      TRUE ~ phenotype),
      prop = paste0(formatC(prop, digits = 3), " (", formatC(low, digits = 3), "-", formatC(high, digits = 3), ")")) %>%
    select(-diabetes_definition,-high,-low) 
  print(all_penetrance_format)
  if(by_gene){
    all_penetrance_format <- all_penetrance_format %>%
      left_join(p_values_subset,by=join_by) %>%
      pivot_wider(names_from = c("data_source"),values_from=c("count","total","prop","Beta (se)","P value","OR")) %>%
      mutate("Restricted no lipid medications" = ifelse(lipidmeds,"Yes","")) %>%
      rename(Condition=condition) %>%
    select("Condition","Gene","total_AMP-T2D-GENES", "prop_AMP-T2D-GENES","Beta (se)_AMP-T2D-GENES","P value_AMP-T2D-GENES","OR_AMP-T2D-GENES",
           "total_UKBB", "prop_UKBB","Beta (se)_UKBB","P value_UKBB","OR_UKBB")
    print(all_penetrance_format$`P value_UKBB`)
  }else{
    all_penetrance_format <- all_penetrance_format %>%
      pivot_wider(names_from = c("Carrier"),values_from=c("count","total","prop")) %>% 
      left_join(p_values_subset,by=join_by) %>%
      pivot_wider(names_from = c("data_source"),values_from=c("count_Carrier","count_Non Carrier","total_Carrier","total_Non Carrier","prop_Carrier","prop_Non Carrier","Beta (se)","P value","OR")) %>%
      mutate("Restricted no lipid medications" = ifelse(lipidmeds,"Yes","")) %>%
      rename(Condition=condition) %>%
      select("Condition","Restricted no lipid medications","total_Non Carrier_AMP-T2D-GENES","count_Non Carrier_AMP-T2D-GENES", "prop_Non Carrier_AMP-T2D-GENES",
             "total_Carrier_AMP-T2D-GENES","count_Carrier_AMP-T2D-GENES","prop_Carrier_AMP-T2D-GENES","Beta (se)_AMP-T2D-GENES","P value_AMP-T2D-GENES","OR_AMP-T2D-GENES",
             "total_Non Carrier_UKBB","count_Non Carrier_UKBB", "prop_Non Carrier_UKBB",
             "total_Carrier_UKBB","count_Carrier_UKBB","prop_Carrier_UKBB","Beta (se)_UKBB","P value_UKBB","OR_UKBB")
  }
  
  all_penetrance_format$Condition = factor(all_penetrance_format$Condition, 
                                           levels = c("High LDL (>190 mg/dl)",
                                                      "High LDL (>190 mg/dl) no lipid meds",
                                                      "High LDL ADJ (>190 mg/dl)",
                                                      "Low LDL (<80 mg/dl)",
                                                      "Low LDL (<80 mg/dl) no lipid meds",
                                                      "Low LDL ADJ (<80 mg/dl)",
                                                      "High HDL (>70 mg/dl)",
                                                      "High triglycerides (>200 mg/dl)",
                                                      "High triglycerides (>200 mg/dl) no lipid meds",
                                                      "High triglycerides ADJ (>200 mg/dl)",
                                                      "Monogenic obesity (BMI >30)",
                                                      "MODY (T2D)",
                                                      "MODY (T2D and pre-diabetes)",
                                                      "MODY extended (T2D)",
                                                      "MODY extended (T2D and pre-diabetes)",
                                                      "Lipodystrophy (T2D)",
                                                      "Lipodystrophy (T2D and pre-diabetes)",
                                                      "Neonatal diabetes (T2D)",
                                                      "Neonatal diabetes (T2D and pre-diabetes)"))
  
  all_penetrance_format <- all_penetrance_format %>% arrange(Condition)
  
  return(all_penetrance_format)
}


final_variant_table <- function(full_carrier_list){
  other_notes <- read.delim("../data/MODY_variant_independent_assesment.txt", header=T, sep="\t", stringsAsFactors=F)
  other_notes <- other_notes %>% select("Variant_b37", "Variant_type", "Assessment_Notes", "Additional_independent_assessment", "Additional_independent_assessment_notes")
  full_variant_list_clean <- full_carrier_list %>%
    mutate(Level=factor(tolower(Level),levels=all_levels)) %>%
    group_by(Variant_b37,Trait,ClinVarOrLOF,Gene,VariationID,AlleleID,LastEvaluated,high_confidence_lab_2017_clinical_significance,Level) %>% 
    mutate(Notes = paste(sort(unique(Notes)),collapse=", ")) %>% 
    ungroup() %>% 
    group_by(Variant_b37,Trait,ClinVarOrLOF,Gene,VariationID,AlleleID,LastEvaluated,high_confidence_lab_2017_clinical_significance,Level,data_source,Notes) %>% 
    dplyr::summarise(N=dplyr::n()) %>% 
    spread(data_source,N) %>%
    rename(Variant_type = ClinVarOrLOF) %>%
    arrange(Trait,Gene,Level,Variant_b37) %>%
    filter((Trait != "MonogenicDM_all") & !((Trait == "MODY Extended")&(Gene %in% c("GCK","HNF1A","HNF1B","PDX1","HNF4A")))) %>%
    left_join(other_notes,by=c("Variant_b37"="Variant_b37","Variant_type"="Variant_type")) %>%
    left_join(full_carrier_list %>% filter(data_source == "UKBB") %>%select("Variant_b37","Variant") %>% unique(), by = c("Variant_b37"="Variant_b37")) 
    

  return(full_variant_list_clean)
}






















penetrance_threshold <- function(data,
                                 trait,
                                 phenotype,
                                 label,
                                 threshold_under=FALSE,
                                 lipidmeds=FALSE,
                                 stratify_vus = FALSE,
                                 levels_to_plot=path_levels,
                                 by_gene = FALSE){
  pheno <- data$pheno
  full_carrier_list <- data$full_carrier_list
  phenotype_var <- enquo(phenotype)
  phenotype_str <- deparse(substitute(phenotype))
  pheno_reduced_ann <- annotate_phenotype_with_carriers(pheno,
                                                        subset(AD_gene_list, Trait==trait)$Gene_name,
                                                        full_carrier_list,
                                                        trait,
                                                        c("combined_id",phenotype_str,"lipidmeds","esp_phenotype"),
                                                        levels_to_plot, collapse_carrier_multiple = TRUE,
                                                        stratify_vus=stratify_vus,
                                                        by_gene=by_gene)

  if(phenotype_str %in% c("t2d","t2d_with_pre")){
    threshold_results <- pheno_reduced_ann %>%
      mutate(above_threshold = (!!phenotype_var == 1))
  }else{
    threshold <- threshold_for_plot$threshold[which(threshold_for_plot$Trait==trait)]
    if(threshold_under){
      threshold_results <- pheno_reduced_ann %>%
        mutate(above_threshold = (!!phenotype_var < threshold))
    }else{
      threshold_results <- pheno_reduced_ann %>%
        mutate(above_threshold = (!!phenotype_var > threshold))
    }
  }
  if(lipidmeds){
    threshold_results <- threshold_results %>%
      filter(lipidmeds == 0)
  }
  
  if(by_gene){
    threshold_results <- threshold_results %>%
      mutate(above_threshold = factor(above_threshold),Gene = factor(Gene)) %>%
      group_by(data_source, Carrier, Gene, above_threshold, .drop=FALSE) %>% 
      dplyr::summarize(count = dplyr::n()) %>% 
      filter(!is.na(above_threshold)) %>% 
      group_by(data_source, Carrier, Gene) %>% 
      mutate(prop = count/sum(count), total = sum(count)) %>% 
    filter((above_threshold=="TRUE") & !is.nan(prop)) #%>%
  }else{
    threshold_results <- threshold_results %>%
      mutate(above_threshold = factor(above_threshold)) %>%
      group_by(data_source, Carrier, above_threshold, .drop=FALSE) %>% 
      dplyr::summarize(count = dplyr::n()) %>% 
      filter(!is.na(above_threshold)) %>% 
      group_by(data_source, Carrier) %>% 
      mutate(prop = count/sum(count), total = sum(count)) %>% 
    filter(above_threshold=="TRUE" & !is.nan(prop)) #%>%
  }
  
  threshold_results$condition = label
  threshold_results$condition_simple = trait
  threshold_results$phenotype = phenotype_str
  if(phenotype_str == "t2d"){
    threshold_results$diabetes_definition = "T2D"
  }else if(phenotype_str == "t2d_with_pre"){
    threshold_results$diabetes_definition = "T2D and pre-diabetes" 
  }else{
    threshold_results$diabetes_definition = NA
  }
  
  if(lipidmeds){
    threshold_results <- threshold_results %>% mutate(condition = paste(condition, "no lipid meds"),
                                                      lipidmeds = TRUE)
  }else{
    threshold_results <- threshold_results %>% mutate(lipidmeds = FALSE)
  }
  
  return(threshold_results)
}


age_of_diagnosis_survival_plot <- function(data, AD_gene_list){
  carrier_subset <- annotate_phenotype_with_carriers(data$pheno,
                                                     subset(AD_gene_list, Trait=="MODY")$Gene_name,
                                                     data$full_carrier_list, "MODY",
                                                     c("combined_id","age","t2d","aod","Age_diabetes_diagnosed_0","Age_diabetes_diagnosed_1","Age_diabetes_diagnosed_2","lipidmeds"),
                                                     path_levels,
                                                     collapse_carrier_multiple = TRUE,
                                                     by_gene = TRUE)
  
  
  bla <- carrier_subset %>%
    filter(Gene %in% c("GCK","HNF1A","Non Carrier")) %>%
    mutate(Gene = factor(Gene,levels = c("GCK","HNF1A","Non Carrier"))) %>%
    select(Carrier,age,t2d,aod,Age_diabetes_diagnosed_0,Age_diabetes_diagnosed_1,Age_diabetes_diagnosed_2,Gene,data_source) %>%
    mutate(aod = case_when(!is.na(aod) ~ as.integer(aod),
                           !is.na(Age_diabetes_diagnosed_0) ~ Age_diabetes_diagnosed_0,
                           !is.na(Age_diabetes_diagnosed_1) ~ Age_diabetes_diagnosed_1,
                           !is.na(Age_diabetes_diagnosed_2) ~ Age_diabetes_diagnosed_2,
                           t2d == 0 ~ as.integer(age)))
  
  
  ukbb_aod_surv <- survfit(Surv(aod, t2d) ~Gene, data = bla %>% filter(data_source == "UKBB"))
  ampt2d_aod_surv <- survfit(Surv(aod, t2d) ~Gene, data = bla %>% filter(data_source == "AMP-T2D-GENES"))

  ukbb_p <- ggsurvplot(ukbb_aod_surv,
                       data = bla,
                       surv.median.line="none",
                       legend.title = "",
                       linetype = c("solid","dashed","solid"),
                       size=1.2,
                       palette = c("#d6604d","#d6604d","#666666"),
                       censor = FALSE,
                       legend.labs = c("GCK","HNF1A","Non Carrier"),
                       conf.int = FALSE)
  
  ukbb_p <- ukbb_p$plot +
    theme(legend.key.width = unit(3.5, "line")) +
    theme(legend.position = "bottom")+
    ggtitle("UKBB") + 
    ylab("Proportion without diabetes") + xlab("Age (years)") +
    theme(plot.title = element_text(size=14,hjust = 0.5,family = "Helvetica"),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13),
          axis.title.y = element_text(size=13),
          axis.title.x = element_text(size=13),
          legend.text=element_text(size=13),
          legend.title=element_blank())

  ampt2d_p <- ggsurvplot(ampt2d_aod_surv,
                         data = bla,
                         surv.median.line="none",
                         legend.title = "",
                         linetype = c("solid","dashed","solid"),
                         size=1.2,
                         palette = c("#d6604d","#d6604d","#666666"),
                         censor = FALSE,
                         legend.labs = c("GCK","HNF1A","Non Carrier"),
                         conf.int = FALSE)
  
  ampt2d_p <- ampt2d_p$plot +
    theme(legend.key.width = unit(3.5, "line")) +
    theme(legend.position = "none") +
    ggtitle("AMP-T2D-GENES") + 
    ylab("Proportion without diabetes") + xlab("Age (years)") +
    theme(plot.title = element_text(size=14,hjust = 0.5,family = "Helvetica"),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13),
          axis.title.y = element_text(size=13),
          axis.title.x = element_text(size=13),
          legend.text=element_text(size=13),
          legend.title=element_blank()) 
  
  return(ggarrange(ampt2d_p, ukbb_p, heights = c(1,1), nrow=2,labels = c("A","B")))
  
}

participant_characteristics <- function(data){
  # Creates a table of 3 columns: AMP cases, AMP controls, UKB
  # Each row contains basic study characteristics: 
  # %female, mean age (age range), mean BMI (SD), LDL cholesterol mean (SD), 
  # HDL cholesterol mean (SD), triglycerides mean (SD), 
  # ancestry breakdown with N white, african american, east Asian, south Asian, hispanic/latino
  data <- data %>% select("sex", "age", "bmi", "ldl_mgdl", "ldl_mgdl_adj", "hdl_mgdl", "tg_mgdl", "tg_mgdl_adj", "Ancestry", "data_source", "t2d")
  data <- data %>% mutate(grouping=case_when(data$data_source == "UKBB" ~ "UKB", data$t2d == 0 ~ "AMP-T2D control", data$t2d == 1 ~ "AMP-T2D case"),
                          Ancestry=case_when(data$data_source == "UKBB" ~ "European", data$data_source != "UKBB" ~ Ancestry),
                          sex=if_else(sex==1, "Male", "Female"))
  sex <- as.data.frame(data %>% group_by(grouping) %>% count(sex) %>% dcast(grouping ~ sex, value.var = "n") %>% t() %>% row_to_names(row_number = 1)) %>% rownames_to_column()
  means_sds <- as.data.frame(data %>% group_by(grouping) %>% summarise(across(c("age", "bmi", "ldl_mgdl", "ldl_mgdl_adj", "hdl_mgdl", "tg_mgdl", "tg_mgdl_adj"), list(mean=~mean(.x, na.rm = TRUE), sd=~sd(.x, na.rm = TRUE)))) %>% t() %>% row_to_names(row_number = 1)) %>% rownames_to_column()
  ancestry <- as.data.frame(data %>% group_by(grouping) %>% count(Ancestry) %>% dcast(grouping ~ Ancestry, value.var = "n") %>% t() %>% row_to_names(row_number = 1)) %>% rownames_to_column()
  
  
  all_data = sex %>% rbind(means_sds) %>% rbind(ancestry)
  
  return(list(sex=sex, means_sds=means_sds, ancestry=ancestry, all = all_data))
}


ascertainment_plots <- function(all_data){
  epacts_unrelated <- read.table("../data/Phenotypes/t2d.epacts.covar.unrelated.ped.txt", header=T, sep = "\t", stringsAsFactors=F, comment.char = "", check.names = FALSE)
  
  pheno <- all_data$pheno_annotate_carrier_more_severe
  pheno <- pheno %>% inner_join(epacts_unrelated %>% select(-t2d), by=c("combined_id"="IID"))
  
  # Get means for individuals acertained on High LDL cholesterol and individuals not ascertained on High LDL cholesterol
  # and a comparison of those means
  pheno_sub <- pheno %>%
    mutate(in_excluded_esp_cohort=esp_phenotype %in% c("LDL_High")) %>% 
    filter(!(esp_phenotype %in% c("BMI_High","BMI_High_Nondiab","BMI_High_diab","BMI_Low","LDL_Low"))) %>%
    select(age, sex, C1, C2, C3, C4, C5,C6,C7,C8,C9,C10,combined_id,Trait,in_excluded_esp_cohort, ldl_mgdl, ldl_mgdl_adj, bmi, esp_phenotype,lipidmeds) %>%
    filter(!is.na(ldl_mgdl)) %>%
    distinct()
  
  high_ldl_pval <- summary(
    lm(
      ldl_mgdl_adj~in_excluded_esp_cohort + age + sex + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, 
      pheno_sub %>% filter(Trait == "High LDL")
      )
    )$coefficients[2,4] 
  
  high_ldl_stats <- pheno_sub %>% filter(Trait == "High LDL") %>% group_by(in_excluded_esp_cohort) %>%
    dplyr::summarize(mn=mean(ldl_mgdl_adj), n=n(),sd=sd(ldl_mgdl_adj)) %>% 
    mutate(se=sd/sqrt(n),LCI=mn-qt(1 - ((1 - 0.95) / 2), n - 1) * se,UCI=mn+qt(1 - ((1 - 0.95) / 2), n - 1) * se)
  
  
  
  pheno_sub <- pheno %>%
    mutate(in_excluded_esp_cohort=esp_phenotype %in% c("LDL_Low")) %>%
    filter(!(esp_phenotype %in% c("BMI_High","BMI_High_Nondiab","BMI_High_diab","BMI_Low","LDL_High"))) %>% 
    select(age, sex, C1, C2, C3, C4, C5,C6,C7,C8,C9,C10,combined_id,Trait,in_excluded_esp_cohort, ldl_mgdl, ldl_mgdl_adj, bmi, esp_phenotype,lipidmeds) %>%
    filter(!is.na(ldl_mgdl)) %>%
    distinct()

  low_ldl_pval <- summary(
    lm(
      ldl_mgdl_adj~in_excluded_esp_cohort + age + sex + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, 
      pheno_sub %>% filter(Trait == "Low LDL")
    )
  )$coefficients[2,4] 
  
  low_ldl_stats <- pheno_sub %>% filter(Trait == "Low LDL") %>%
    group_by(in_excluded_esp_cohort) %>% 
    dplyr::summarize(mn=mean(ldl_mgdl_adj), n=n(),sd=sd(ldl_mgdl_adj)) %>% 
    mutate(se=sd/sqrt(n),LCI=mn-qt(1 - ((1 - 0.95) / 2), n - 1) * se,UCI=mn+qt(1 - ((1 - 0.95) / 2), n - 1) * se)
  
  
  
  pheno <- all_data$pheno_annotate_carrier_more_severe %>%
    mutate(in_excluded_esp_cohort=esp_phenotype %in% c("LDL_High")) %>% 
    filter(!(esp_phenotype %in% c("BMI_High","BMI_High_Nondiab","BMI_High_diab","BMI_Low","LDL_Low")))
  pheno <- pheno %>% inner_join(epacts_unrelated %>% select(-t2d), by=c("combined_id"="IID"))
  pheno_sub <- pheno %>%
    filter(Trait == "High LDL") %>% 
    select(age, sex, C1, C2, C3, C4, C5,C6,C7,C8,C9,C10,combined_id,Trait,in_excluded_esp_cohort, ldl_mgdl, ldl_mgdl_adj, bmi, esp_phenotype,lipidmeds,Variant) %>% 
    filter(!is.na(ldl_mgdl)) %>% 
    distinct()
  pheno_non_carrier <- all_data$pheno %>%
    filter(!(combined_id %in% pheno_sub$combined_id) & !(esp_phenotype %in% c("BMI_High","BMI_High_Nondiab","BMI_High_diab","BMI_Low","LDL_Low","LDL_High"))) %>% 
    inner_join(epacts_unrelated %>% 
                 select(-t2d), by=c("combined_id"="IID")
    ) %>% 
    mutate(Trait="High LDL", in_excluded_esp_cohort="Non-carriers", Variant=NA) %>% 
    select(age, sex, C1, C2, C3, C4, C5,C6,C7,C8,C9,C10,combined_id,Trait,in_excluded_esp_cohort, ldl_mgdl, ldl_mgdl_adj, bmi, esp_phenotype,lipidmeds,Variant) %>% 
    filter(!is.na(ldl_mgdl))
  pheno_sub_w_non_carrier <- pheno_sub %>% mutate(in_excluded_esp_cohort=ifelse(in_excluded_esp_cohort, "Carriers,\nascertained", "Carriers,\nnot ascertained")) %>% bind_rows(pheno_non_carrier)
  pheno_sub_w_non_carrier <- pheno_sub_w_non_carrier %>% mutate(in_excluded_esp_cohort=factor(in_excluded_esp_cohort, levels = c("Non-carriers", "Carriers,\nnot ascertained", "Carriers,\nascertained")))
  

  p3 <- ggplot(pheno_sub_w_non_carrier %>% filter(Trait == "High LDL"), aes(x=in_excluded_esp_cohort, y=ldl_mgdl)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=3,aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.position = "none",
          legend.title.align=0.5) +
    xlab("") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("black","royal blue2","red 2"), name="Ascertained on\nhigh LDL")
  
  p3_adj <- ggplot(pheno_sub_w_non_carrier %>% filter(Trait == "High LDL"), aes(x=in_excluded_esp_cohort, y=ldl_mgdl_adj)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=3,aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.position = "none",
          legend.title.align=0.5) +
    xlab("") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("#666666", "royal blue2","red 2"), name="Ascertained on\nhigh LDL")
  
  p3_adj_oth <- ggplot(pheno_sub_w_non_carrier %>% filter(Trait == "High LDL"), aes(x=in_excluded_esp_cohort, y=ldl_mgdl_adj)) +
    geom_quasirandom(aes(color=in_excluded_esp_cohort),size=3, groupOnX=TRUE) + 
    #stat_summary(fun.data=mean_cl_normal, size=0.8,geom="errorbar")  +
    stat_summary(fun.data=mean_cl_normal, size=1,geom="errorbar")  +
    stat_summary(fun.data=mean_cl_normal, size=3, geom="point")  +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.position = "none",
          legend.title.align=0.5) +
    xlab("") + ylab("LDL cholesterol (mg/dL)") +
    scale_color_manual(values=c("#666666", "royal blue2","red 2"), name="Ascertained on\nhigh LDL")
  
  
  
  variants_in_esp <- (pheno_sub %>% filter(in_excluded_esp_cohort))$Variant
  variants_no_esp <- (pheno_sub %>% filter(!in_excluded_esp_cohort))$Variant
  pheno_sub_in_both <- pheno_sub %>% filter((Variant %in% variants_in_esp) & (Variant %in% variants_no_esp))
  high_ldl_share_pval <- summary(lm(ldl_mgdl_adj~in_excluded_esp_cohort + age + sex + C1 + C2 + C3 + C4 + C5, pheno_sub_in_both %>% filter(Trait == "High LDL")))$coefficients[2, 4]
  high_ldl_share_stats <- pheno_sub_in_both %>% group_by(in_excluded_esp_cohort) %>% dplyr::summarize(mn=mean(ldl_mgdl_adj), n=n(), sd=sd(ldl_mgdl_adj)) %>% mutate(se=sd/sqrt(n), LCI=mn-qt(1 - ((1 - 0.95) / 2), n - 1) * se, UCI=mn+qt(1 - ((1 - 0.95) / 2), n - 1) * se)
  
  pheno_sub_in_both <- pheno_sub_in_both %>% mutate("Variant_In_ESP_Cohort"=paste(Variant, in_excluded_esp_cohort, sep="-"))
  pheno_sub_in_both <- pheno_sub_in_both %>% mutate(Variant = case_when(Variant=="19-11213450-G-A"~"LDLR p.Glu101Lys",
                                                                        Variant=="19-11217344-T-A"~"LDLR p.Asp266Glu",
                                                                        Variant=="19-11227604-G-A"~"LDLR p.Gly592Glu",
                                                                        Variant=="2-21229160-C-T"~"APOB p.Arg3527Gln"))
  
  
  p1 <- ggplot(pheno_sub_in_both %>% filter(Trait == "High LDL"), aes(x=Variant, y=ldl_mgdl)) + geom_jitter(size=3, aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title.align=0.5) +
    xlab("Shared LDL-raising variant") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("royal blue2","red 2"),name="Ascertained on\nhigh LDL",
                        breaks=c(FALSE,TRUE),
                        labels=c("No","Yes"))
  
  p1_adj <- ggplot(pheno_sub_in_both %>% filter(Trait == "High LDL"), aes(x=Variant, y=ldl_mgdl_adj)) + geom_jitter(size=3, aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title.align=0.5) +
    xlab("Shared LDL-raising variant") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("royal blue2","red 2"),name="Ascertained on\nhigh LDL",
                        breaks=c(FALSE,TRUE),
                        labels=c("No","Yes"))
  
  
  pheno <- all_data$pheno_annotate_carrier_more_severe %>% mutate(in_excluded_esp_cohort=esp_phenotype %in% c("LDL_Low")) %>% filter(!(esp_phenotype %in% c("BMI_High", "BMI_High_Nondiab", "BMI_High_diab", "BMI_Low", "LDL_High")))
  pheno <- pheno %>% filter(Trait %in% c("Low LDL"))

  pheno <- pheno %>% inner_join(epacts_unrelated %>% select(-t2d), by=c("combined_id"="IID"))
  pheno_sub <- pheno %>% filter(Trait == "Low LDL") %>% select(age, sex, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, combined_id, Trait, in_excluded_esp_cohort, ldl_mgdl, ldl_mgdl_adj, bmi, esp_phenotype, lipidmeds, Variant) %>% filter(!is.na(ldl_mgdl)) %>% distinct()
  pheno_non_carrier <- all_data$pheno %>% filter(!(combined_id %in% pheno_sub$combined_id) & !(esp_phenotype %in% c("BMI_High", "BMI_High_Nondiab", "BMI_High_diab", "BMI_Low", "LDL_Low", "LDL_High"))) %>% inner_join(epacts_unrelated %>% select(-t2d), by=c("combined_id"="IID")) %>% mutate(Trait="Low LDL", in_excluded_esp_cohort="Non-carriers", Variant=NA) %>% select(age, sex, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, combined_id, Trait, in_excluded_esp_cohort, ldl_mgdl, ldl_mgdl_adj, bmi, esp_phenotype, lipidmeds, Variant) %>% filter(!is.na(ldl_mgdl))
  pheno_sub_w_non_carrier <- pheno_sub %>% mutate(in_excluded_esp_cohort=ifelse(in_excluded_esp_cohort, "Carriers,\nascertained", "Carriers,\nnot ascertained")) %>% bind_rows(pheno_non_carrier)
  pheno_sub_w_non_carrier <- pheno_sub_w_non_carrier %>% mutate(in_excluded_esp_cohort=factor(in_excluded_esp_cohort, levels = c("Non-carriers", "Carriers,\nnot ascertained", "Carriers,\nascertained")))
  
  p4 <- ggplot(pheno_sub %>% filter(Trait == "Low LDL") %>% mutate(in_excluded_esp_cohort=ifelse(in_excluded_esp_cohort, "Yes", "No")), aes(x=in_excluded_esp_cohort, y=ldl_mgdl)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=3,aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.position = "none",
          legend.title.align=0.5) +
    xlab("") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("royal blue2","red 2"),name="Ascertained on\nlow LDL")
  
  p4_adj <- ggplot(pheno_sub %>% filter(Trait == "Low LDL") %>% mutate(in_excluded_esp_cohort=ifelse(in_excluded_esp_cohort, "Yes", "No")), aes(x=in_excluded_esp_cohort, y=ldl_mgdl_adj)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(size=3,aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.position = "none",
          legend.title.align=0.5) +
    xlab("") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("royal blue2","red 2"),name="Ascertained on\nlow LDL")
  
  p4_adj_oth <- ggplot(pheno_sub_w_non_carrier %>% filter(Trait == "Low LDL"), aes(x=in_excluded_esp_cohort, y=ldl_mgdl_adj)) +
    geom_quasirandom(aes(color=in_excluded_esp_cohort),size=3, groupOnX=TRUE) + 
    stat_summary(fun.data=mean_cl_normal, size=1,geom="errorbar")  +
    stat_summary(fun.data=mean_cl_normal, size=3, geom="point")  +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.position = "none",
          legend.title.align=0.5) +
    xlab("") + ylab("LDL cholesterol (mg/dL)") +
    scale_color_manual(values=c("#666666", "royal blue2","red 2"), name="Ascertained on\nlow LDL")
  
  variants_in_esp <- (pheno_sub %>% filter(in_excluded_esp_cohort))$Variant
  variants_no_esp <- (pheno_sub %>% filter(!in_excluded_esp_cohort))$Variant
  pheno_sub_in_both <- pheno_sub %>% filter((Variant %in% variants_in_esp) & (Variant %in% variants_no_esp))
  
  low_ldl_share_pval <- summary(lm(ldl_mgdl_adj~in_excluded_esp_cohort + age + sex + C1 + C2 + C3 + C4 + C5, pheno_sub_in_both %>% filter(Trait == "Low LDL")))$coefficients[2, 4]
  low_ldl_share_stats <- pheno_sub_in_both %>% group_by(in_excluded_esp_cohort) %>% dplyr::summarize(mn=mean(ldl_mgdl_adj), n=n(), sd=sd(ldl_mgdl_adj)) %>% mutate(se=sd/sqrt(n), LCI=mn-qt(1 - ((1 - 0.95) / 2), n - 1) * se, UCI=mn+qt(1 - ((1 - 0.95) / 2), n - 1) * se)
  
  
  pheno_sub_in_both <- pheno_sub_in_both %>% mutate("Variant_In_ESP_Cohort"=paste(Variant, in_excluded_esp_cohort, sep="-"))
  pheno_sub_in_both <- pheno_sub_in_both %>% mutate(Variant = case_when(Variant=="1-55512222-C-G"~"PCSK9 p.Tyr142Ter"))
  p2 <- ggplot(pheno_sub_in_both %>% filter(Trait == "Low LDL"), aes(Variant, ldl_mgdl)) +
    geom_jitter(size=3,aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title.align=0.5) +
    xlab("Shared LDL-lowering variant") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("royal blue2","red 2"),name="Ascertained on\nlow LDL",
                        breaks=c(FALSE,TRUE),
                        labels=c("No","Yes"))
  
  p2_adj <- ggplot(pheno_sub_in_both %>% filter(Trait == "Low LDL"), aes(Variant, ldl_mgdl_adj)) +
    geom_jitter(size=3,aes(color=in_excluded_esp_cohort), width=0.2) +
    theme_classic() + 
    theme(plot.title = element_text(size=13,hjust=0.5),
          strip.background = element_blank(), 
          strip.placement = "outside",
          strip.text.x = element_text(size=10),
          strip.text.y = element_blank(),
          axis.text.y = element_text(size=12),
          axis.text.x = element_text(size=12),
          axis.title.y = element_text(size=12),
          axis.title.x = element_text(size=12),
          legend.text=element_text(size=12),
          legend.title.align=0.5) +
    xlab("Shared LDL-lowering variant") + ylab("LDL cholesterol (mg/dL)") +
    scale_colour_manual(values=c("royal blue2","red 2"),name="Ascertained on\nlow LDL",
                        breaks=c(FALSE,TRUE),
                        labels=c("No","Yes"))
  
  ggarrange(p3,p1,p4,p2,ncol=2,nrow=2,widths=c(1,2))
  ggarrange(p3_adj,p1_adj,p4_adj,p2_adj,ncol=2,nrow=2,widths=c(1,2))
  p_all_adj <- ggarrange(p3_adj_oth, p1_adj, p4_adj_oth, p2_adj, ncol=2, nrow=2, widths=c(1, 2))
  
  high_ldl_stats <- high_ldl_stats %>%
    mutate(in_excluded_esp_cohort=if_else(in_excluded_esp_cohort, "LDL_ascertained" , "not_LDL_ascertained")) %>% 
    pivot_wider(names_from = in_excluded_esp_cohort, values_from = c(mn, n, sd, se, LCI, UCI)) %>%
    add_column(Pvalue=high_ldl_pval) %>%
    add_column(Condition="High LDL")
  
  low_ldl_stats <- low_ldl_stats %>%
    mutate(in_excluded_esp_cohort=if_else(in_excluded_esp_cohort, "LDL_ascertained" , "not_LDL_ascertained")) %>% 
    pivot_wider(names_from = in_excluded_esp_cohort, values_from = c(mn, n, sd, se, LCI, UCI)) %>%
    add_column(Pvalue=low_ldl_pval) %>%
    add_column(Condition="Low LDL")
  
  high_ldl_share_stats <- high_ldl_share_stats %>%
    mutate(in_excluded_esp_cohort=if_else(in_excluded_esp_cohort, "LDL_ascertained" , "not_LDL_ascertained")) %>% 
    pivot_wider(names_from = in_excluded_esp_cohort, values_from = c(mn, n, sd, se, LCI, UCI)) %>%
    add_column(Pvalue=high_ldl_share_pval) %>%
    add_column(Condition="High LDL - Restricted to shared variants")
  
  low_ldl_share_stats <- low_ldl_share_stats %>%
    mutate(in_excluded_esp_cohort=if_else(in_excluded_esp_cohort, "LDL_ascertained" , "not_LDL_ascertained")) %>% 
    pivot_wider(names_from = in_excluded_esp_cohort, values_from = c(mn, n, sd, se, LCI, UCI)) %>%
    add_column(Pvalue=low_ldl_share_pval) %>%
    add_column(Condition="Low LDL - Restricted to shared variants")
  
  ascertainment_stats <- high_ldl_stats %>% union(low_ldl_stats) %>% union(high_ldl_share_stats) %>% union(low_ldl_share_stats)
  
  
  return(
    list(
      ascertainment_plot=p_all_adj, 
      ascertainment_stats=ascertainment_stats
      )
    )
}
