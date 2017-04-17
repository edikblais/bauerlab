
options(java.parameters = "-Xmx8g")
options(stringsAsFactors = FALSE)
library(reshape2)
library(Biobase)
library(plyr)
library(dplyr)
library(ggplot2)
library(limma)

eset.pdata = bind_rows(list(
  paste0(eset.folder, "bauerlab_pdata.txt") %>% 
    read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F, fill = T) %>% bldf %>%
    select(-starts_with("sample_file_")),
  paste0(eset.folder, "tcga_pdata.txt.gz") %>% 
    read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F, fill = T) %>% bldf %>%
    mutate(pdatabase_id = "tcga", patient_id = PATIENT),
  paste0(eset.folder, "tcgasmk_pdata.txt.gz") %>% 
    read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F, fill = T) %>% bldf %>%
    mutate(pdatabase_id = "tcga", patient_id = PATIENT, dataset_id = "tcgasmk"),
  paste0(eset.folder, "firesmk_pdata.txt.gz") %>% 
    read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F, fill = T) %>% bldf %>%
    mutate(pdatabase_id = "firehose"))) %>% bldf

sample.info = bind_rows(list(
  bauerlab_efit_setup(eset.contrasts, eset.design, eset.pdata, eset.folder, "eset", "gene")))
eset.info = sample.info %>% 
  select(organism_id, dataset_id, starts_with("eset_")) %>% distinct

efit.info = sample.info %>% select(organism_id, dataset_id, contrast_id, contrast_string, contrast_abbrev,
                                   starts_with("eset_"), starts_with("efit_")) %>% distinct

eset.foretinib = eset.info %>% 
  bauerlab_eset_file_check %>% filter(eset_ok) %>%
  filter(eset_id %in% c("chronic")) %>% 
  filter(eset_type == "eset") %>% distinct %>% 
  bauerlab_tbl_slice("eset_root") %>% lapply(bauerlab_eset_load, T)
efit.foretinib = efit.info %>% 
  bauerlab_efit_file_check %>% filter(efit_ok) %>%
  filter(eset_id %in% c("chronic")) %>% 
  filter(eset_type == "eset") %>% 
  bauerlab_efit_load()

foretinib.contrasts = c("acTL","chTL","exTL","lmDrug","lmTime","syResist")
foretinib.list = efit.foretinib %>% 
  mutate(efit_id = gsub("chronic_","",efit_id)) %>% 
  mutate(patient_id = gsub("(PDX|T366|T608).*","\\1",efit_id)) %>% 
  mutate(efit_id = gsub("PDX|T366|T608","",efit_id)) %>% 
  filter(efit_id %in% foretinib.contrasts) %>% 
  mutate(ftest = ftest_fdr, fpval = ftest_pval) %>% 
  arrange(fdr, pval) %>%
  group_by(efit_id,eset_id) %>% mutate(index = 1:n()) %>% ungroup %>% 
  arrange(sign(logfc), fdr, pval) %>%
  group_by(efit_id,eset_id) %>% mutate(indexdn = 1:n()) %>% ungroup %>% 
  arrange(desc(sign(logfc)), fdr, pval) %>%
  group_by(efit_id,eset_id) %>% mutate(indexup = 1:n()) %>% ungroup %>% 
  bauerlab_tbl_slice("patient_id") %>% 
  lapply(bauerlab_efit_melt,x.columns = c("logfc", "fdr", "pval", "ftest", "fpval", "index", "indexdn", "indexup"))

foretinib.table366 = efit.foretinib %>% 
  mutate(efit_id = gsub("chronic_","",efit_id)) %>% 
  filter(grepl("T366", efit_id)) %>% 
  filter(efit_id %in% c("T366acTL","T366chTL")) %>% 
  mutate(efit_id = c("T366acTL" = "acute","T366chTL" = "chronic")[efit_id]) %>% 
  mutate(ftest = ftest_fdr, fpval = ftest_pval) %>% 
  arrange(fdr, pval) %>%
  group_by(efit_id,eset_id) %>% mutate(index = 1:n()) %>% ungroup %>% 
  arrange(sign(logfc), fdr, pval) %>%
  group_by(efit_id,eset_id) %>% mutate(indexdn = 1:n()) %>% ungroup %>% 
  arrange(desc(sign(logfc)), fdr, pval) %>%
  group_by(efit_id,eset_id) %>% mutate(indexup = 1:n()) %>% ungroup %>%
  arrange(efit_id) %>% 
  bauerlab_efit_melt(x.columns = c("ave", "logfc", "fdr", "pval", "ftest", "fpval", "index", "indexdn","indexup"))

foretinib.table366 %>% rename(ave_exprs = acute_ave) %>% select(-chronic_ave, -ave_ave) %>% 
  write.table(gzfile("foretinib_expression_changes_t366.txt.gz"), sep = "\t", quote = F, row.names = F)

foretinib.conditions = c("T366", "T366acTL", "T366chTL")


foretinib.genes = data_frame(gene_symbol = c(
  "AXL","GAS6", "MET","HGF", "MST1R","MST1", "KDR","VEGFA","VEGFC","FIGF","PDGFRA","PDGFA","PDGFB")) %>% 
  mutate(gene_index = 1:n())


foretinib.boxdata = eset.foretinib %>% 
  bauerlab_eset_expression(x.features = foretinib.genes, x.melt = T) %>% 
  bind_rows %>% filter(eset_type == "eset") %>% 
  filter(eset_id %in% c("chronic")) %>% 
  filter(patient_id %in% c("T366")) %>% 
  mutate(condition_id = sample_label) %>%
  filter(sample_abbrev %in% c("T366", "T366acTL", "T366chTL")) %>%
  mutate(condition_label = c("T366" = "Control", "T366acTL" = "T+L (Acute)", "T366chTL" = "T+L (Chronic)")[sample_abbrev]) %>%
  select(eset_id, sample_filename, patient_id,
         condition_id, condition_label,
         feature_array, feature_id, probeset_id,
         gene_id, gene_symbol, gene_name, expression_value = value) %>%
  arrange(gene_symbol, eset_id, condition_id) %>% 
  group_by(eset_id) %>% 
  mutate(condition_control = condition_id == (condition_id[1])) %>% ungroup %>%
  group_by(eset_id, gene_symbol) %>% 
  mutate(control_average = mean(expression_value[condition_control])) %>% ungroup %>%
  mutate(relative_value = expression_value - control_average) %>%
  group_by(eset_id, gene_symbol, condition_id) %>% 
  mutate(condition_average = mean(expression_value),
         relative_average = mean(relative_value)) %>% ungroup %>%
  left_join(foretinib.genes) %>% 
  arrange(gene_index, gene_symbol, eset_id, condition_id)

foretinib.boxdata %>% 
  write.table(file = "2017_03_07_bauerlab_foretinib_expression_data.txt", 
              sep = "\t", quote = F, row.names = F)

allgenes.foretinib = eset.foretinib$chronic_gene_eset_ %>% 
  fData %>% as.tbl %>% select(gene_symbol) %>% distinct %>% 
  mutate(gene_index = 1:n())

foretinib.allboxdata = eset.foretinib %>% 
  bauerlab_eset_expression(x.features = allgenes.foretinib, x.melt = T) %>% 
  bind_rows %>% filter(eset_type == "eset") %>% 
  filter(eset_id %in% c("chronic")) %>% 
  filter(patient_id %in% c("T366")) %>% 
  mutate(condition_id = sample_label) %>%
  filter(sample_abbrev %in% c("T366", "T366acTL", "T366chTL")) %>%
  mutate(condition_label = c("T366" = "Control", "T366acTL" = "T+L (Acute)", 
                             "T366chTL" = "T+L (Chronic)")[sample_abbrev]) %>%
  select(eset_id, sample_filename, patient_id,
         condition_id, condition_label,
         feature_array, feature_id, probeset_id,
         gene_id, gene_symbol, gene_name, expression_value = value) %>%
  arrange(gene_symbol, eset_id, condition_id) %>% 
  group_by(eset_id) %>% 
  mutate(condition_control = condition_id == (condition_id[1])) %>% ungroup %>%
  group_by(eset_id, gene_symbol) %>% 
  mutate(control_average = mean(expression_value[condition_control])) %>% ungroup %>%
  mutate(relative_value = expression_value - control_average) %>%
  group_by(eset_id, gene_symbol, condition_id) %>% 
  mutate(condition_average = mean(expression_value),
         relative_average = mean(relative_value)) %>% ungroup %>%
  arrange(eset_id, gene_symbol, condition_id, sample_filename)

foretinib.allboxdata %>% 
  write.table(file = "2017_03_07_bauerlab_foretinib_expression_alldata.txt", 
              sep = "\t", quote = F, row.names = F)


fdr.cutoff = 0.1
expression.changes = foretinib.list %>% lapply(function(x) setNames(x, gsub("_logfc","",names(x)))) %>%
  bind_rows(.id = "patient") %>% mutate(
    significant_ftest = acTL_ftest < fdr.cutoff, 
    significant_trt_acute = (acTL_fdr < fdr.cutoff) & significant_ftest,
    significant_ctl_chronic = (exTL_fdr < fdr.cutoff) & significant_ftest,
    significant_trt_chronic = (chTL_fdr < fdr.cutoff) & significant_ftest,
    direction_trt_acute = sign(acTL) * ifelse(significant_trt_chronic, 1,0),
    direction_ctl_chronic = sign(exTL) * ifelse(significant_ctl_chronic, 1,0),
    direction_trt_chronic = sign(chTL) * ifelse(significant_trt_chronic, 1,0),
    significant_upregulated = paste0(ifelse(direction_trt_acute > 0,"A",""),
                                     ifelse(direction_ctl_chronic > 0,"R",""),
                                     ifelse(direction_trt_chronic > 0,"C","")),
    significant_dnregulated = paste0(ifelse(direction_trt_acute < 0,"A",""),
                                     ifelse(direction_ctl_chronic < 0,"R",""),
                                     ifelse(direction_trt_chronic < 0,"C","")),
    significant_differences = trimws(gsub("[\\/]$","",gsub("^[\\/]","",trimws(
      paste0(significant_upregulated, " / ", significant_dnregulated))))),
    overlap_acute_resistance = grepl("AR", significant_differences),
    overlap_acute_chronic = grepl("A[R]?C", significant_differences),
    overlap_chronic_resistance = grepl("RC", significant_differences),
    overlap_acute_chronic_resistance = grepl("ARC", significant_differences),
    significant_overlap = overlap_acute_resistance | overlap_acute_chronic | overlap_acute_chronic_resistance,
    significant_acute_only = (!significant_overlap) & (significant_trt_acute),
    significant_chronic_only = (!significant_overlap) & (significant_trt_chronic),
    significant_category = paste0(ifelse(significant_overlap, "overlap",""),
                                  ifelse(significant_acute_only, "acute",""),
                                  ifelse(significant_chronic_only, "chronic","")),
    significant_classification = ifelse(significant_overlap, "independent", ifelse(
      significant_acute_only, "acute", ifelse(significant_chronic_only, "chronic", "unaffected"))),
    expectation_max = pmax(acTL + exTL, acTL, exTL), 
    expectation_min = pmin(acTL + exTL, acTL, exTL),
    expectation = ifelse(abs(expectation_max) > abs(expectation_min), expectation_max, expectation_min),
    expectation_exceeded = abs(chTL) > abs(expectation),
    observation = chTL,
    observation_within = (observation <= expectation_max) & (observation >= expectation_min),
    observation_beyond = (observation > pmax(expectation_max, 0)) & (observation < pmin(expectation_min, 0)),
    expectation_propensity = ifelse(expectation > 0, 1, -1),
    observation_propensity = ifelse(observation > 0, 1, -1),
    treatment_propensity = ifelse(abs(expectation) > abs(observation), expectation_propensity, observation_propensity),
    significant_drug_effect = (lmDrug_fdr < fdr.cutoff) & significant_ftest,
    significant_time_effect = (lmTime_fdr < fdr.cutoff) & significant_ftest,
    significant_interaction_effect = (syResist_fdr < fdr.cutoff) & significant_ftest,
    treatment_drug = sign(lmDrug) * ifelse(significant_drug_effect, 1, 0),
    treatment_time = sign(lmTime) * ifelse(significant_time_effect, 1, 0),
    treatment_interaction = sign(syResist) * sign(treatment_propensity) * ifelse(significant_interaction_effect, 1, 0),
    treatment_synergistic = treatment_interaction > 0,
    treatment_antagonistic = treatment_interaction < 0,
    treatment_maineffects = abs(treatment_drug) + abs(treatment_time),
    treatment_factors = abs(treatment_drug) + abs(treatment_time) + abs(treatment_interaction),
    treatment_additive = (treatment_maineffects > 1) & (sign(treatment_drug) == sign(treatment_time)),
    treatment_subtractive = (treatment_maineffects > 1) & (sign(treatment_drug) != sign(treatment_time)),
    treatment_dampened = (treatment_maineffects > 0) & treatment_antagonistic & (abs(expectation) > abs(observation)),
    treatment_augmented = (treatment_maineffects > 0) & treatment_synergistic & (abs(expectation) < abs(observation)),
    treatment_independent = (treatment_interaction == 0) & (treatment_subtractive == 0),#!treatment_dampened & 
    significant_factors = paste0(ifelse(treatment_synergistic,"S",""),
                                 ifelse(treatment_antagonistic,"A",""),
                                 ifelse(treatment_drug,"D",""),
                                 ifelse(treatment_time,"T","")),
    category = ifelse(!treatment_independent, 
                      paste0(ifelse(treatment_synergistic, paste0("S"),""),
                             ifelse(treatment_antagonistic, paste0("A"),"")),
                      paste0(ifelse(treatment_drug != 0, paste0("D"),""),
                             ifelse(treatment_time != 0, paste0("T"),""))),
    classification = ifelse(!treatment_independent, 
                            paste0(ifelse(treatment_synergistic, paste0("synergistic"),""),
                                   ifelse(treatment_antagonistic, paste0("antagonistic"),"")),
                            ifelse(treatment_maineffects > 0, "independent", "unaffected")),
    classification = ifelse(!treatment_independent, 
                            paste0(ifelse(treatment_synergistic, paste0("synergistic"),""),
                                   ifelse(treatment_antagonistic, paste0("antagonistic"),"")),
                            ifelse(treatment_maineffects > 0, "independent", "unaffected")),
    treatment_classification = ifelse(treatment_dampened,"dampened", ifelse(treatment_augmented, "augmented",classification)),
    pclass = significant_classification,
    gclass = classification,
    ggclass = treatment_classification)


expression.changes %>% count(paste0(
  ifelse(significant_drug_effect, ifelse(treatment_drug > 0, "D+", "T-"), ""),
  ifelse(significant_time_effect, ifelse(treatment_time > 0, "T+", "T-"), ""),
  ifelse(significant_interaction_effect, ifelse(treatment_interaction > 0, "i+", "i-"), "")))


