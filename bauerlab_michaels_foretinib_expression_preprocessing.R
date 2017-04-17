
options(java.parameters = "-Xmx4g")
options(stringsAsFactors = FALSE)
library(oligo)
library(genefilter)
library(hgu133plus2.db)
library(hugene10sttranscriptcluster.db)
library(hugene20sttranscriptcluster.db)
library(pd.hg.u133.plus.2)
library(pd.hugene.1.0.st.v1)
library(pd.hugene.2.0.st)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(xlsx)

source("bauerlab_helper.R")

eset.folder = c("C:/Users/edikb/Google Drive/results/bauerlab/")
cels.folder = c("C:/Users/edikb/Google Drive/affydata/")

sample.info = paste0(eset.folder, "bauerlab_pdata.txt") %>% 
  read.table(sep = "\t", quote = "", header = T, stringsAsFactors = F, check.names = F, fill = T) %>% bldf %>%
  mutate(sample_file_full = paste0(cels.folder, sample_id),
         sample_file_exists = sample_file_full %>% file.exists) %>%
  mutate(Smoker = ifelse(!is.na(Smoking) & (Smoking == "formersmoker"), "smoker", Smoker))

sample.info %>% count(sample_unique,dataset_id) %>% ungroup %>% arrange(desc(n))
sample.info %>% count(dataset_id,sample_file_exists) %>% data.frame
nrow(sample.info %>% filter(!sample_file_exists)) # 0

eset.selection = sample.info$dataset_id %>% unique
eset.version = "gene"

eset.check = data_frame(dataset_id = c(sample.info$dataset_id,eset.selection)) %>% distinct %>% 
  mutate(core_filename = paste0(eset.folder, dataset_id,"_",eset.version,"_core_exprs.txt.gz"),
         eset_filename = paste0(eset.folder, dataset_id,"_",eset.version,"_eset_exprs.txt.gz"),
         nset_filename = paste0(eset.folder, dataset_id,"_",eset.version,"_nset_exprs.txt.gz"),
         core_status = file.exists(core_filename),
         eset_status = file.exists(eset_filename),
         nset_status = file.exists(nset_filename))
eset.check %>% count(core_status, eset_status, nset_status)
eset.check %>% filter(!eset_status)

eset.status = eset.selection %>% lapply(bauerlab_eset_preprocess, sample.info, x.version = "gene", x.folder = eset.folder)

eset.status %>% bind_rows %>% count(dataset_status)

qc.status = eset.selection %>% lapply(bauerlab_eset_quality, x.version = "gene", x.type = "core", x.folder = eset.folder)

