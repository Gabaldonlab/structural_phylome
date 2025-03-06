library(tidyverse)
library(glue)

n_sample <- 500
min_seqs <- 4

params <- yaml::read_yaml("config/params.yaml")
params_dataset <- yaml::read_yaml("config/Hsap_draft.yaml")
tax <- read_delim(glue("results/data/meta/",
                       params_dataset$homology_dataset, "_uniprot_genomes.tsv"), 
                  show_col_types=F)

# some sequences do not have structures so to simplify our lives 
# we only get those ortthogroups with at least 4 existying files
db_fasta <- Biostrings::readAAStringSet("results/homology/hsap_euka/db/all_seqs_fsdb_ss.fa")

dbs <- read_delim(paste0(params$data_dir, "ids/", tax$Proteome_ID, "_ids.tsv"),
           delim = "\t", id = "proteome", show_col_types = FALSE) %>% 
  mutate(proteome = gsub("_ids.tsv", "", basename(proteome)), 
         target=paste0("AF-", Entry, "-F1")) %>% 
  select(target, OMA, eggNOG) %>% 
  filter(target %in% names(db_fasta)) %>% 
  # as one id can have multiple eggnog we need to unnest
  mutate(eggNOG=strsplit(eggNOG, ";"),
         OMA=gsub(";", "", OMA)) %>%
  unnest(eggNOG)

eggNOGs <- dbs %>% 
  select(target, eggNOG) %>% 
  filter(!is.na(eggNOG)) %>% 
  group_by(eggNOG) %>% 
  summarise(dist=n_distinct(target)) %>% 
  filter(dist>min_seqs) %>% 
  pull(eggNOG) %>% 
  sample(n_sample, replace = FALSE)

dbs %>% 
  filter(eggNOG %in% eggNOGs) %>% 
  select(eggNOG, target) %>% 
  arrange(eggNOG) %>% 
  write_delim("data/seeds/eggnog_benchmark.txt", col_names = FALSE, delim = "\t")

OMAs <- dbs %>% 
  select(target, OMA) %>% 
  filter(!is.na(OMA)) %>% 
  group_by(OMA) %>% 
  summarise(dist=n_distinct(target)) %>% 
  filter(dist>min_seqs) %>% 
  pull(OMA) %>% 
  sample(n_sample, replace = FALSE)

dbs %>% 
  filter(OMA %in% OMAs) %>% 
  distinct(target, OMA) %>% 
  select(OMA, target) %>% 
  arrange(OMA) %>% 
  write_delim("data/seeds/oma_benchmark.txt", col_names = FALSE, delim = "\t")

