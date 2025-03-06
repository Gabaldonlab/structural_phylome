suppressPackageStartupMessages(library(tidyverse))

min_size <- snakemake@params[["min_seqs"]]
dataset <- snakemake@wildcards[["db"]]

tax <- read_delim(snakemake@input[["meta"]], 
                  show_col_types=F)

in_db <- read_delim(snakemake@input[["idmap"]], 
                    show_col_types = FALSE,col_names = c("entry", "taxid"))

# remove all those that are not in your dbs anyway!
ids_df <- read_delim(snakemake@input[["ids"]], show_col_types = FALSE) 
colnames(ids_df) <- str_to_lower(colnames(ids_df))

ids_df <- ids_df %>% 
  mutate(entry = paste0("AF-", entry, "-F1")) %>% 
  left_join(in_db) %>% 
  filter(!is.na(taxid)) %>% 
  select(entry, match(dataset, colnames(ids_df))) %>% 
  rename(dataset = 2) %>% 
  filter(!is.na(dataset)) %>% 
  rowwise() %>%
  mutate(dataset = strsplit(dataset, ";")) %>%
  unnest(dataset)


# get groups with at least X members
groups_df <- ids_df %>% 
  count(dataset) %>% 
  filter(n > min_size)

ids_df %>%
  right_join(groups_df) %>%  
  pull(entry) %>% 
  unique() %>% 
  writeLines(snakemake@output[[1]])

