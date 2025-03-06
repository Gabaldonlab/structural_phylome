library(ape)

# For our sampling some species are not in eggnog or OMA, we can prune them 
# from the sptree.
params <- yaml::read_yaml("config/params.yaml")
params_dataset <- yaml::read_yaml("config/Hsap_draft.yaml")
tax <- read_delim(glue("results/data/meta/",
                       params_dataset$homology_dataset, "_uniprot_genomes.tsv"), 
                  show_col_types=F)
taxmap <- read_delim(glue(homodir, "/db/taxidmap"), 
                     delim = "\t", col_names = c("label", "Tax_ID"),
                     show_col_types = FALSE) %>% 
  left_join(tax, by = "Tax_ID") %>% 
  dplyr::select(label, Proteome_ID, mnemo)

euka_tree <- read.tree(params_dataset$species_tree)

dbs <- read_delim(paste0(params$data_dir, "ids/", tax$Proteome_ID, "_ids.tsv"),
                  delim = "\t", id = "proteome", show_col_types = FALSE) %>% 
  mutate(proteome = gsub("_ids.tsv", "", basename(proteome))) %>% 
  mutate(target=paste0("AF-", Entry, "-F1")) %>% 
  select(-proteome, -Entry)

missing_df <- dbs %>%
  right_join(taxmap, by=c("target"="label")) %>%
  group_by(mnemo) %>%
  select(mnemo, eggNOG, OMA) %>% 
  summarise_all(~ sum(is.na(.))/n())

oma_missing <- missing_df[missing_df$OMA==1,]$mnemo
eggnog_missing <- missing_df[missing_df$eggNOG==1,]$mnemo

oma_tree <- TreeTools::DropTip(euka_tree, oma_missing)
write.tree(oma_tree, "data/sptrees/homo_oma.nwk")
eggnog_tree <- TreeTools::DropTip(euka_tree, eggnog_missing)
write.tree(eggnog_tree, "data/sptrees/homo_eggnog.nwk")

# get QFO species tree and do some small modifications for the 
# Foldtree benchmark
t <- treeio::read.nhx("resources/species_tree_QFO2020.nx")
t <- TreeTools::Preorder(t@phylo)
t$tip.label <- gsub("_.*", "", t$tip.label)
t <- TreeTools::Subtree(t, getMRCA(t, tip = c("HUMAN", "GIAIC")))
t <- TreeTools::AddTip(t, where = "CANFA", label = "CANLF")
t <- TreeTools::AddTip(t, getMRCA(t, tip = c("DANRE", "ORYLA")), label = "LEPOC")
t <- TreeTools::AddTip(t, getMRCA(t, tip = c("DROME", "BOMMO")), label = "TRICA")
t <- TreeTools::AddTip(t, where = "TRIVA", label = "TRIV3")
t <- TreeTools::AddTip(t, where = "SCHMA", label = "HELRO")
t <- TreeTools::AddTip(t, where = "SORBI", label = "MAIZE")
write.tree(t, "data/sptrees/QFO.nwk")

