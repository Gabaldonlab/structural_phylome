suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(patchwork))
source("workflow/scripts/functions.R")
theme_set(theme_classic())

# read species tree
sptree <- read.tree(snakemake@input[["sptree"]])

taxidmap <- read_delim(snakemake@input[["taxidmap"]], 
                       col_names = c("target", "Tax_ID"), show_col_types = FALSE)

# trees Dataframe
df_trees <- read_delim(snakemake@input[["trees"]],
                       show_col_types = FALSE,
                       col_names = c("gene", "target", "alphabet", "model", "tree")) %>% 
  mutate(model=factor(model, levels=models))

# ML trees with selected model
df_ml <- read_delim(snakemake@input[["mltrees"]], delim = "\t", 
                    show_col_types=FALSE) %>% 
  # mutate(fn=gsub("_mltrees.txt", "", basename(fn))) %>% 
  # separate(fn, c("gene", "target", "alphabet", "model"), sep = "_") %>% 
  mutate(model=factor(gsub("3DI", "3Di", gsub("GTR20", "GTR", gsub("\\+.*", "", Model))), levels=models), 
         rate=sub("\\+", "", gsub("GTR20|3DI|LG|AF|LLM", "", Model)))

# Disco similarity to sptree
disco_fls <- c(snakemake@input[["disco"]])
# disco_fls <- list.files("results/phylogeny/hsap_1kseeds/reco/disco/",
#                         pattern = "*output*", full.names = T)

# Astral-pro results
clades <- read_delim(snakemake@input[["table"]], show_col_types = FALSE)

apro_files <- snakemake@input[["apro_trees"]]
# apro_files <- list.files("results/phylogeny/hsap_1kseeds/reco/apro//",
#                         pattern = "*support*", full.names = T)
apro_trees <- read.tree(text = sapply(apro_files, readLines))
names(apro_trees) <- sapply(str_split(basename(apro_files), "_"), 
                            function(x) paste0(x[2], "_", x[4]))

# 1. Compute astral-pro quartet
nodes_clades <- fortify(sptree) %>% 
  # filter(!isTip) %>% 
  mutate(ordered = rank(y)) %>% 
  mutate(label_int = fct_reorder(label, y)) %>%
  select(label_int, ordered)

trees_df <- fortify(apro_trees) 
nodes_df <- trees_df %>% 
  group_by(.id) %>% 
  mutate(ordered = rank(y)) %>% 
  left_join(nodes_clades) %>% 
  separate(.id, c("target", "model")) %>% 
  filter(!isTip, label!="") %>% 
  mutate(label = gsub("\\[|\\]|\\'", "", str_replace_all(label, "[a-z]{1,2}[0-9]=", ""))) %>% 
  separate(label, c("pp1", "pp2", "pp3", "f1", "f2", "f3", "q1", "q2", "q3"), 
           ";", convert = TRUE) %>% 
  left_join(nodes_clades) %>% 
  mutate(freq=f1+f2+f3, 
         model=factor(model, levels=models))

plot_apro <- nodes_df %>%
  # filter(model!="FTPY") %>% 
  ggplot(aes(model, q1, fill=target)) + 
  geom_boxplot() + 
  scale_fill_manual(values = palettes_method) +
  labs(x="", y="Quartet support", fill="Target sets") +
  theme(legend.position = "none")


disco_rf <- get_disco_rf(disco_fls, sptree)

ts <- read.tree(text = df_trees$tree)
names(ts) <- paste(df_trees$gene, df_trees$target, df_trees$alphabet, df_trees$model, sep = "_")
n_tips <- sapply(ts, function(x) length(x$tip.label)) %>% 
  enframe() %>%
  separate(name, c("gene", "target", "alphabet", "model"), sep = "_") %>%
  rename("n_tips"="value")

# Reconciliation notung
reco <- read_delim(snakemake@input[["reco"]], 
                   show_col_types = FALSE, 
                   col_names = c("gene", "target", "alphabet", "model", "dups", "losses")) %>% 
  # separate(gene, c("gene", "target", "alphabet", "model"), sep = "_") %>% 
  # left_join(var_df) %>% 
  left_join(n_tips) %>%
  # left_join(n_sps) %>% 
  mutate(model=factor(model, levels=models),
         n_events = (dups+losses)/n_tips) 

# scores <- read_delim(snakemake@input[["scores"]], show_col_types = FALSE) %>% 
#   left_join(reco) %>% 

# Ranger results
plot_DL <- reco %>% 
  # filter(model!="FTPY") %>% 
  ggplot(aes(model, n_events, fill=target)) +
  geom_boxplot(outlier.size = .1) +
  labs(y="(D+L)/# Tips", x="",fill="Target sets") +
  # ggthemes::scale_fill_solarized() +
  scale_fill_manual(values = palettes_method) +
  theme(legend.position = "none")

plot_rf <- disco_rf %>% 
  filter(ntips>10) %>%
  # filter(ntips>=length(sptree$tip.label)/2) %>%
  separate(bn, c("seed", "target", "alphabet", "model")) %>% 
  mutate(model=factor(model, levels=models)) %>% 
  ggplot(aes(model, rf, fill=target)) + 
  geom_boxplot() + 
  scale_fill_manual(values = palettes_method) +
  labs(x="", y="Norm. RF to species tree (only if > 10 tips)", fill="Target sets")  +
  theme(legend.position = "none")

rank_plot <- reco %>% 
  group_by(gene, target) %>% 
  mutate(rank = factor(dense_rank(n_events), levels = seq(0,9,0.5))) %>% 
  ggplot(aes(rank, fill=model, alpha=target)) +
  # geom_bar(position = "dodge") + 
  geom_histogram(stat="count", position = "dodge") +
  # facet_grid(target~.) +
  scale_y_continuous(expand = expansion(c(0,0))) +
  scale_fill_manual(values = palettes_model) + 
  scale_alpha_manual(values = seq(1, 0.3, length.out=4))

reco_plot <- ((plot_rf | plot_apro | plot_DL) / (rank_plot)) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom") 

ggsave(snakemake@output[["reco"]], reco_plot, width = 8, height = 6)


# nodes_df %>% 
#   group_by(node, target) %>% 
#   mutate(rank = factor(dense_rank(-q1), levels = seq(0,9,0.5))) %>% 
#   ggplot(aes(rank, fill=model, alpha=target)) +
#   geom_histogram(stat="count", position = "dodge") +
#   # facet_grid(target~.) +
#   scale_y_continuous(expand = expansion(c(0,0))) +
#   scale_fill_manual(values = palettes_model) + 
#   scale_alpha_manual(values = seq(1, 0.3, length.out=4))

