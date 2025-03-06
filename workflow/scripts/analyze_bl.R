suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidyverse))
source("workflow/scripts/functions.R")

theme_set(theme_classic())

# all trees
df_trees <- read_delim(snakemake@input[["trees"]], delim = "\t", show_col_types = FALSE,
                       col_names = c("gene", "targets", "alphabet", "model", "tree")) %>% 
  mutate(model=factor(model, levels=models)) %>% 
  filter(targets=="common")

# Read all trees
ts <- read.tree(text = df_trees$tree)
names(ts) <- paste(df_trees$gene, df_trees$targets, df_trees$alphabet, df_trees$model, sep = "_")

fs <- read_delim(snakemake@input[["fs"]], 
                 col_names = fs_columns, show_col_types = FALSE) %>% 
  filter(query %in% paste0("AF-", df_trees$gene, "-F1")) 
  # filter(evalue<snakemake@params[["eval_both"]], 
  #        qcov*100>snakemake@params[["coverage"]], 
  #        tcov*100>snakemake@params[["coverage"]])

blast <- read_delim(snakemake@input[["blast"]], 
                    col_names = blast_columns, show_col_types = FALSE) %>% 
  filter(query %in% paste0("AF-", df_trees$gene, "-F1")) #%>% 
  # filter(evalue<snakemake@params[["eval_both"]], 
  #        length/qlen*100>snakemake@params[["coverage"]], 
  #        length/slen*100>snakemake@params[["coverage"]])


df <- full_join(blast, fs,
                by = c("query", "target"), suffix=c("_blast", "_fs")) %>%
  mutate(pident_fs = 100*pident_fs, 
         qcov_fs = 100*qcov_fs,
         singleton = case_when(is.na(evalue_fs) ~ "only_blast",
                               is.na(evalue_blast) ~ "only_fs", 
                               .default = "common"))

df_red <- select(df, query, target, singleton) 

res <- NULL
for (idx in 1:nrow(df_trees)){
  t <- ts[[idx]]
  seed <-  t$tip.label[grepl(df_trees$gene[idx], t$tip.label)]
  # cophe <- cophenetic(t)
  # total <- sum(cophe)
  # bl <- cophe[seed,]
  nn <- as.matrix(adephylo::distTips(t, method = "nNodes"))[seed, ]
  max_dist <- max(nn)

  row <- enframe(nn, name = "target", value="nn_dist") %>% 
    # left_join(enframe(bl, name = "target", value = "bl_dist"), by="target") %>% 
    mutate(query=seed, 
           max_dist=max_dist,
           # n_splits=TreeTools::NSplits(t),
           # total_length=total,
           targets=df_trees$targets[idx],
           alphabet=df_trees$alphabet[idx],
           model=df_trees$model[idx]) %>% 
    left_join(df_red, by=c("query", "target"))
  res <- bind_rows(res, row)
}

# for what this should be grouped?
plot_dist <- res %>% 
  # filter(query=="AF-Q01718-F1") %>% 
  group_by(query, singleton, model) %>%
  summarise(median_dist=median(nn_dist/max_dist), n=n()) %>% 
  ggplot(aes(median_dist, model, fill=singleton)) + 
  geom_violin(alpha=0.8, draw_quantiles = 0.5) +
  # facet_grid(model~.) + 
  scale_x_continuous(expand = c(0,0.01), limits = c(0,1)) +
  theme(panel.border = element_rect(fill = NA)) + 
  scale_fill_manual(values = palette_singleton) + 
  theme(legend.position = "bottom")

ggsave(snakemake@output[[1]], plot_dist, width = 5, height = 8)