suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
source("workflow/scripts/functions.R")

theme_set(theme_classic())

taxidmap <- read_delim(snakemake@input[["taxidmap"]], 
                       col_names = c("target", "Tax_ID"), show_col_types=F)


# tax <- read_delim(snakemake@input[["groups"]], show_col_types=F)
# # read input table
# lvls_tax <- rev(c(unique(tax$Clade), "archaea", "bacteria"))
# 
# table_columns <- c('Proteome_ID', 'Tax_ID', 'count1', 'count2', 'count3', 
#                    'genome', 'source', 'species', 'mnemo')

table <- read_delim(snakemake@input[["table"]], show_col_types=F, 
                    delim = "\t") 
  # left_join(tax, by = "Tax_ID") %>% 
  # mutate(Clade = ifelse(is.na(Clade), mnemo, Clade),
  #        Clade = factor(Clade, levels = lvls_tax))


# seeds <- readLines("results/test/ids/UP000005640_common.ids")


fs <- read_delim(snakemake@input[["fs"]], 
                 col_names = fs_columns, show_col_types = FALSE) %>% 
  # filter(query %in% seeds) %>%
  mutate(evalue = ifelse(evalue==0, 1e-180, evalue)) %>%
  # filter(query!=target) %>% 
  left_join(taxidmap) %>% 
  filter(evalue<snakemake@params[["eval_both"]], 
         qcov*100>snakemake@params[["coverage"]], 
         tcov*100>snakemake@params[["coverage"]]) %>% 
  group_by(query) %>% 
  filter(row_number()<=snakemake@params[["max_seqs"]]) %>% 
  ungroup()

blast <- read_delim(snakemake@input[["blast"]], 
                    col_names = blast_columns, show_col_types = FALSE) %>% 
  # filter(query %in% seeds) %>%
  mutate(evalue = ifelse(evalue==0, 1e-180, evalue)) %>% 
  filter(evalue<snakemake@params[["eval_both"]], 
         length/qlen*100>snakemake@params[["coverage"]], 
         length/slen*100>snakemake@params[["coverage"]]) %>% 
  group_by(query) %>% 
  filter(row_number()<=snakemake@params[["max_seqs"]]) %>% 
  ungroup()
# $11<{params.eval_both} && $4/$15*100>{params.coverage} && $4/$16*100


df <- full_join(blast, fs,
          by = c("query", "target"), suffix=c("_blast", "_fs")) %>%
  mutate(pident_fs = 100*pident_fs, 
         qcov_fs = 100*qcov_fs,
         singleton = case_when(is.na(evalue_fs) ~ "only_blast",
                               is.na(evalue_blast) ~ "only_fs", 
                               .default = "common"))

# first of all let's explore common properties:
# bitscore, evalue, pident, prop mismatch and prop gaps
plot_common <- df %>% 
  select(query, singleton, 
         contains(c("pident", "evalue", "mismatch", "gap", "length"))) %>% 
  mutate(evalue_blast=-log10(evalue_blast),
         evalue_fs=-log10(evalue_fs),
         propgaps_blast=gapopen_blast/length_blast,
         propgaps_fs=gapopen_fs/length_fs,
         propmismatch_blast=mismatch_blast/length_blast,
         propmismatch_fs=mismatch_fs/length_fs) %>% 
  # head(10000) %>%
  select(-starts_with(c("mismatch", "length", "gap"))) %>% 
  pivot_longer(cols = !c(query, singleton)) %>% 
  separate(name, c("what", "method"), "_") %>% 
  ggplot(aes(value, color=singleton, linetype=method)) + 
  geom_density() +
  facet_wrap(~what, scales = "free") + 
  scale_color_manual(values = palette_singleton) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_linetype_manual(values = c(2,1)) + 
  labs(subtitle = "Common properties distribution") + 
  theme(legend.position = "bottom")


plot_fs <- df %>% 
  select(query, singleton, lddt, alntmscore, rmsd, prob) %>% 
  pivot_longer(cols = !c(query, singleton)) %>% 
  ggplot(aes(value, color=singleton)) + 
  geom_density() +
  facet_wrap(~name, scales = "free") + 
  scale_color_manual(values = palette_singleton) + 
  scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(subtitle = "Foldseek properties distribution") + 
  theme(legend.position = "none")


plot_prop <- group_by(df) %>% 
  count(query, singleton) %>% 
  group_by(query) %>% 
  mutate(nn=sum(n), prop=n/nn) %>% 
  ggplot(aes(singleton, prop, color=singleton)) + 
  geom_boxplot() +
  labs(subtitle = "Proportion of singletons") + 
  scale_color_manual(values = palette_singleton)+ 
  theme(legend.position = "none")


# plot_point <- group_by(df) %>% 
#   count(query, singleton) %>% 
#   pivot_wider(id_cols = query, names_from = singleton, values_from = n) %>% 
#   ggplot(aes(only_blast, only_fs, size=common, fill=common)) + 
#   geom_point(pch=21, alpha=.8) + 
#   # coord_fixed() +
#   # facet_grid(~proceed) +
#   scale_fill_viridis_c()+ 
#   labs(subtitle = "Scatter of singletons") + 
#   theme(legend.position = "bottom")

out <- (plot_prop | (plot_fs/plot_common)) + 
  plot_layout(widths = c(1,3))
ggsave(snakemake@output[[1]], out, width = 10, height = 6)

# it would be cool to check if positions overlap? 
# if not, are they different from the common ones? i.e. another local hit
