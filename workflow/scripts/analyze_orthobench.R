suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
source("workflow/scripts/functions.R")
theme_set(theme_classic())

ids_df <- read_delim(snakemake@input[["ids"]], 
                     show_col_types = FALSE, delim = "\t")
# ids_df <- read_delim(paste0("results/data/ids/", tax$Proteome_ID, "_ids.tsv"))

eval <- snakemake@params[["evalue"]]
# eval <- 1e-3
min_size <- snakemake@params[["min_seqs"]]
# min_size <- 5
dataset <- snakemake@wildcards[["db"]]
# dataset <- "eggnog"

seed_fa <- Biostrings::readAAStringSet(snakemake@input[["fa"]])
# seed_fa <- Biostrings::readAAStringSet("results/ortho_bench/eggnog/db/db.fa")

length_df <- Biostrings::width(seed_fa)
names(length_df) <- names(seed_fa)

colnames(ids_df) <- str_to_lower(colnames(ids_df))

true_set <- ids_df %>% 
  mutate(entry = paste0("AF-", entry, "-F1")) %>% 
  filter(entry %in% names(length_df)) %>%
  select(entry, match(dataset, colnames(ids_df))) %>% 
  rename(dataset = 2) %>% 
  rowwise() %>%
  mutate(dataset = strsplit(dataset, ";")) %>%
  unnest(dataset)

# fs_bench <- read_delim("results/ortho_bench/eggnog//hits_fs.tsv", col_names = fs_columns) %>%
#   filter(evalue < eval, query != target) %>%
#   mutate(method = "fs", qcov = qcov * 100)
fs_bench <- read_delim(snakemake@input[["fs"]], col_names = fs_columns, 
                       show_col_types = FALSE, delim = "\t") %>% 
  filter(evalue < eval, query != target) %>% 
  mutate(method = "fs", qcov = qcov * 100)

# blast_bench <- read_delim("results/ortho_bench/eggnog/hits_blast.tsv", col_names = blast_columns)  %>%
#   filter(evalue < eval, query != target) %>%
#   mutate(method = "blast")
blast_bench <- read_delim(snakemake@input[["blast"]], col_names = blast_columns, 
                          show_col_types = FALSE, delim = "\t")  %>% 
  filter(evalue < eval, query != target) %>% 
  mutate(method = "blast")

comb_df <- true_set %>% 
  # head(1000) %>% 
  group_by(dataset) %>%
  filter(n() > min_size) %>% 
  summarize(pairs = list(combn(entry, 2, simplify = FALSE)), .groups = "drop") %>%
  unnest(pairs) %>% 
  mutate(entry1 = map_chr(pairs, 1), entry2 = map_chr(pairs, 2)) %>% 
  select(dataset, entry1, entry2)  %>%
  bind_rows(rename(., entry1 = entry2, entry2 = entry1)) %>% 
  left_join(select(blast_bench, query, target, qcov, method), 
          by = c("entry1"="query", "entry2"="target")) %>% 
  left_join(select(fs_bench, query, target, qcov, method), 
            by = c("entry1"="query", "entry2"="target"), 
            suffix = c("_blast", "_fs"))


accuracy_df <- comb_df %>% 
  group_by(dataset) %>% 
  left_join(enframe(length_df, name="entry1", value="length")) %>% 
  # left_join(enframe(length_df, name="entry2", value="length_t")) %>% 
  summarise(recall_blast = sum(!is.na(qcov_blast))/n(),
            recall_fs = sum(!is.na(qcov_fs))/n(), 
            n = n_distinct(entry1),
            avg_len = mean(length),
            avg_qcov_blast = mean(qcov_blast, na.rm = T),
            avg_qcov_fs = mean(qcov_fs, na.rm = T))


qcov_box <- comb_df %>% 
  pivot_longer(cols = c(qcov_blast, qcov_fs), values_to = "qcov") %>% 
  mutate(class = case_when(is.na(method_blast) & is.na(method_fs) ~ NA,
                           is.na(method_fs) ~ "only Blast",
                           is.na(method_blast) ~ "only Foldseek",
                           .default = "Common")) %>% 
  filter(!is.na(class)) %>% 
  ggplot(aes(qcov, class, fill=class)) + 
  geom_violin() +
  # facet_wrap(~name) +
  geom_text(data = . %>% count(class), 
            aes(x = 15, label = paste0("n: ", n/2)),
            nudge_y = -0.2, hjust = 0, nudge_x = 1) +
  scale_fill_manual(values = palette_singleton_p) + 
  labs(y="", x = "Query coverage")  + 
  # scale_x_continuous(expand = expansion(add = c(5, 35))) +
  theme(legend.position = "none")

recall_plot <- ggplot(accuracy_df, aes(recall_blast, recall_fs)) + 
  # geom_point(pch=21) + 
  geom_bin2d(binwidth = .05) + 
  stat_bin2d(geom = "text", aes(label = after_stat(count)), 
             binwidth = .05, color="white", size=2) +
  geom_abline(slope = 1) +
  scale_fill_viridis_c() +
  labs(x = "Recall Blast", y = "Recall Foldseek") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  theme(legend.position = "none")
# 
# recall_length <- accuracy_df %>% 
#   pivot_longer(contains("recall")) %>% 
#   mutate(name=gsub(".*._", "", name)) %>% 
#   ggplot(aes(avg_len, value)) + 
#   geom_bin2d(binwidth = c(25, 0.05)) + 
  # geom_smooth() +
  # stat_bin2d(geom = "text", aes(label = after_stat(count)), 
  #            binwidth = .05, color="white", size=2) +
  # facet_wrap(~name, nrow = 2) + 
  # scale_fill_viridis_c() +
  # labs(x="Average length", y="Recall")

# outplot <- (recall_plot | qcov_box) / (recall_length) + 
#   plot_annotation(title = paste(str_to_title(dataset), " benchmark"), 
#                   tag_level = "A") + 
#   plot_layout(heights = c(1, 0.7))

outplot <- (recall_plot | qcov_box ) + 
  plot_annotation(#title = paste(str_to_title(dataset), " benchmark"), 
                  tag_level = "A") + 
  plot_layout(widths = c(1, 0.5)) &
  theme(plot.tag = element_text(size = 10, face = "bold", family = "Helvetica"))

# ggsave("draft/figures/Fig_recall.pdf", outplot, width = 6.5, height = 4)
ggsave(snakemake@output[[1]], outplot, width = 6.5, height = 4)


  