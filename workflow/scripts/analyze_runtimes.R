suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
source("workflow/scripts/functions.R")
theme_set(theme_bw())

# all_fls <- list.files("results/phylogeny/hsap_1kseeds/benchmarks/", full.names = T, 
#            include.dirs = F, pattern = "*txt", recursive = T) 
# 
# stats <- read_delim("tmp.stats") %>% 
#   separate(file, c("gene", "targets", "alphabet"), sep = "_") %>%
#   mutate(alphabet = gsub(".alg", "", alphabet)) %>% 
#   separate(alphabet, c("alphabet", "clean")) %>% 
#   filter(!is.na(clean))
# 
# tmp <- read_delim(all_fls[1])
# df_tmp <- sapply(all_fls, function(x) readLines(x)[-1]) %>% 
#   enframe() %>% 
#   separate(value, colnames(tmp), convert = TRUE) 
# 
# df <- df_tmp %>% 
#   mutate(method = basename(dirname(name)), name = basename(name),
#          max_rss_gb = max_rss/1024) %>% 
#   separate(name, c("seed", "gene", "targets", "alphabet", "model"), sep = "_") %>% 
#   mutate(targets = gsub(".txt", "", targets),
#          alphabet = gsub(".txt", "", alphabet),
#          model = gsub(".txt", "", 
#                       case_when(method=="foldtree" ~ NA,
#                                 method=="fastme" ~ NA, 
#                                 .default = model))) %>% 
#   left_join(stats)
# 
# 
# # Alignment
# # df %>% 
# #   filter(method %in% c("mafft", "foldmason")) %>% 
# #   select(targets, alphabet, `h:m:s`, num_seqs, avg_len, method, max_rss_gb) %>% 
# #   ggplot(aes(avg_len, `h:m:s`, color=targets)) + 
# #   geom_point() + 
# #   facet_grid(~method)
# 
# 
# 
# df <- read_delim(snakemake@input[["time"]], 
#                  show_col_types = FALSE) %>% 
#   mutate(basename=gsub(".txt", "", basename),
#          step=basename(dirname),
#          max_rss_gb = max_rss/1024) %>% 
#   separate(basename, c("seed", "gene", "method", "alphabet", "model"))
# 
# 
# df %>%
#   filter(method %in% c("fastme", "iqtree")) %>%
#   select(targets, model, `h:m:s`, num_seqs, avg_len, method, max_rss_gb) %>% 
#   mutate(model = ifelse(is.na(model), "FM", model)) %>% 
#   ggplot(aes(num_seqs*avg_len, `h:m:s`, color=model)) + 
#   geom_point() + 
#   geom_smooth()
#   scale_color_manual(values = palettes_model)

aln_stats <- read_delim(snakemake@input[["aln"]], show_col_types = FALSE) %>%
  mutate(trim = ifelse(grepl("clean", file), "clean", "aln"),
         file = gsub("\\..*", "", file),
         gap_prop = sum_gap/(num_seqs*min_len)) %>% 
  separate(file, into = c("gene", "method", "alphabet"), "_")

# is blast slower than fs
plot_homo <- filter(df, step=="homology") %>% 
    mutate(global=is.na(method), method=ifelse(is.na(method), gene, method)) %>% 
  ggplot(aes(method, y=`h:m:s`, fill=method)) + 
  geom_jitter(data = . %>% filter(!global), size=2, pch=21, width = .1) +
  geom_point(data = . %>% filter(global), size=4, pch=21) + 
  labs(subtitle = "Homology search", y = "Time", x="Method") + 
  scale_fill_manual(values = palettes_method) + 
  theme(legend.position = "none")


# is the alignment different?
plot_aln <- filter(df, step=="mafft") %>%
  left_join(filter(aln_stats, trim=="aln")) %>% 
  ggplot(aes(num_seqs*min_len, `h:m:s`, shape=method, color=alphabet)) +
  geom_point(size=1) + 
  # geom_smooth(se=F) +
  labs(y="Time", subtitle = "Mafft --auto") +
  scale_color_manual(values = palettes_alphabet) + 
    labs(x="seqs*length") +
  theme(legend.position = "bottom")

# is iqtree slower than ft?
# plot_memtree <- filter(df, step %in% c("foldtree", "iqtree", "quicktree")) %>%
#   ggplot(aes(`h:m:s`, color=step)) +
#   geom_boxplot()
#   geom_point(alpha=.4, size=1) + 
#   labs(y="Time", x="Mem (GB)", subtitle="Foldtree vs IQTree") + 
#   scale_color_manual(values = palettes_step) + 
#   theme(legend.position = "bottom")
format_hm <- function(sec) stringr::str_sub(format(sec), end = -8L)
# foldseek time and memory usage
plot_fs <- filter(df, step=="foldseek") %>%
  left_join(filter(aln_stats, trim=="clean")) %>% 
  ggplot(aes(method, `h:m:s`, fill=method)) + 
  geom_boxplot() + 
  scale_fill_manual(values = palettes_method) + 
  labs(x="Targets", y="Time", subtitle = "Foldseek all v all") + 
  scale_y_time(labels = format_hm) + 
  theme(legend.position = "none")

# if so how much does it change between alphabets?
plot_iqtree <- filter(df, step=="iqtree") %>%
  left_join(filter(aln_stats, trim=="clean")) %>% 
  ggplot(aes(log(min_len*num_seqs), `h:m:s`, color=model)) + 
  geom_point(alpha=0.4, size=1) +
  geom_smooth(se=F) +
  labs(x="Log(seqs*length)", y="Time", subtitle = "IQtree runtimes") + 
  scale_color_manual(values = palettes_model) +
  facet_grid(~method) +
  theme(legend.position = "bottom")

final_plot <- plot_homo + plot_aln + plot_fs + plot_iqtree + 
  plot_layout(widths = c(1, 3)) + 
  plot_annotation(tag_levels = c("A"))

ggsave(snakemake@output[[1]], final_plot, width = 10, height = 8, dpi = 300)

