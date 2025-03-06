suppressPackageStartupMessages(library(tidyverse))
library(patchwork)
suppressPackageStartupMessages(library(cowplot))
source("workflow/scripts/functions.R")
theme_set(theme_bw())

taxidmap <- read_delim(snakemake@input[["taxidmap"]], 
                       col_names = c("target", "Tax_ID"), show_col_types=FALSE)

# tax <- read_delim(snakemake@input[["groups"]], show_col_types=F)
# # read input table
# lvls_tax <- rev(c(unique(tax$Clade), "archaea", "bacteria"))

# table_columns <- c('Proteome_ID', 'Tax_ID', 'count1', 'count2', 'count3', 
#                    'genome', 'source', 'species', 'mnemo')

table <- read_delim(snakemake@input[["table"]], show_col_types=FALSE, 
                    delim = "\t") 
  # left_join(tax, by = "Tax_ID") %>% 
  # mutate(Clade = ifelse(is.na(Clade), mnemo, Clade),
  #        Clade = factor(Clade, levels = lvls_tax))

blast_self <- read_delim(snakemake@input[["self_blast"]], show_col_types=FALSE, 
                        col_names = blast_columns) %>% 
  filter(query==target) %>% 
  mutate(blast_selfbit=bitscore) %>% 
  select(query, blast_selfbit)

fs_self <- read_delim(snakemake@input[["self_fs"]], show_col_types=FALSE, 
                      col_names = fs_columns) %>% 
  filter(query==target) %>% 
  mutate(fs_selfbit=bitscore) %>% 
  select(query, fs_selfbit)

self <- full_join(blast_self, fs_self)

fs <- read_delim(snakemake@input[["fs"]], col_names = fs_columns, show_col_types = FALSE) %>% 
  mutate(evalue = ifelse(evalue==0, 1e-180, evalue)) %>%
  filter(query!=target) %>% 
  left_join(taxidmap)
blast <- read_delim(snakemake@input[["blast"]], col_names = blast_columns, show_col_types = FALSE) %>% 
  mutate(evalue = ifelse(evalue==0, 1e-180, evalue)) %>%
  filter(query!=target)

#queries <- sample(unique(blast$query), 200)
#fs <- filter(fs, query %in% queries) 
#blast <- filter(blast, query %in% queries)

evals <- c(Inf, 1e-2, 1e-3, 1e-5)
min_covs <- c(0, 30, 80)


# fs_brhs <- read_delim(snakemake@input[["fs_brh"]], col_names = c("query", "target"), show_col_types = FALSE) %>% 
#   mutate(method="blast_brh")
#   # group_by(query) %>%
#   # summarise(fs_brhs=list(target))

# blast_brhs <- read_delim(snakemake@input[["blast_brh"]], col_names = c("query", "target"), show_col_types = FALSE) %>% 
#   mutate(method="fs_brh")
#   # group_by(query) %>%
#   # summarise(blast_brhs=list(target))

# brhs <- rbind(blast_brhs, fs_brhs) %>% 
#   group_by(query, target, method) %>% 
#   mutate(n=1) %>% 
#   pivot_wider(names_from = method, values_fill = 0, values_from = n) %>% 
#   mutate(common = ifelse(fs_brh+blast_brh==2, 1, 0))

# grouped_brhs <- brhs %>%
#   group_by(query) %>% 
#   summarise(fs_brh = sum(fs_brh), blast_brh=sum(blast_brh), common_brh=sum(common))

df <- full_join(select(blast, query, target, pident, evalue, bitscore, qcov),
                select(fs, query, target, pident, evalue, bitscore, qcov),
                by = c("query", "target"), suffix=c("_blast", "_fs")) %>%
  left_join(self) %>%
  mutate(pident_fs = 100*pident_fs, 
         qcov_fs = 100*qcov_fs,
         singleton = case_when(is.na(evalue_fs) ~ "only_blast",
                               is.na(evalue_blast) ~ "only_fs", .default = "common"),
         blast_norm = bitscore_blast/blast_selfbit, fs_norm = bitscore_fs/fs_selfbit) %>% 
  select(query, target, pident_blast, pident_fs,  
         qcov_blast, qcov_fs, evalue_blast, evalue_fs, blast_norm, fs_norm, singleton)


norm_bit_p <- ggplot(df, aes(blast_norm, fs_norm)) + 
  stat_bin_2d(breaks=seq(0, 1,0.02)) + 
  scale_fill_gradientn(colors=wespal) +
  coord_fixed(expand = 0) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none", panel.grid = element_blank()) 
norm_bit_p <- add_marginal(df, norm_bit_p, "blast_norm", "fs_norm", "singleton", palette_singleton)

eval_p <- ggplot(df, aes(-log10(evalue_blast), -log10(evalue_fs))) + 
  stat_bin_2d(breaks=seq(0,20,0.4)) + 
  scale_fill_gradientn(colors=wespal) +
  geom_hline(yintercept = -log10(evals), linetype=2) + 
  geom_vline(xintercept = -log10(evals), linetype=2) + 
  coord_fixed(expand = 0) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0, 20)) +
  theme(legend.position = "none", panel.grid = element_blank()) 
eval_p <- add_marginal(mutate(df, evalue_blast=-log10(evalue_blast), evalue_fs=-log10(evalue_fs)), 
             eval_p, "evalue_blast", "evalue_fs", "singleton", palette_singleton)

pident_p <- ggplot(df, aes(pident_blast, pident_fs)) + 
  stat_bin_2d(breaks=seq(0,100,2)) + 
  scale_fill_gradientn(colors=wespal) +
  # scale_fill_viridis_c(option = "C") +
  coord_fixed(expand = 0) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none", panel.grid = element_blank())
pident_p <- add_marginal(df, pident_p, "pident_blast", "pident_fs", "singleton", palette_singleton)

qcov_p <- ggplot(df, aes(qcov_blast, qcov_fs)) + 
  stat_bin_2d(breaks=seq(0,100,2)) + 
  geom_hline(yintercept = min_covs, linetype=2) + 
  geom_vline(xintercept = min_covs, linetype=2) + 
  scale_x_continuous(limits = c(0, 100)) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_fill_gradientn(colors=wespal) +
  coord_fixed(expand = 0) +
  theme(legend.position = "none", panel.grid = element_blank()) 
qcov_p <- add_marginal(df, qcov_p, "qcov_blast", "qcov_fs", "singleton", palette_singleton)

df_overlap <- tibble()

for (cov in min_covs){
  for (eval1 in evals){
    for (eval2 in evals){
      
      print(paste(cov, "-", eval1, "-", eval2))

      a <- filter(blast, evalue<eval1, qcov >= cov)
      b <- filter(fs, evalue<eval2, qcov >= cov/100)

      to_add <- full_join(select(a, query, target, pident, evalue, bitscore, qcov),
                          select(b, query, target, pident, evalue, bitscore, qcov),
                          by = c("query", "target"), suffix=c("_blast", "_fs")) %>% 
        mutate(singleton = case_when(is.na(evalue_fs) ~ "only_blast",
                           is.na(evalue_blast) ~ "only_fs", .default = "common")) %>%
        group_by(query, singleton) %>% 
        count() %>% 
        pivot_wider(id_cols = query, names_from = singleton, 
                    values_from = n, values_fill = 0, names_prefix = "total_") %>% 
        mutate(total_common = ifelse("total_common" %in% names(.), total_common, 0),
               total_only_blast = ifelse("total_only_blast" %in% names(.), total_only_blast, 0),
               total_only_fs = ifelse("total_only_fs" %in% names(.), total_only_fs, 0)) %>% 
        mutate(eval_blast = factor(eval1), 
               eval_fs = factor(eval2), 
               cov = cov,
               total = total_common+total_only_blast+total_only_fs,
               jaccard = total_common/total,
               overlap_blast = total_common/(total_common+total_only_blast),
               overlap_fs = total_common/(total_common+total_only_fs)) %>% 
        select(query, eval_blast, eval_fs, cov, jaccard, overlap_blast, overlap_fs, total)

      df_overlap <- bind_rows(df_overlap, to_add)
    }
  }
}


plot_intersection <- df_overlap %>% 
  filter(total>2) %>% 
  select(eval_blast, eval_fs, cov, jaccard, overlap_blast, overlap_fs) %>% 
  pivot_longer(cols = !c(eval_blast, eval_fs, cov), names_to = "intersection") %>%
  mutate(cov=paste("cov:", cov), intersection = factor(intersection, levels = c("overlap_blast","jaccard","overlap_fs"))) %>% 
  group_by(eval_blast, eval_fs, intersection, cov) %>% 
  summarise(mean_overlap=mean(value, na.rm = T), sd_overlap=sd(value, na.rm = T)) %>% 
  ggplot(aes(eval_blast, eval_fs)) + 
  geom_tile(aes(fill=mean_overlap), color="black") +
  geom_text(aes(label=round(mean_overlap, 2)), color="white") +
  # geom_text(aes(label=paste(round(mean_overlap, 2), "±", round(sd_overlap, 2))), color="white") +
  scale_fill_viridis_c(limits = c(0,1), option = "B") +
  facet_grid(cov~intersection) +
  labs(x="Blast E-value", y="Foldseek E-value") + 
  coord_fixed(expand = 0) +
  theme_minimal() + 
  theme(legend.position = "none")


# taxonomy plot
# taxo_df <- tibble()

# for (eval in evals) {
#   for (cov in min_covs) {
#     print(paste(eval, "-",cov))
#     add_bs <- blast %>%
#       filter(evalue < eval, qcov > cov) %>% 
#       left_join(table, by = c("staxids"="Tax_ID")) %>% 
#       group_by(query, Clade) %>% 
#       count() %>% 
#       group_by(Clade) %>% 
#       count(n>0) %>% 
#       mutate(n=n/length(unique(blast$query)), evalue=eval, cov=cov, method="blast", brh="All")
#     add_bs_brh <- blast %>%
#       filter(evalue < eval, qcov > cov) %>% 
#       inner_join(blast_brhs, by = c("query", "target")) %>% 
#       left_join(table, by = c("staxids"="Tax_ID")) %>% 
#       group_by(query, Clade) %>% 
#       count() %>% 
#       group_by(Clade) %>% 
#       count(n>0) %>% 
#       mutate(n=n/length(unique(blast$query)), evalue=eval, cov=cov, method="blast", brh="BRH")
#     add_fs <- fs %>%
#       filter(evalue < eval, qcov > cov/100) %>% 
#       left_join(table, by="Tax_ID") %>% 
#       group_by(query, Clade) %>% 
#       count() %>% 
#       group_by(Clade) %>% 
#       count(n>0) %>% 
#       mutate(n=n/length(unique(fs$query)), evalue=eval, cov=cov, method="fs", brh="All")
#     add_fs_brh <- fs %>%
#       filter(evalue < eval, qcov > cov/100) %>% 
#       inner_join(fs_brhs, by = c("query", "target")) %>% 
#       left_join(table, by="Tax_ID") %>% 
#       group_by(query, Clade) %>% 
#       count() %>% 
#       group_by(Clade) %>% 
#       count(n>0) %>% 
#       mutate(n=n/length(unique(fs$query)), evalue=eval, cov=cov, method="fs", brh="BRH")
#     taxo_df <- bind_rows(taxo_df, add_bs, add_bs_brh, add_fs, add_fs_brh)
#   }
# }



# taxo_plot <- ggplot(taxo_df) + 
#   geom_tile(aes(y=Clade, method, fill=n)) + 
#   # facet_grid(evalue~cov) +
#   ggh4x::facet_nested(brh~evalue+cov, scale = "free_y") +
#   # scale_fill_gradientn(colours = wespal) + 
#   scale_fill_viridis_c(option = "D") +
#   coord_cartesian(expand=0)


final_plot <- (pident_p | qcov_p | eval_p | norm_bit_p) / (plot_intersection)# + taxo_plot)

ggsave(snakemake@output[["eda"]], final_plot, width = 14, height = 10, dpi = 300)


blast_p <- blast %>% 
  # filter(evalue<min_evalue) %>% 
  mutate('Inf' = evalue < Inf, 'e<1e-2' = evalue<1e-2, 
         'e<1e-3' = evalue<1e-3, 'e<1e-5' = evalue<1e-5) %>% 
  pivot_longer(cols = c('Inf', 'e<1e-2', 'e<1e-3', 'e<1e-5')) %>% 
  group_by(query, name) %>% 
  # mutate(included = ifelse(row_number()<max_seqs, TRUE, FALSE)) %>%
  summarise(n=sum(value)) %>% 
  ggplot(aes(n, color=name)) + 
  stat_ecdf() +
  scale_color_manual(values = viridisLite::cividis(4)) +
  labs(subtitle = "Blast #hits based on evalue filter") + 
  geom_vline(xintercept = snakemake@params[["max_seqs"]])

fs_p <- fs %>% 
  # filter(evalue<min_evalue) %>% 
  mutate('Inf' = evalue < Inf, 'e<1e-2' = evalue<1e-2, 
         'e<1e-3' = evalue<1e-3, 'e<1e-5' = evalue<1e-5) %>% 
  pivot_longer(cols = c('Inf', 'e<1e-2', 'e<1e-3', 'e<1e-5')) %>% 
  group_by(query, name) %>% 
  # mutate(included = ifelse(row_number()<max_seqs, TRUE, FALSE)) %>%
  summarise(n=sum(value)) %>% 
  ggplot(aes(n, color=name)) + 
  stat_ecdf() +
  scale_color_manual(values = viridisLite::cividis(4)) +
  labs(subtitle = "Foldseek #hits based on evalue filter") + 
  geom_vline(xintercept = snakemake@params[["max_seqs"]])


saturation_plot <- (blast_p | fs_p) + plot_layout(guides = "collect")
ggsave(snakemake@output[["saturation"]], saturation_plot, width = 7, height = 5, dpi = 300)

