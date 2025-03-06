a <- df %>% 
  # head(10000) %>% 
  mutate(qcov = log10(qcov_Blast/qcov_Foldseek),
         # evalue_ratio = evalue_Blast/evalue_Foldseek,
         pident = log10(pident_Blast/pident_Foldseek)) %>% 
  # filter(pident_ratio < 1)
  select(c("qcov", "pident")) %>% 
  pivot_longer(cols = everything())

ggplot(a, aes(value, color=name)) +
  geom_density() + 
  scale_x_continuous(limits = c(-1, 1)) +
  geom_vline(xintercept = 0) + 
  scale_color_manual(values = c("#F2CB05", "#0378A6")) + 
  labs(x="Log10(Blast/Foldseek)") +
  theme(legend.position = "bottom", legend.title = element_blank())


df_ml %>% 
  filter(model!="LG") %>% 
  # filter(gene=="AF-A0A075B6U4-F1") %>% 
  group_by(gene, targets) %>% 
  mutate(rank = factor(dense_rank(BIC))) %>% 
  select(BIC, model, rank) %>% 
  ungroup() %>% 
  count(model, rank) %>% 
  ggplot(aes(rank, model, fill=n)) + 
  # geom_point(aes(size=n), pch=21) +
  geom_tile() +
  geom_text(aes(label=n), color="white") +
  scale_fill_viridis_c() + 
  coord_fixed(expand = F) +
  labs(x="Rank", y="Model") + 
  theme(legend.position = "none")


dl <- read_delim("results/old/test/reco/UP000005640_DTL.tsv") %>% 
  rename("gene"="id", "target"="targets")

df_trees <- read_delim("results/old/test/trees/UP000005640_unrooted_trees.txt",
                       show_col_types = FALSE,
                       col_names = c("gene", "target", "alphabet", "model", "tree")) %>% 
  rowwise() %>% 
  mutate(model=factor(model, levels=models),
         n_tips = length(read.tree(text = tree)$tip.label)) 

a <- left_join(df_trees, dl) %>% 
  filter(!is.na(reco_score))

lm <- lm(dups+losses ~ n_tips, data = a)
a$res <- lm$residuals

tmp_plot <- ggplot(a, aes(y = res, x = n_tips, fill=target)) + 
  geom_point(pch=21) +
  facet_wrap(~target) + 
  geom_abline(slope = 0) + 
  scale_fill_manual(values = palettes_method) + 
  labs(x="Num. of Tips", 
       y="Residuals of: D+L ~ Num. of Tips") + 
  theme(legend.position = "none")

ggsave("draft/figures/review_hetero_TMP.png", tmp_plot,
       width = 4, height = 3)