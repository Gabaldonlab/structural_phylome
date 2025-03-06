suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(patchwork))
# suppressPackageStartupMessages(library(ggpubr))
source("workflow/scripts/functions.R")
theme_set(theme_classic())

# read species tree
sptree <- read.tree(snakemake@input[["sptree"]])

taxidmap <- read_delim(snakemake@input[["taxidmap"]], 
                       col_names = c("target", "Tax_ID"), show_col_types = FALSE)

# trees Dataframe
df_trees <- read_delim(snakemake@input[["trees"]], show_col_types = FALSE,
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
# disco_fls <- list.files("results/phylogeny/eggnog_benchmark/reco/disco/", full.names = T, pattern = "out*")

# Astral-pro results
clades <- read_delim(snakemake@input[["table"]], show_col_types = FALSE)

apro_files <- snakemake@input[["apro_trees"]]
# apro_files <- list.files("results/phylogeny/eggnog_benchmark/reco/apro/", full.names = T, pattern = "supp*")
apro_files <- apro_files[sapply(apro_files, function(x) file.size(x))>0]

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
  labs(x="", y="Quartet support", fill="Target sets")


disco_rf <- NULL
for (file in disco_fls){
  a <- read.tree(file)
  rf <- TreeDist::RobinsonFoulds(sptree, a, normalize = T)
  ntips <- sapply(a, function(x) length(x$tip.label))
  disco_rf <- bind_rows(disco_rf, 
                        tibble(ntips=ntips, rf=rf, bn=rep(basename(file), length(rf))))
}

pp1_plot <- nodes_df %>% 
    ggplot(aes(q1, label_int, fill=model, color=model, group=model)) +
    geom_linerange(aes(x=x, xmin=xmin, xmax=xmax, color=model),
                   data = . %>% 
                     group_by(model, label_int) %>% 
                     summarise(x=mean(q1), xmin=min(q1), xmax=max(q1)),
                   position = position_dodge(width=.7), color="grey70"
                   ) +
    geom_point(aes(shape=target),
               position = position_dodge(width=.7)) +
  scale_fill_manual(values=palettes_model) +
  scale_color_manual(values = palettes_model_d) +
  scale_shape_manual(values = c(21,22,23,24)) +
  guides(color = guide_legend(override.aes = list(pch = 21, size=3)),
         shape = guide_legend(override.aes = list(size=3))) +
  theme(legend.position = "bottom", legend.box="vertical")

sptree_plot <- ggtree(sptree) + 
  geom_tiplab() + 
  geom_nodelab(color="firebrick", 
               nudge_x = -.1, hjust = 1, 
               nudge_y = .2, size=2)+
  scale_x_continuous(expand = c(.2,.2))
  
apro_plot <- (sptree_plot | pp1_plot) + plot_layout(widths = c(1,1.5))
ggsave(snakemake@output[["apro"]], apro_plot, width = 12, height = 7)

# Read all trees
ts <- read.tree(text = df_trees$tree)
names(ts) <- paste(df_trees$gene, df_trees$target, df_trees$alphabet, df_trees$model, sep = "_")

# Compute var root2tip distance
var_df <- sapply(ts, function(x) 
  var(adephylo::distRoot(phytools::midpoint_root(x)))) %>% 
  enframe() %>%
  separate(name, c("gene", "target", "alphabet", "model")) %>%
  rename("variance_r2t"="value")

n_tips <- sapply(ts, function(x) length(x$tip.label)) %>% 
  enframe() %>%
  separate(name, c("gene", "target", "alphabet", "model")) %>%
  rename("n_tips"="value")

# n_sps <- sapply(ts, function(x) 
#   length(unique(taxidmap[taxidmap$target %in% x$tip.label, ]$Tax_ID))) %>% 
#   enframe() %>%
#   separate(name, c("gene", "target", "alphabet", "model")) %>%
#   rename("n_sps"="value")

# get BS distribution
df_bs <- fortify(ts) %>%
  filter(!isTip) %>%
  mutate(support=as.numeric(label)) %>%
  separate(.id, c("gene", "target", "alphabet", "model")) %>%
  filter(!is.na(support)) %>% 
  mutate(model=factor(model, levels=models))

# plot bootstrap
plot_bs <- ggplot(df_bs, aes(model, support, fill=target)) +
  # geom_violin(position = position_dodge(width = 1)) +
  geom_boxplot() +
  # ggdist::stat_pointinterval(position = position_dodge(width = 1), color="black") +
  scale_fill_manual(values = palettes_method) +
  theme_classic() +
  labs(y="Boostrap", x="") +
  theme(legend.position = "bottom")

# First of all let's check if they have the same heterogeneity model
plot_rate <- df_ml %>%
  # select(-tree) %>%
  mutate(Freq=grepl("F", rate), Inv=grepl("I", rate),
         Gamma=grepl("G", rate), FreeRate=grepl("R", rate)) %>%
  select(gene, targets, model, Freq, Inv, Gamma, FreeRate) %>%
  pivot_longer(cols = c(Freq,Inv,Gamma,FreeRate)) %>%
  inner_join(.,., by=c("gene", "targets", "name"), relationship = "many-to-many") %>%
  filter(as.numeric(model.x) > as.numeric(model.y)) %>%
  mutate(comp=paste0(model.x, " vs ", model.y)) %>%
  group_by(name, model.x, model.y) %>%
  count(value.x, value.y) %>%
  ggplot(aes(value.x, value.y, fill=n)) +
  geom_tile() +
  geom_text(aes(label=n), color="white") +
  scale_fill_viridis_c(option = "C") +
  ggh4x::facet_nested(name~model.x+model.y) +
  coord_fixed(expand = F) +
  labs(x="Is same?") +
  theme(legend.position = "none", axis.title.y = element_blank())


# Compare trees with RF distance
rf_df <- filter(df_trees, target=="common") %>% 
  select(-alphabet) %>% 
  inner_join(x=., y = ., by = c("gene", "target")) %>%
  filter(!is.na(tree.x), !is.na(tree.y)) %>%
  rowwise() %>%
  mutate(RF = phangorn::RF.dist(normalize = TRUE, check.labels = TRUE,
                                read.tree(text = tree.x),
                                read.tree(text = tree.y))) %>% 
         # n_tips = length(read.tree(text = tree.x)$tip.label)) %>%
  select(-tree.x, -tree.y)


rf_text <- group_by(rf_df, model.x, model.y) %>% 
  filter(as.integer(model.x)<=as.integer(model.y)) %>% 
  summarise(value=median(RF, na.rm=T))

# plot the heatmap of dist distribution and median value
rf_plot <- rf_df %>% 
  ggplot(aes(RF)) +
  geom_tile(aes(x=0.5, y=0.5, fill=value), data = rf_text) +
  geom_text(aes(x=0.5, y=.5, label=round(value, 2)), color="white", data = rf_text) +
  geom_rug(linewidth=.1, data= . %>% filter(as.integer(model.x)>as.integer(model.y))) +
  geom_density(aes(y = after_stat(scaled)), data= . %>% filter(as.integer(model.x)>as.integer(model.y))) + 
  scale_fill_distiller(palette="YlGnBu", direction = 1, limits = c(0,1)) + 
  scale_x_continuous(breaks = c(0,.5,1)) +
  coord_fixed(xlim = c(0,1), expand=0) +
  # facet_grid(method~method2, switch = "both") +
  facet_grid(model.x~model.y) +
  # ggh4x::facet_nested(targets+model~targets2+model2, nest_line = element_line(linetype = 1)) +
  labs(x="RF Distance") +
  theme_classic() +
  theme(legend.position = "none", 
        panel.background = element_rect(fill='transparent'),
        strip.placement = "outside",
        strip.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # strip.text = element_text(size=6),
        axis.text.x = element_text(size=6),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())


# Reconciliation notung
reco <- read_delim("results/phylogeny/eggnog_benchmark/reco/UP000005640_notung.tsv", 
                     show_col_types = FALSE, col_names = c("gene", "dups", "losses")) %>% 
  separate(gene, c("gene", "target", "alphabet", "model")) %>% 
  left_join(var_df) %>% 
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


plot_varr2t <- reco %>% 
  # filter(model %in% c("3Di", "LG", "GTR")) %>% 
 ggplot(aes(variance_r2t, model, fill=model)) +
  geom_violin(scale = "width", color="black") +
  scale_x_continuous(limits = c(0,.5), expand = expansion(mult=0, add = c(0,0.05))) +
  labs(x="Variance Root-to-tip distance") +
  scale_fill_manual(values = palettes_model) +
  theme(legend.position = 'none',
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,0,0,0))


reco_plot <- (plot_bs|plot_DL)/(plot_rf|plot_apro)

model <- ((plot_rate/rf_plot) | plot_varr2t) +
  plot_layout(heights = c(1,4,1))

ggsave(snakemake@output[["reco"]], reco_plot, width = 10, height = 8)
ggsave(snakemake@output[["model"]], model, width = 15, height = 12)

# plot_3di <- filter(df_ml, model!="LG") %>% 
#   pivot_wider(id_cols = c(gene, targets), values_from = c(rate,BIC), names_from = model) %>% 
#   # print() %>% 
#   mutate(diff_BIC=BIC_3Di-BIC_GTR, better_3Di=diff_BIC<0) %>% 
#   ggplot(aes(diff_BIC)) +
#   # geom_density() + 
#   ggdist::stat_slab(aes(fill = after_stat(level)), .width = c(.5, .66, .95, 1)) +
#   scale_fill_brewer() +
#   # geom_text(data = . %>%
#   #             mutate(x=quantile(diff_BIC, .1)) %>% 
#   #             group_by(targets, better_3Di, x) %>% 
#   #             summarise(n=n()) %>% 
#   #             pivot_wider(names_from = better_3Di, values_from = n) %>% 
#   #             mutate(prop_better=`TRUE`/(`TRUE`+`FALSE`)), 
#   #           mapping = aes(x=x, y=0.5, 
#   #                         label=paste0(round(prop_better*100, 2),"%"))) +
#   geom_vline(xintercept = 0) + 
#   coord_cartesian(expand = 0) + 
#   facet_grid(~targets, scales="free_x") +
#   labs(x="BIC 3Di - BIC GTR", y="Density", subtitle = "GTR vs 3Di BIC") + 
#   theme(legend.position = "bottom")

