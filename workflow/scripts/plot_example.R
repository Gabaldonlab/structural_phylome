suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggmsa))
library(patchwork)

source("workflow/scripts/functions.R")
theme_set(theme_classic())

seeds <- readLines(snakemake@input[["ids"]])

seed_sp <- snakemake@params[["seed"]]
seed_dir <- gsub("trees", "seeds", dirname(snakemake@input[["trees"]]))

taxmap <- read_delim(snakemake@input[["taxidmap"]], 
                     delim = "\t", col_names = c("label", "Tax_ID"),
                     show_col_types = FALSE) %>% 
  left_join(read_delim(snakemake@input[["table"]],
                       show_col_types = FALSE), by = "Tax_ID") %>% 
  dplyr::select(label, Proteome_ID, mnemo)

sptree <- read.tree(snakemake@input[["sptree"]])
translate <- deframe(distinct(dplyr::select(taxmap, Proteome_ID, mnemo)))
sptree$tip.label <- names(translate)[match(translate, sptree$tip.label)]

species_dist_df <- as.matrix(adephylo::distTips(sptree, method = "nNodes"))[seed_sp, ] %>% 
  enframe(name = "Proteome_ID", value = "dist")

gff <- read_delim(c(snakemake@input[["gff"]]), 
                  delim = "\t", comment = "#", 
                  col_names = c("seqid", "source", "type",
                                "start", "end", "score",
                                "strand", "phase", "rest"),
                  show_col_types = FALSE) %>% 
  filter(type %in% c(interesting_types, "Chain"))

reco <- read_delim(snakemake@input[["reco"]], 
                   show_col_types = FALSE, col_names = c("gene", "dups", "losses")) %>% 
  separate(gene, c("id", "targets", "alphabet", "model")) 

plot_list = list()
plot_domains = list()

for (seed in seeds) {
  print(seed)
  fls <- list.files(paste0(seed_dir, "/", seed_sp, "/", seed), full.names = TRUE)
  ids <- fls[grepl("(blast|fs).*ids$", fls)]

  x <- sapply(ids, readLines, simplify = FALSE)
  names(x) <- gsub(".*_|\\.ids", "", basename(ids))
  
  singleton_df <- list("common" = intersect(x$fs, x$blast), 
                 "only_fs" = setdiff(x$fs, x$blast), 
                 "only_blast" = setdiff(x$blast, x$fs)) %>% 
    enframe(name = "singleton", value = "label") %>% 
    unnest(cols = label)
  
  species_df <- singleton_df %>% 
    left_join(taxmap, by="label") %>% 
    left_join(species_dist_df, by="Proteome_ID") %>% 
    group_by(singleton, mnemo) %>% 
    summarise(furthest=max(dist))
  furthest <- species_df %>% slice_max(furthest,with_ties = FALSE)
  subtitle <- paste(length(unique(species_df$mnemo)), "species with max distance:", 
        paste0(furthest$singleton, ":", furthest$furthest, collapse = ", "))
  
  singleton_df_red <- tibble(targets= c("blast", "fs", "union", "common"), 
                             n=c(length(x$blast), length(x$fs), 
                                 length(union(x$blast, x$fs)),
                                 length(intersect(x$blast, x$fs))))

  conf_df <- tibble(label=union(x$blast, x$fs)) %>% 
    left_join(taxmap, by="label") %>% 
    mutate(confidence_file=paste0(snakemake@params[["struct_dir"]],
                                  Proteome_ID,
                                  "/confidence/",
                                  label,
                                  "-confidence_v4.json.gz")) %>% 
    rowwise() %>% 
    mutate(conf=list(jsonlite::fromJSON(gzfile(confidence_file))$confidenceScore))
  
  plot_venn <- ggvenn::ggvenn(x, auto_scale = T, set_name_size = 0,
                              fill_color = c("#F29F05", "#267365"))

  alns <- fls[grepl("union.*clean$", fls)]

  seqs <- Biostrings::readAAStringSet(fls[grepl("union.*_aa.seqs", fls)])
  length <- nchar(as.character(seqs[[paste0("AF-", seed, "-F1")]]))
  prot <- as.data.frame(seqs@ranges) %>% 
    mutate(seqid=gsub("AF-|-F1", "", names))
  
  plot_seed <- gff %>%
    filter(seqid==seed, type!="Chain") %>% 
    rowwise() %>% 
    mutate(Note=gsub(";..*", "", gsub("Note=", "", rest))) %>% 
    # filter(!type %in% c("Helix", "Beta strand", "Turn")) %>% 
    ggplot() + 
    geom_segment(aes(x=1, xend=width,
                     y=seqid, yend=seqid), data = filter(prot, seqid==seed)) +
    geom_segment(aes(x=start, xend=end,
                     y=seqid, yend=seqid, color=type), 
                 linewidth = 3) + 
   ggrepel::geom_text_repel(aes(label = Note, x=(start+end)/2, y=seqid), nudge_y = -.05) + 
    theme_void() + 
    theme(legend.position = "bottom")
  # 3Di aln
  aln_3di <- Biostrings::readAAMultipleAlignment(gsub(".clean$", "", alns[grepl("3Di.alg", alns)]))
  aln_3di_trimmed <- Biostrings::readAAMultipleAlignment(alns[grepl("3Di.alg", alns)])
  # aln_3di_mat <- as.matrix(aln_3di@unmasked)
  # pos_df_3di <- tibble(pos=1:ncol(aln_3di_mat), 
  #                      cons=nrow(aln_3di_mat) - colSums(aln_3di_mat=="-" | aln_3di_mat=="X"))
  
  idx_df <- get_idxs_df(seqs, aln_3di, aln_3di_trimmed)

  scaled_conf_df <- conf_df %>% 
    dplyr::select(-confidence_file) %>% 
    unnest(conf) %>% 
    group_by(label) %>% 
    mutate(pos=row_number()) %>% 
    left_join(idx_df, by=c("label", "pos")) %>% 
    filter(is_in_trimmed) %>% 
    group_by(pos_trim) %>% 
    summarise(mn=mean(conf)) %>% 
    mutate(scaled_mn=length(names(seqs))*mn/100)
  
  plot_3di_aln <- ggmsa(aln_3di_trimmed, seq_name = FALSE, 
                        font = NULL, show.legend = TRUE, border = NA) +
    labs(title = paste(seed, "- 3Di union hits cleaned aln")) + 
    geom_line(mapping = aes(x=pos_trim, y=scaled_mn), data = scaled_conf_df, inherit.aes = F) +
    coord_cartesian() + 
    # geom_seqlogo() + 
    theme_void() + 
    theme(axis.text.x = element_text())
  
  aln_aa_trimmed <- Biostrings::readAAMultipleAlignment(alns[grepl("aa.alg", alns)])
  # aln_aa_mat <- as.matrix(aln_aa@unmasked)
  # pos_df_aa <- tibble(pos=1:ncol(aln_aa_mat), 
  #                      cons=nrow(aln_aa_mat) - colSums(aln_aa_mat=="-" | aln_aa_mat=="X"))
  
  plot_aa_aln <- ggmsa(aln_aa_trimmed, seq_name = FALSE, font = NULL, show.legend = TRUE, border = NA) +
    labs(title = paste(seed, "- aa union hits cleaned aln")) + 
    # geom_line(mapping = aes(x=pos, y=cons), data = pos_df_aa, inherit.aes = F) +
    coord_cartesian() + 
    # geom_seqlogo() + 
    theme_void() + 
    theme(axis.text.x = element_text())

  trees_files <- fls[grepl(".nwk$", fls)]
  trees_files <- trees_files[!grepl("_rnm.nwk$", trees_files)]

  trees <- read.tree(text = sapply(trees_files, readLines))
  names(trees) <- gsub(paste0(seed,"_|(3Di|aa)_|.nwk$"), "", basename(trees_files))
  
  tree_plot <- fortify(trees) %>% 
    separate(.id, into = c("targets", "model")) %>% 
    left_join(singleton_df, by="label") %>% 
    mutate(singleton = ifelse(label==paste0("AF-", seed, "-F1"), "seed", singleton),
           seed = seed,
           model = factor(model, levels=models)) %>% 
    # left_join(scores, by=c("targets", "model", "seed"="id")) %>%
    ggtree() + 
    geom_tippoint(aes(color=singleton)) + 
    # geom_text(mapping = aes(label=label), 
    #           data = . %>% 
    #             group_by(targets, model) %>% 
    #             summarise(label=paste0(mean(dups), "|", mean(losses), "|", mean(score)),
    #                       y=sum(isTip)-1, x=mean(branch.length))) +
    facet_grid(targets~model, scales = "free") + 
    scale_color_manual(values = c(seed="black", palette_singleton))

  a <- plot_dom(trees[["union_LG"]], seed, gff, prot, singleton_df, "LG")
  b <- plot_dom(trees[["union_FT"]], seed, gff, prot, singleton_df, "FT")
  plot_tree_domain <- (a[[2]] + a[[1]] + b[[2]] + b[[1]]) + plot_layout(nrow = 1)

  plot_notung <- reco %>% 
    filter(id==seed) %>% 
    left_join(singleton_df_red) %>% 
    ggplot(aes(targets, (dups+losses)/n, color=model)) + 
    # geom_boxplot(show.legend = FALSE) +
    geom_point(shape="-", alpha=0, size=5) + 
    geom_text(aes(label=model, y=a), data = . %>% group_by(model, targets) %>% 
                  summarise(a=min((dups+losses)/n)), 
              position = position_dodge(width = .5)) +
    scale_color_manual(values = palettes_model) + 
    guides(color=guide_legend(nrow=2,byrow=TRUE, override.aes = list(alpha=1))) +
    theme(legend.position = "none", legend.justification = "left")

  plot_conf <- conf_df %>% 
    dplyr::select(-confidence_file) %>% 
    unnest(conf) %>% 
    group_by(label) %>% 
    mutate(pos=row_number()) %>% 
    left_join(idx_df, by=c("label", "pos")) %>% 
    left_join(singleton_df, by="label") %>% 
    ggplot(aes(is_in_trimmed, conf, fill=singleton)) + 
    geom_violin(alpha=.9) +     
    scale_y_continuous(limits = c(0,100)) + 
    labs(y="Mean LDDT per site") +
    scale_fill_manual(values = palette_singleton) + 
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    theme(legend.position = "bottom", legend.title=element_blank(),
          legend.key.size = unit(0.3, "cm"))

  title <- filter(gff, seqid==seed, type=="Chain") %>%
    mutate(Note=gsub(";..*", "", gsub("Note=", "", rest))) %>% 
    pull(Note) %>% simplify()
  title <- ifelse(is.null(title), seed, paste(seed, "-", title))

  first_row <- (plot_seed | plot_venn) + plot_layout(widths = c(4.5,1))
  bottom_row <- ((plot_conf / plot_notung) | tree_plot) + plot_layout(widths = c(1, 5))
  outplot <- (first_row / plot_3di_aln / plot_aa_aln / bottom_row) + 
    plot_layout(heights = c(.5,1,1,2)) +
    plot_annotation(title = title, subtitle = subtitle, tag_levels = "A")
  
  plot_list[[seed]] = outplot
  plot_domains[[seed]] = plot_tree_domain
}

pdf(snakemake@output[[1]], width = 16, height = 16)
for (i in names(plot_list)) {
  print(plot_list[[i]])
  print(plot_domains[[i]])
}
dev.off()
