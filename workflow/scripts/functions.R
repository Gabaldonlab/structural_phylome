interesting_types <- c("Topological domain", "Domain", 
                       "Region", "Helix", "Zinc finger",
                       "Beta strand", "Turn", "Coiled coil", "Transmembrane")

blast_columns <- c("query","target","pident","length","mismatch","gapopen",
                   "qstart","qend","sstart","send","evalue","bitscore",
                   "qcov","qcovhsp","qlen","slen","staxids")

fs_columns <- c("query","target","pident","length","mismatch","gapopen",
                "qstart","qend","sstart","send","evalue","bitscore",
                "lddt","alntmscore","rmsd","prob","qcov","tcov")

models <- c("C10", "C60","LG", "FM", "FMG","QT", "part", "FT",
            "FTPY", "3Di", "GTR", "AF", "LLM")

palettes_model <- c("#c71818", "#c71818", "#BF360C", "#FFA000", "#FFA000", "#ffc400",
                    "purple", "#D4E157","#646B00",
                    "#75C323", "#40A79D", "#33ceff", "#3396ff")
names(palettes_model) <- models


df_model <- tibble(model = factor(models), 
       data = c(rep("AA", 6), "Mixed", rep("3Di", 6)),
       algorithm = c(rep("ML", 3), rep("Distance", 3), "ML", 
                     rep("Distance", 2), rep("ML", 4))) %>% 
  filter(!model %in% c("QT", "FTPY") )
df_model$data <- factor(df_model$data, levels = c("AA", "Mixed", "3Di"))
df_model$algorithm <- factor(df_model$algorithm, levels = c("ML", "Distance"))

palettes_model_d <- colorspace::darken(palettes_model, amount = 0.6)
names(palettes_model_d) <- models
# "#FF0000" "#00A08A" "#F2AD00" "#F98400" "#5BBCD6"

palettes_method <- c("#00A08A", "#d63a3a", "#ffaa0e", "#6f63bb")
names(palettes_method) <- c("common", "blast", "fs", "union")
palettes_method_d <- colorspace::darken(palettes_method, amount = 0.3)
names(palettes_method_d) <- c("common", "blast", "fs", "union")
palette_singleton <- c("#00A08A", "#d63a3a", "#ffaa0e", "#6f63bb")
names(palette_singleton) <- c("common", "only_blast", "only_fs", "union")

palette_singleton_p <- c("#00A08A", "#d63a3a", "#ffaa0e", "#6f63bb")
names(palette_singleton_p) <- c("Common", "only Blast", "only Foldseek", "Union")

palette_singleton_s <- c("#00A08A", "#d63a3a", "#ffaa0e", "#6f63bb")
names(palette_singleton_s) <- c("Common", "Only Bp", "Only Fs", "Union")

palette_brh <- c("#00A08A","#d63a3a", "#ffaa0e")
names(palette_brh) <- c("common","blast_brh", "fs_brh")

palettes_alphabet <- c("#327AD9", "#BF3475", "#F29F05")
names(palettes_alphabet) <- c("3Di", "aa", "fident")
wespal <- wesanderson::wes_palette("Zissou1", 100, type = "continuous")

palettes_step <- c( "#F0E442", "#009E73", "#CC79A7")
names(palettes_step) <- c("foldtree", "iqtree", "quicktree")

palette_domains <- c("#31b7bcff", "#95c11fff")
names(palette_domains) <- c("Pfam", "3d")


add_marginal <- function(df, plot, x, y, fill, palette){
  xdens <- axis_canvas(plot, axis = "x")+
    geom_density(data = df, aes(x = .data[[x]], fill = .data[[fill]]),
                 alpha = 0.7, size = 0.2) +
    coord_cartesian(expand=0) + 
    scale_fill_manual(values = palette)
  ydens <- axis_canvas(plot, axis = "y", coord_flip = TRUE) +
    geom_density(data = df, aes(x = .data[[y]], fill = .data[[fill]]),
                 alpha = 0.7, size = 0.2)+
    coord_flip(expand = 0) +
    scale_fill_manual(values = palette)
  
  p1 <- insert_xaxis_grob(plot, xdens, grid::unit(.2, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  ggdraw(p2)
}


# Function to find identical columns and return their indexes
find_identical_columns <- function(matrix1, matrix2) {
  identical_columns <- apply(matrix1, 2, function(col1) any(apply(matrix2, 2, function(col2) all(col1 == col2))))
  return(which(identical_columns))
}


get_idxs_df <- function(seqs, aln, trimmed) {
  kept_columns <- find_identical_columns(as.matrix(aln@unmasked), as.matrix(trimmed@unmasked)) - 1
  df_idxs <- NULL
  for (seq in names(seqs)) {
    str_untrimmed <- as.character(aln@unmasked[[seq]])
    gaps <- stringr::str_locate_all(str_untrimmed, '-')[[1]]
    idxs <- c()
    aln_idxs <- c()
    
    for (idx in 1:nchar(str_untrimmed)){
      if (!idx %in% gaps[,1]) {
        to_remove <- sum(gaps[,1]<idx)
        aln_idxs <- c(aln_idxs, idx)
        idxs <- c(idxs, idx-to_remove)
      }
    }
    
    df_seed <- tibble(label=seq, pos=idxs, pos_aln=aln_idxs) %>% 
      left_join(tibble(pos_aln=kept_columns, pos_trim=1:length(kept_columns)), by="pos_aln") %>% 
      mutate(is_in_trimmed=pos_aln %in% kept_columns)
    
    df_idxs <- bind_rows(df_idxs, df_seed)
  }
  return(df_idxs)
}


plot_dom <- function(tree, seed, gff, prot, singleton_df, title="") {

  tree <- ape::root(tree, paste0("AF-", seed, "-F1"), resolve.root=TRUE)
  tree$edge.length <- (tree$edge.length/max(diag(ape::vcv(tree))))

  lbls <- fortify(tree, ladderize = FALSE) %>% 
    filter(isTip) %>% 
    arrange(y) %>% 
    pull(label)

  a <- gff %>%
    filter(seqid %in% gsub("AF-|-F1", "", tree$tip.label), type=="Domain") %>% 
    mutate(Note=gsub(";..*", "", gsub("Note=", "", rest)),
           seqid=factor(seqid, levels=gsub("AF-|-F1", "", lbls), ordered = TRUE)) %>% 
    ggplot(aes(x=start, xend=end, y=seqid, yend=seqid)) + 
    geom_segment(aes(x=1, xend=width,
                     y=seqid, yend=seqid), 
                 data = prot %>% mutate(seqid=factor(seqid, levels=gsub("AF-|-F1", "", lbls), ordered = TRUE))) +
    geom_segment(aes(x=start, xend=end,
                     y=seqid, yend=seqid, color=Note),  
                 linewidth = 2) + 
    geom_text(aes(label = Note, x=(start+end)/2, y=seqid), size=2) +
    theme_void() +
    labs(title = title) +
    theme(legend.position = "none")
  
  b <- fortify(tree, ladderize = FALSE) %>% 
    left_join(singleton_df, by="label") %>% 
    mutate(singleton = ifelse(label==paste0("AF-", seed, "-F1"), "seed", singleton)) %>%
    ggtree(ladderize = FALSE) + 
    # geom_tiplab() +
    geom_tippoint(aes(color=singleton)) +
    scale_color_manual(values = c(seed="black", palette_singleton)) + 
    # scale_x_continuous(limits = c(0,3)) +
    theme(legend.position = "none") 
  
  return (list(a, b))
}


get_disco_rf <- function(disco_fls, sptree) {
  disco_rf <- NULL
  for (file in disco_fls){
    a <- read.tree(file)
    rf <- TreeDist::RobinsonFoulds(sptree, a, normalize = T)
    nsps <- sapply(a, function(x) length(unique(x$tip.label)))
    disco_rf <- bind_rows(disco_rf, 
                          tibble(nsps=nsps, rf=rf, bn=rep(basename(file), length(rf))))
  }
  disco_rf <- disco_rf %>% 
    mutate(bn = gsub("disco_", "", gsub(".nwk", "", bn))) %>% 
    separate(bn, into = c("targets", "alphabet", "model"), sep = "_")
  return(disco_rf)  
}


get_bs_df <- function(trees) {
  df <- sapply(trees, function(x) as.numeric(x$node.label)) %>% 
    enframe(value = "support") %>% 
    unnest(cols = c(support)) %>% 
    separate(name, c("gene", "target", "alphabet", "model"), sep = "_") %>% 
    filter(!is.na(support))
  return(df)
}

# tree_stats <- function(trees, taxidmap) {
#   taxid_lookup <- setNames(taxidmap$Tax_ID, taxidmap$target)
#   df <- lapply(trees, function(x) {
#     c(n_tips = length(x$tip.label),
#       n_taxa = length(unique(taxid_lookup[x$tip.label])),
#       root_var = var(adephylo::distRoot(phytools::midpoint_root(x))) )
#   }) %>% 
#     do.call(rbind, .) %>%  # Combine results into a data frame
#     as.data.frame() %>% 
#     rownames_to_column("name") %>% 
#     separate(name, c("gene", "target", "alphabet", "model"), sep = "_") 
#   return(df)
# }



get_apro_stats <- function(apro_trees) {
  trees_df <- ggtree::fortify(apro_trees) %>% 
    group_by(.id) %>% 
    mutate(ordered = rank(y)) %>% 
    separate(.id, c("targets", "model")) %>% 
    filter(!isTip, label!="") %>% 
    mutate(label = gsub("\\[|\\]|\\'", "", str_replace_all(label, "[a-z]{1,2}[0-9]=", ""))) %>% 
    separate(label, c("pp1", "pp2", "pp3", "f1", "f2", "f3", "q1", "q2", "q3"), 
             ";", convert = TRUE) %>% 
    mutate(freq=f1+f2+f3, 
           model=factor(model, levels=models))
  return(trees_df)
}

get_rf_df <- function(df) {
  rf_df <- filter(df, target=="common") %>% 
    select(-alphabet) %>% 
    inner_join(x=., y = ., by = c("gene", "target")) %>%
    filter(!is.na(tree.x), !is.na(tree.y)) %>%
    rowwise() %>%
    mutate(RF =  TreeDist::RobinsonFoulds(read.tree(text = tree.x),  
                                          read.tree(text=tree.y),
                                          normalize = T)) %>%
    select(-tree.x, -tree.y)
  return(rf_df)
}


distance_to_seed <- function(df_trees) {
  ts <- read.tree(text = df_trees$tree)
  names(ts) <- paste0(df_trees$gene, "_", df_trees$model)

  process_row <- function(idx, ts) {
    nm <- str_split(names(ts)[idx], pattern = "_", simplify = T)
    seed <- nm[1]
    if (seed %in% ts[[idx]]$tip.label) {
      nn <- as.matrix(adephylo::distTips(ts[[idx]], method = "nNodes"))[seed, ]
      enframe(nn, name = "target", value = "nn_dist") %>%
        mutate(query = seed, model = nm[2])  
    }
  }
  df <- sapply(seq_along(ts), function(x) process_row(x, ts), simplify = FALSE) %>% 
    bind_rows()
  return(df)
}


stat_test <- function(df, alt="g", paired=T) {
  df <- droplevels(df)

  stat.test <- df %>% 
    t_test(y ~ model, p.adjust.method = "holm", 
           paired=paired,
           alternative=alt, ref.group = "LG") %>%
    add_xy_position(x = "model", fun = "max") %>% 
    left_join(df_model, by = c("group2"="model")) %>% 
    mutate(group2=paste0(group2,".",data))
    # mutate(model=factor(group2, levels = models)) %>% 
    # arrange(model) %>% 
  return(stat.test)
}
