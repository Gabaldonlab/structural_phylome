suppressPackageStartupMessages(library(tidyverse))
source("workflow/scripts/functions.R")

ids <- read_delim(snakemake@input[["seed_file"]], show_col_types = FALSE,
                  delim = "\t", col_names = c("group", "target"))
ids_todo <- unique(ids$group)
oprefix <- dirname(dirname(snakemake@output[[1]][1]))

if (ncol(ids)>1){
  print("The orthogroups are already set: only union set will be computed")
  for (id in ids_todo){
    hits <- filter(ids, group==id) %>% 
      pull(target)
    writeLines(hits, paste0(oprefix, "/", id, "/", id, "_common.ids"))
    } 
  } else {
  print("proceding with phylomeDB approach")
  blast <- read_delim(snakemake@input[["blast"]],
                      show_col_types=FALSE, col_names = blast_columns) %>% 
    mutate(type="blast")
  
  fs <- read_delim(snakemake@input[["fs"]],
                      show_col_types=FALSE, col_names = fs_columns) %>% 
    mutate(type="fs")
  
  df <- full_join(select(blast, query, target, type), select(fs, query, target, type),
            by = c("query", "target"), suffix=c("_blast", "_fs")) %>% 
    mutate(singleton = case_when(is.na(type_fs) ~ "only_blast",
                          is.na(type_blast) ~ "only_fs", .default = "common"))

  for (id in ids_todo){
    blast_hits <- filter(df, query==id) %>% 
      filter(type_blast=="blast", target!=id) %>% 
      slice_head(n = 150) %>% 
      pull(target)
    blast_hits <- c(id, blast_hits)
    fs_hits <- filter(df, query==id) %>% 
      filter(type_fs=="fs", target!=id) %>% 
      slice_head(n = snakemake@params[["max_seqs"]]-1) %>% 
      pull(target)
    fs_hits <- c(id, fs_hits)

    if ("blast" %in% snakemake@params[["modes"]]) {
      writeLines(blast_hits, paste0(oprefix,  "/", id, "/", id, "_blast.ids"))
    }
    if ("fs" %in% snakemake@params[["modes"]]) {
      writeLines(fs_hits, paste0(oprefix,  "/", id, "/", id, "_fs.ids"))
    }
    if ("common" %in% snakemake@params[["modes"]]) {
      common <- intersect(blast_hits, fs_hits)
      writeLines(common, paste0(oprefix,  "/", id, "/", id, "_common.ids"))
    }
    if ("union" %in% snakemake@params[["modes"]]) {
      union <- union(blast_hits, fs_hits)
      writeLines(union, paste0(oprefix,  "/", id, "/", id, "_union.ids"))
    }
  }
}

