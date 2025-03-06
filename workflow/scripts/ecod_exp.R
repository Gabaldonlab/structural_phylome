ecod <- read_delim("~/Downloads/ecod.develop292.domains.txt", 
                   comment = "#", delim = " ", col_names = FALSE)
df_filt$target
indb <- gsub("AF-", "", gsub("-F1", "_F1", unique(names(db_fasta))))
length(intersect(indb, ecod$X1))/length(indb)
