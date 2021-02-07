
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
func_dir <- paste0(project_dir, "scripts/functions/")
ref_dir <- paste0(project_dir, "refs/")

RNA_name <- "hgsoc_repeats/RNA-seq-final"
RNA_dir <- paste0(home_dir, "projects/", RNA_name, "/")
DE_dir <- paste0(RNA_dir, "results/DE/without_primary_ascites/site/",
  "primary_vs_FT/tables/")

out_dir <- paste0(results_dir, "bp_assoc_RT/")
Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))

sig_p <- 0.05


########################################################################
### 0.  Load packages ###
########################################################################

library(GenomicRanges)
library(tibble)
library(dplyr)


########################################################################
### 1. Load retrotransposon breakpoint-assoc genes and top DE genes ###
########################################################################

RT_bp_genes <- readRDS(
  paste0(Robject_dir, "retrotransposon_breakpoint_assoc_genes.Rdata"),
)

up_genes <- read.table(
  paste0(DE_dir, "primary_vs_FT_upregulated_non_repeat.txt"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
down_genes <- read.table(
  paste0(DE_dir, "primary_vs_FT_downregulated_non_repeat.txt"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
DE_genes <- rbind(
  up_genes[up_genes$FDR < 0.05,],
  down_genes[down_genes$FDR < 0.05,]
)
DE_genes <- subset(DE_genes, select=c(logFC, FDR, symbol))


########################################################################
### 2. Find retrotransposon breakpoint-assoc genes which are DE ###
########################################################################

DE_RT_bp_genes <- lapply(RT_bp_genes, function(x) {

  DE <- x[x$symbol %in% DE_genes$symbol]
  mcols(DE) <- merge(mcols(DE), DE_genes, by="symbol")

  return(DE)

})
DE_RT_bp_genes <- unlist(as(DE_RT_bp_genes, "GRangesList"))
DE_RT_bp_genes$ensembl_id <- gsub("^.*\\.", "", names(DE_RT_bp_genes))
names(DE_RT_bp_genes) <- NULL

# remove unecessary cols:
mcols(DE_RT_bp_genes) <- mcols(DE_RT_bp_genes)[
  ,colSums(is.na(mcols(DE_RT_bp_genes))) < nrow(mcols(DE_RT_bp_genes))
]

write.table(
  as.data.frame(DE_RT_bp_genes),
  paste0(table_dir, "retrotransposon_breakpoint_assoc_DE_genes.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

gene_symbols <- unique(DE_RT_bp_genes$symbol)


########################################################################
### 3. Find retrotransposon breakpoint-assoc genes which are DE and 
# reported as cancer genes in OncoKB ###
########################################################################

oncokb <- read.table(
  paste0(ref_dir, "oncokb_cancer_genes.txt"),
  sep = "\t",
  header = T,
  stringsAsFactors = F
)

oncokb_genes <- gene_symbols[gene_symbols %in% oncokb$symbol]

# add custom cancer genes:
cancer_genes <- c(
  oncokb_genes,
  "CXCR2", "TOX", "UBE2C", "FGF18"
)

cancer_assoc_gr <- DE_RT_bp_genes[
  DE_RT_bp_genes$symbol %in% cancer_genes
]

# remove NA columns:
mcols(cancer_assoc_gr) <- mcols(cancer_assoc_gr)[
  ,colSums(is.na(mcols(cancer_assoc_gr))) < nrow(mcols(cancer_assoc_gr))
]


########################################################################
### 4. Identify bp types and mates ###
########################################################################

# load all breakpoints:
all_bp <- readRDS(paste0(Robject_dir, "all_breakpoints.Rdata"))
all_bp <- unlist(as(all_bp, "GRangesList"))

# fetch info for bps in cancer_assoc_gr:
for (i in 1:2) {

  info_gr <- all_bp[
    all_bp$id %in% eval(
      parse(
        text = paste0("cancer_assoc_gr$bp", i, "_id")
      )
    ), 
  ]
  info_df <- data.frame(
    id = as.character(info_gr$id),
    info = info_gr$bp_info,
    mate_id = gsub(
      ";.*$", "",
      gsub("^.*MATEID", "MATEID", info_gr$bp_info)
    ),
    stringsAsFactors = F
  )

  # fetch SV class:
  temp_class <- strsplit(
    as.character(info_df$info), ";"
  )
  info_df$class <- unlist(lapply(temp_class, function(x) x[[3]]))
  info_df$class[grep("MAPQ", info_df$class)] <- NA

  # remove info column:
  info_df <- subset(info_df, select = -info)

  colnames(info_df) <- c(
    paste0("bp", i, "_id"), 
    paste0("bp", i, "_mate_id"),
    paste0("bp", i, "_class")
  )

  # merge with cancer_assoc_gr
  mcols(cancer_assoc_gr) <- merge(
    mcols(cancer_assoc_gr),
    info_df,
    by = paste0("bp", i, "_id"),
    all = T
  )

}

# format mate id column:
cancer_assoc_gr$bp1_mate_id <- gsub(
  "MATEID=", "", 
  cancer_assoc_gr$bp1_mate_id
)
cancer_assoc_gr$bp2_mate_id <- gsub(
  "MATEID=", "", 
  cancer_assoc_gr$bp2_mate_id
)

# load low confidence breakpoints:
low_conf <- readRDS(
  paste0(Robject_dir, "all_low_confidence_breakpoints.Rdata")
)

# isolate all those from samples of interest:
low_conf <- low_conf[names(low_conf) %in% cancer_assoc_gr$sample]
low_conf <- unlist(as(low_conf, "GRangesList"))

total_bp <- c(all_bp, low_conf)
total_bp <- total_bp[!duplicated(total_bp)]
total_bp$id <- as.character(total_bp$id)
total_bp$bp_info <- as.character(total_bp$bp_info)
total_bp$joining_nt <- as.character(total_bp$joining_nt)
total_bp$join <- as.character(total_bp$join)

# find chr, position and class for mates of first bp:
m <- match(cancer_assoc_gr$bp1_mate_id, total_bp$id)

cancer_assoc_gr$bp1_mate_chr <- seqnames(total_bp)[m]
cancer_assoc_gr$bp1_mate_pos <- total_bp$bp[m]

bp1_mate_class <- gsub("^.*HOMSEQ", "HOMSEQ", total_bp$bp_info[m])
bp1_mate_class <- gsub("^.*INSERTION", "INSERTION", bp1_mate_class)
bp1_mate_class <- gsub("^.*IMPRECISE", "IMPRECISE", bp1_mate_class)
bp1_mate_class <- gsub(";.*$", "", bp1_mate_class)
cancer_assoc_gr$bp1_mate_class <- gsub("^.*MAPQ.*$", "NONE", bp1_mate_class)

# annotate breakpoint-mate joins:
cancer_assoc_gr$bp1_joining_nt <- total_bp$joining_nt[m]
cancer_assoc_gr$bp1_join <- total_bp$join[m]

# find chr, position and class for mates of second bp:
m <- match(cancer_assoc_gr$bp2_mate_id, total_bp$id)

cancer_assoc_gr$bp2_mate_chr <- as.character(seqnames(total_bp))[m]
cancer_assoc_gr$bp2_mate_pos <- total_bp$bp[m]

bp2_mate_class <- gsub("^.*HOMSEQ", "HOMSEQ", total_bp$bp_info[m])
bp2_mate_class <- gsub("^.*INSERTION", "INSERTION", bp2_mate_class)
bp2_mate_class <- gsub("^.*IMPRECISE", "IMPRECISE", bp2_mate_class)
bp2_mate_class <- gsub(";.*$", "", bp2_mate_class)
cancer_assoc_gr$bp2_mate_class[1:length(m)] <- gsub("^.*MAPQ.*$", "NONE", bp2_mate_class)

# annotate breakpoint-mate joins:
cancer_assoc_gr$bp2_joining_nt <- total_bp$joining_nt[m]
cancer_assoc_gr$bp2_join <- total_bp$join[m]

# reorder mcols:
cancer_assoc <- data.frame(
  sample = cancer_assoc_gr$sample,
  symbol = cancer_assoc_gr$symbol,
  ensembl_id = cancer_assoc_gr$ensembl_id,
  logFC = cancer_assoc_gr$logFC,
  FDR = cancer_assoc_gr$FDR,
  chr = seqnames(cancer_assoc_gr),
  bp1_id = cancer_assoc_gr$bp1_id,
  bp1_pos = cancer_assoc_gr$bp1_pos,
  bp1_class = cancer_assoc_gr$bp1_class,
  bp1_mate_id = cancer_assoc_gr$bp1_mate_id,
  bp1_mate_chr = cancer_assoc_gr$bp1_mate_chr,
  bp1_mate_pos = cancer_assoc_gr$bp1_mate_pos,
  bp1_mate_class = cancer_assoc_gr$bp1_mate_class,
  bp1_joining_nt = cancer_assoc_gr$bp1_joining_nt,
  bp1_join = cancer_assoc_gr$bp1_join,
  bp2_id = cancer_assoc_gr$bp2_id,
  bp2_pos = cancer_assoc_gr$bp2_pos,
  bp2_class = cancer_assoc_gr$bp2_class,
  bp2_mate_id = cancer_assoc_gr$bp2_mate_id,
  bp2_mate_chr = cancer_assoc_gr$bp2_mate_chr,
  bp2_mate_pos = cancer_assoc_gr$bp2_mate_pos,
  bp2_mate_class = cancer_assoc_gr$bp2_mate_class,
  bp2_joining_nt = cancer_assoc_gr$bp2_joining_nt,
  bp2_join = cancer_assoc_gr$bp2_join
)

write.table(
  cancer_assoc,
  paste0(table_dir, "final_cancer_gene_and_retrotransposon_assoc_breakpoints.txt"),
  sep = "\t",
  quote = F,
  col.names = T,
  row.names = F
)

save.image(paste0(Robject_dir, "identified_RT_bp_assoc_DE_genes_image.Rdata"))

