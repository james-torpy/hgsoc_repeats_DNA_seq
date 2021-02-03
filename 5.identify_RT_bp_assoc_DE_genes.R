
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

write.table(
  as.data.frame(DE_RT_bp_genes),
  paste0(table_dir, "retrotransposon_breakpoint_assoc_DE_genes.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

final_gene_gr <- DE_RT_bp_genes
mcols(final_gene_gr) <- subset(
  mcols(final_gene_gr),
  select = c(symbol, sample, logFC, FDR)
)

gene_symbols <- unique(final_gene_gr$symbol)


########################################################################
### 3. Find retrotransposon breakpoint-assoc genes which are DE and 
# reported CNA-associated by TGCA (from 489 HGSOC samples) ###
########################################################################

TGCA_cna <- read.table(
  paste0(ref_dir, "TGCA/OV_CNV/CNA_genes.txt"),
  header = T,
  sep = "\t",
  stringsAsFactors = F
)
TGCA_cna_assoc <- gene_symbols[gene_symbols %in% TGCA_cna$Gene]
TGCA_cna_assoc_gr <- DE_RT_bp_genes[
  DE_RT_bp_genes$symbol %in% TGCA_cna_assoc
]

# add whether gain or loss:
for (i in 1:length(TGCA_cna_assoc_gr)) {
  
  temp_info <- TGCA_cna[TGCA_cna$Gene == TGCA_cna_assoc_gr$symbol[i],]

  if (nrow(temp_info) == 1) {
    TGCA_cna_assoc_gr$CNA[i] <- temp_info$CNA
    TGCA_cna_assoc_gr$CNA_freq[i] <- temp_info$Freq
    TGCA_cna_assoc_gr$CNA2[i] <- NA
    TGCA_cna_assoc_gr$CNA2_freq[i] <- NA
  } else {
    TGCA_cna_assoc_gr$CNA[i] <- temp_info$CNA[1]
    TGCA_cna_assoc_gr$CNA_freq[i] <- temp_info$Freq[1]
    TGCA_cna_assoc_gr$CNA2[i] <- temp_info$CNA[2]
    TGCA_cna_assoc_gr$CNA2_freq[i] <- temp_info$Freq[2]
  }

}

write.table(
  as.data.frame(TGCA_cna_assoc_gr),
  paste0(table_dir, "retrotransposon_breakpoint_assoc_DE_genes_TGCA_CNA_assoc.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

cancer_cna <- TGCA_cna[TGCA_cna$OncoKB_gene == "Yes",]
cancer_cna_assoc <- gene_symbols[gene_symbols %in% cancer_cna$Gene]
cancer_cna_assoc_gr <- TGCA_cna_assoc_gr[
  TGCA_cna_assoc_gr$symbol %in% cancer_cna_assoc
]
write.table(
  as.data.frame(cancer_cna_assoc_gr),
  paste0(table_dir, "retrotransposon_breakpoint_assoc_DE_cancer_genes_TGCA_CNA_assoc.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)

cancer_cna_assoc_df <- data.frame(
  symbol = cancer_cna_assoc_gr$symbol,
  sample = cancer_cna_assoc_gr$sample,
  logFC = cancer_cna_assoc_gr$logFC,
  FDR = cancer_cna_assoc_gr$FDR,
  CNA = cancer_cna_assoc_gr$CNA,
  CNA_freq = cancer_cna_assoc_gr$CNA_freq,
  CNA2 = cancer_cna_assoc_gr$CNA2,
  CNA2_freq = cancer_cna_assoc_gr$CNA2_freq
)

write.table(
  cancer_cna_assoc_df,
  paste0(table_dir, "retrotransposon_breakpoint_assoc_DE_cancer_genes_TGCA_CNA_assoc_simple.txt"),
  sep = "\t",
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE
)





