
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")

in_dir <- paste0(results_dir, "bp_assoc_RT/tables/")

out_path <- paste0(results_dir, "bp_assoc_RT/gviz/")
Robject_dir <- paste0(out_path, "Rdata/")
plot_dir <- paste0(out_path, "plots/")


########################################################################
### 0.  Load packages ###
########################################################################

library(Gviz)
library(rtracklayer)
library(GenomicRanges)


########################################################################
### 1. Load data and genes to snapshot ###
########################################################################

RT <- import(
  paste0(
  	in_dir, 
  	"AOCS-114/AOCS-114_breakpoint_assoc_retrotransposons.sorted.gtf"
  )
)
bp <- import(
  paste0(
  	in_dir, 
  	"AOCS-114/AOCS-114_retrotransposon_assoc_breakpoints.sorted.gtf"
  )
)
bp_windows <- import(
  paste0(
  	in_dir, 
  	"AOCS-114/AOCS-114_retrotransposon_assoc_breakpoint_windows.sorted.gtf"
  )
)

ss_genes <- c("UBE2C")