
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/")

svaba_dir <- paste0(results_dir, "svaba/")
manta_dir <- paste0(results_dir, "manta/")
col_dir <- paste0(home_dir, "R/colour_palettes/")

RNA_name <- "hgsoc_repeats/RNA-seq-final"
RNA_dir <- paste0(home_dir, "projects/", RNA_name, "/")
DE_dir <- paste0(RNA_dir, "results/DE/without_primary_ascites/")

cancer_vs_ctl_dir <- paste0(
  DE_dir, "site/primary_vs_FT/tables/"
)
recurrent_vs_primary_dir <- paste0(
  DE_dir, "site/recurrent_vs_primary/tables/"
)
resistant_vs_sensitive_dir <- paste0(
  DE_dir, "drug_response/resistant_vs_sensitive/tables/"
)

out_dir <- paste0(results_dir, "bp_assoc_RT/")

Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))

overlap_window <- 10000
active_RT <- c(
  "L1PA2", "L1HS", "AluYk2", "AluYd8", "AluYb8", "AluYh3", "AluYh9"
)


########################################################################
### 0.  Load packages and functions ###
########################################################################

library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)

create_svaba_granges <- dget(paste0(func_dir, "create_svaba_granges.R"))
create_manta_granges <- dget(paste0(func_dir, "create_manta_granges.R"))
return_top_DE <- dget(paste0(func_dir, "return_top_DE.R"))


########################################################################
### 1.  Find consensus breakpoints in genomes ###
########################################################################

# fetch sample names:
sample_names <- grep(
  "sub", 
  list.files(svaba_dir, pattern = "AOCS"),
  invert = T,
  value = T
)
sample_names <- sample_names[
  sample_names %in% list.files(manta_dir, pattern = "AOCS")
]

if (!file.exists(paste0(Robject_dir, "all_breakpoints.Rdata"))) {

  for (s in 1:length(sample_names)) {
  
    print(
      paste0(
        "Fetching and verifying breakpoints for ", sample_names[s]
      )
    )
  
    svaba_gr_filt <- create_svaba_granges(
      samp_name = sample_names[s],
      in_dir = svaba_dir,
      SV_type = "filtered",
      bp_window = 100
    )
  
    svaba_gr <- create_svaba_granges(
      samp_name = sample_names[s],
      in_dir = svaba_dir,
      SV_type = "unfiltered",
      bp_window = 100
    )
  
    manta_gr_filt <- create_manta_granges(
      samp_name = sample_names[s],
      in_dir = manta_dir,
      SV_type = "filtered",
      bp_window = 100
    )
  
    manta_gr <- create_manta_granges(
      samp_name = sample_names[s],
      in_dir = manta_dir,
      SV_type = "unfiltered",
      bp_window = 100
    )
    
    # save bps that occur within 50 bp window from filtered svaba
    # and unfiltered manta:
    svaba_olaps <- findOverlaps(svaba_gr_filt, manta_gr)
    verified_svaba <- svaba_gr_filt[unique(queryHits(svaba_olaps))]
    
    # save bps that occur within 50 bp window from filtered manta
    # and unfiltered svaba:
    manta_olaps <- findOverlaps(manta_gr_filt, svaba_gr)
    verified_manta <- manta_gr_filt[unique(queryHits(manta_olaps))]

    combined <- c(verified_svaba, verified_manta)
    combined <- combined[!duplicated(combined)]
    
    # combine:
    if (s==1) {
      all_bps <- list(combined)
      names(all_bps)[s] <- sample_names[s] 
      print(
        paste0(
          "Total verified breakpoints = ", length(all_bps[[s]])
        )
      )
    } else {
      all_bps[[s]] <- combined
      names(all_bps)[s] <- sample_names[s] 
      print(
        paste0(
          "Total verified breakpoints = ", length(all_bps[[s]])
        )
      )
    }
  
    writeLines("\n")
  
  }
  
  saveRDS(all_bps, paste0(Robject_dir, "all_breakpoints.Rdata"))

} else {
  all_bps <- readRDS(paste0(Robject_dir, "all_breakpoints.Rdata"))
}


########################################################################
### 2. Load retrotransposon coordinates ###
########################################################################

# load repeat gtf:
if (!file.exists(paste0(Robject_dir, "repeat_gtf.Rdata"))) {

  repeat_gtf <- rtracklayer::import(
    paste0(genome_dir, "repeats.hg38.gtf")
  )
  
  # remove alternate chrs:
  repeat_gtf <- repeat_gtf[
    grep("K|G|M", seqnames(repeat_gtf), invert = T)
  ]

  saveRDS(repeat_gtf, paste0(Robject_dir, "repeat_gtf.Rdata"))

} else {
  repeat_gtf <- readRDS(paste0(Robject_dir, "repeat_gtf.Rdata"))
}


###### remove duplicate entries #####
#repeat_df <- as.data.frame(repeat_gtf)
#repeat_df <- subset(
#  repeat_df, 
#  select = -c(score, phase, feature, analysis)
#)
#repeat_df <- repeat_df %>% distinct()
#
#repeat_gtf <- GRanges(
#  seqnames = Rle(repeat_df$seqnames),
#  ranges = IRanges(
#    start = repeat_df$start, 
#    end = repeat_df$end
#  ),
#  strand = Rle("*"),
#  source = "repbase",
#  type = repeat_df$type,
#  class = repeat_df$class,
#  class_desc = repeat_df$class_desc
#)
#export(repeat_gtf, paste0(genome_dir, "repeats.hg38.gtf"))
#######

# load retrotransposons DE cancer vs ctls:
top_RT_symbols <- return_top_DE(
  in_dir = cancer_vs_ctl_dir,
  up_filename = "primary_vs_FT_upregulated_retrotransposon.txt",
  down_filename = "primary_vs_FT_downregulated_retrotransposon.txt"
)

# load retrotransposons DE primary vs recurrent cancer:
top_prim_vs_rec_RT_symbols <- return_top_DE(
  in_dir = recurrent_vs_primary_dir,
  up_filename = "recurrent_vs_primary_upregulated_retrotransposon.txt",
  down_filename = "recurrent_vs_primary_downregulated_retrotransposon.txt"
)

# load retrotransposons DE drug resistant vs sensitive cancer:
top_resistant_vs_sensitive_RT_symbols <- return_top_DE(
  in_dir = resistant_vs_sensitive_dir,
  up_filename = "resistant_vs_sensitive_GIN_upregulated_retrotransposon.txt",
  down_filename = "resistant_vs_sensitive_GIN_downregulated_retrotransposon.txt"
)

top_RT_symbols <- unique(
  c(
    top_RT_symbols, 
    top_prim_vs_rec_RT_symbols, 
    top_resistant_vs_sensitive_RT_symbols,
    active_RT
  )
)

# manually add that DE for known vs unknown GIN:
top_RT_symbols <- c(top_RT_symbols, "HERVK11D-int")

# take coordinates of top RTs:
top_RT <- repeat_gtf[repeat_gtf$type %in% top_RT_symbols]

# remove duplicate entries of L1HS/L1PA2:
top_RT <- top_RT[-which(duplicated(top_RT)),]


########################################################################
### 3. Find overlaps ###
########################################################################

# find overlaps and record details in breakpoint granges:
for (i in 1:length(all_bps)) {

  # number names of ranges to prevent bug:
  names(all_bps[[i]]) <- 1:length(all_bps[[i]])
  names(top_RT) <- 1:length(top_RT)

  olaps <- findOverlapsOfPeaks(all_bps[[i]], top_RT, maxgap = 10000)

  png(paste0(plot_dir, names(all_bps)[i], "_overlap_venn.png"))
    venn <- makeVennDiagram(
      olaps,
      fill=c("#009E73", "#F0E442"), # circle fill color
      col=c("#D55E00", "#0072B2"), #circle border color
      cat.col=c("#D55E00", "#0072B2")
    )
  dev.off()

  olap_peaks <- olaps$overlappingPeaks[["all_bps..i..///top_RT"]]

  sample_RT_bp <- GRanges(
    seqnames = Rle(olap_peaks$seqnames),
    ranges = IRanges(
      start = olap_peaks[,3], 
      end = olap_peaks[,4]
    ),
    strand = Rle("*"),
    bp = olap_peaks$bp,
    proximal_RT = olap_peaks$type,
    RT_start = olap_peaks[,11],
    RT_end = olap_peaks[,12],
    dist_from_bp = olap_peaks$shortestDistance,
    sample = olap_peaks$sample
  )

  if (i==1) {
    all_RT_bp <- list(sample_RT_bp)
  } else {
    all_RT_bp[[i]] <- sample_RT_bp
  }

}

# combine as one:
RT_bp <- unlist(as(all_RT_bp, "GRangesList"))

# split according to RT:
RT_bp_spl <- split(RT_bp, as.character(RT_bp$proximal_RT))

# check number of bp associations for each:
bp_assoc_no <- unlist(lapply(RT_bp_spl, length))
bp_assoc_no_per_genome <- round(bp_assoc_no/length(all_RT_bp))

# save retrotransposon assoc breakpoints:
saveRDS(
  RT_bp, 
  paste0(Robject_dir, "retrotransposon_assoc_breakpoints.Rdata")
)
export(RT_bp, paste0(table_dir, "retrotransposon_assoc_breakpoints.gtf"))

# save A0CS_063 e.g.:
AOCS_063_RT_bp <- RT_bp[
  RT_bp$sample == "AOCS-063" & RT_bp$proximal_RT %in% c(
    "L1MD3", "L1HS", "L1PA2", "AluYb8", "AluYh3"
  )
]
# add 10000 either side of ranges:
AOCS_063_RT_bp_windows <- AOCS_063_RT_bp
ranges(AOCS_063_RT_bp_windows) <- IRanges(
  start = start(AOCS_063_RT_bp)-10000,
  end = end(AOCS_063_RT_bp)+10000
)

export(
  AOCS_063_RT_bp_windows, 
  paste0(table_dir, "AOCS_063_retrotransposon_assoc_breakpoint_windows.gtf")
)

# swap ranges:
bp_RT <- RT_bp
ranges(bp_RT) <- IRanges(start = RT_bp$RT_start, end = RT_bp$RT_end)
bp_RT$bp_start <- start(RT_bp)
bp_RT$bp_end <- end(RT_bp)
bp_RT$RT_start <- NULL
bp_RT$RT_end <- NULL

# save breakpoint assoc retrotransposons:
export(bp_RT, paste0(table_dir, "breakpoint_assoc_retrotransposons.gtf"))

# save A0CS_063 e.g.:
AOCS_063_bp_RT <- bp_RT[
  bp_RT$sample == "AOCS-063" & bp_RT$proximal_RT %in% c(
    "L1MD3", "L1HS", "L1PA2", "AluYb8", "AluYh3"
  )
]
export(
  AOCS_063_bp_RT, 
  paste0(table_dir, "AOCS_063_breakpoint_assoc_retrotransposons.gtf")
)

# save all AOCS_063 breakpoints:
export(all_bps[[1]], paste0(table_dir, "AOCS_063_breakpoints.gtf"))

# save assoc nos:
bp_assoc_no_df <- data.frame(
  total = bp_assoc_no,
  per_genome = bp_assoc_no_per_genome
)

bp_assoc_no_df <- dplyr::arrange(
  bp_assoc_no_df, 
  desc(total), 
  desc(per_genome), 
)

write.table(
  bp_assoc_no_df,
  paste0(table_dir, "retrotransposon_assoc_breakpoint_no.txt"),
  quote = F,
  sep = "\t",
  row.names = T,
  col.names = T
)



