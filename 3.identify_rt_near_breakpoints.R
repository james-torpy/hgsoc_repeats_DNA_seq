
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
cancer_vs_ctl_dir <- paste0(
  RNA_dir, "/results/DE/site/primary_vs_FT/tables/"
)
recurrent_vs_primary_dir <- paste0(
  RNA_dir, "/results/DE/site/recurrent_vs_primary/tables/"
)

out_dir <- paste0(results_dir, "bp_assoc_RT/")

Robject_dir <- paste0(out_dir, "Rdata/")
system(paste0("mkdir -p ", Robject_dir))
table_dir <- paste0(out_dir, "tables/")
system(paste0("mkdir -p ", table_dir))
plot_dir <- paste0(out_dir, "plots/")
system(paste0("mkdir -p ", plot_dir))

overlap_window <- 10000
active_RT <- c("L1PA2", "L1HS", "AluYk2", "AluYd8", "AluYb8", "AluYh3")


########################################################################
### 0.  Load packages and functions ###
########################################################################

library(GenomicRanges)

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
  
    # extand ranges to n kb each side of bp:
    ranges(verified_svaba) <- IRanges(
      start = verified_svaba$bp - overlap_window,
      end = verified_svaba$bp + overlap_window
    )
    
    # save bps that occur within 50 bp window from filtered manta
    # and unfiltered svaba:
    manta_olaps <- findOverlaps(manta_gr_filt, svaba_gr)
    verified_manta <- manta_gr_filt[unique(queryHits(manta_olaps))]
  
    ranges(verified_manta) <- IRanges(
      start = verified_manta$bp - overlap_window,
      end = verified_manta$bp + overlap_window
    )
    
    # combine:
    if (s==1) {
      all_bps <- list(c(verified_svaba, verified_manta))
      names(all_bps)[s] <- sample_names[s] 
      print(
        paste0(
          "Total verified breakpoints = ", length(all_bps[[s]])
        )
      )
    } else {
      all_bps[[s]] <- c(verified_svaba, verified_manta)
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
### 2. Find retrotransposon overlaps ###
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

top_both_symbols <- top_prim_vs_rec_RT_symbols[
  top_prim_vs_rec_RT_symbols %in% top_RT_symbols
]

top_RT_symbols <- unique(
  c(
    top_RT_symbols, 
    top_both_symbols, 
    active_RT
  )
)

# take coordinates of top RTs:
top_RT <- repeat_gtf[repeat_gtf$type %in% top_RT_symbols]

# find overlaps and record details in breakpoint granges:
all_bp_RT <- lapply(all_bps, function(x) {

  x$proximal_RT <- NA
  x$RT_start <- NA
  x$RT_end <- NA

  # find overlaps between DE RTs and bp windows:
  olaps <- suppressWarnings(findOverlaps(x, top_RT))

  # add info of overlapping RTs:
  x$proximal_RT[queryHits(olaps)] <- as.character(
    top_RT$type[subjectHits(olaps)]
  )
  x$RT_start[queryHits(olaps)] <- start(top_RT)[subjectHits(olaps)]
  x$RT_end[queryHits(olaps)] <- end(top_RT)[subjectHits(olaps)]

  # add distances of each RT from bps:
  x$start_dist_from_bp <- abs(x$RT_start-x$bp)
  x$end_dist_from_bp <- abs(x$RT_end-x$bp)
  x$dist_from_bp <- pmin(x$start_dist_from_bp, x$end_dist_from_bp)
  x$start_dist_from_bp <- NULL
  x$end_dist_from_bp <- NULL

  # remove bps without proximal RTs:
  x <- x[!is.na(x$proximal_RT)]
  
  return(x)

})

# add sample ids and combine all:
for (i in 1:length(all_bp_RT)) {
  all_bp_RT[[i]]$sample <- names(all_bp_RT)[i]
}
bp_RT <- unlist(as(all_bp_RT, "GRangesList"))

# split according to RT:
bp_RT_spl <- split(bp_RT, bp_RT$proximal_RT)

# check number of bp associations for each:
bp_assoc_no <- unlist(lapply(bp_RT_spl, length))
bp_assoc_no_per_genome <- round(bp_assoc_no/length(all_bp_RT))

# save data:
saveRDS(
  bp_RT, 
  paste0(Robject_dir, "breakpoint_assoc_retrotransposons.Rdata")
)
bp_RT_df <- data.frame(
  sample = bp_RT$sample,
  chr = seqnames(bp_RT),
  bp_pos = bp_RT$bp,
  proximal_RT = bp_RT$proximal_RT,
  RT_start = bp_RT$RT_start,
  RT_end = bp_RT$RT_end,
  dist_from_bp = bp_RT$dist_from_bp
)
write.table(
  bp_RT_df,
  paste0(table_dir, "breakpoint_assoc_retrotransposons.txt"),
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

bp_assoc_no_df <- data.frame(
  total = bp_assoc_no,
  per_genome = bp_assoc_no_per_genome
)

bp_assoc_no_df <- bp_assoc_no_df[
  order(bp_assoc_no_df$per_genome, decreasing = T),
]

write.table(
  bp_assoc_no_df,
  paste0(table_dir, "retrotransposon_assoc_breakpoint_no.txt"),
  quote = F,
  sep = "\t",
  row.names = T,
  col.names = T
)



