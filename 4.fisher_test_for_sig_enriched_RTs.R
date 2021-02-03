
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/") 

in_dir <- paste0(results_dir, "bp_assoc_RT/")
Robject_dir <- paste0(in_dir, "/Rdata/")
table_dir <- paste0(in_dir, "/tables/")

svaba_dir <- paste0(results_dir, "svaba/")
manta_dir <- paste0(results_dir, "manta/")

create_svaba_granges <- dget(paste0(func_dir, "create_svaba_granges.R"))
create_manta_granges <- dget(paste0(func_dir, "create_manta_granges.R"))

overlap_window <- 10000
all_breakpoint_source <- "TGCA"
incl_non_hg38 <- FALSE

if (all_breakpoint_source == "TGCA") {

  if (incl_non_hg38) {
    TGCA_dir <- paste0(ref_dir, "TGCA/OV_CNV/CNV_data/all_genome_builds/")
    out_dir <- paste0(in_dir, "vs_TGCA_breakpoints_all_genome_builds/")
  } else {
    TGCA_dir <- paste0(ref_dir, "TGCA/OV_CNV/CNV_data/hg38_only/")
    out_dir <- paste0(in_dir, "vs_TGCA_breakpoints_all_genome_builds/")
  }
  
} else {
  out_dir <- paste0(in_dir, "vs_low_conf_breakpoints/")
}


########################################################################
### 0. Load packages and functions ###
########################################################################

library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(dplyr)
library(tibble)

create_svaba_granges <- dget(paste0(func_dir, "create_svaba_granges.R"))
remove_overlaps <- dget(paste0(func_dir, "remove_overlaps.R"))


########################################################################
### 1. Load repetitive elements ###
########################################################################

# This is a list, with each memeber representing a GRanges object of 
# location:
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
repeats <- GRangesList(split(repeat_gtf, repeat_gtf$type))

# load all repeats within proximity to breakpoints:
bp_assoc_no_df <- read.table(
  paste0(table_dir, "retrotransposon_assoc_breakpoint_no.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = TRUE
)
repeats <- repeats[rownames(bp_assoc_no_df),]


########################################################################
### 2. Load breakpoints ###
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

if (all_breakpoint_source == "TGCA") {

  # load TGCA breakpoints:
  all_possible_breakpoints <- readRDS(
    paste0(TGCA_dir, "TGCA_OV_breakpoints.Rdata")
  )

  # add window either side:
  ranges(all_possible_breakpoints) <- IRanges(
    start = start(all_possible_breakpoints) - (overlap_window+50),
    end = end(all_possible_breakpoints) + (overlap_window+50)
  )

} else {

  # fetch low confidence breakpoints present in both svaba and manta
  # outputs:
  if (
    !file.exists(
      paste0(Robject_dir, "consensus_low_conf_breakpoints.Rdata")
    )
  ) {
  
    for (s in 1:length(sample_names)) {
      
      print(
        paste0(
          "Fetching and verifying breakpoints for ", sample_names[s]
        )
      )
    
      svaba_gr <- create_svaba_granges(
        samp_name = sample_names[s],
        in_dir = svaba_dir,
        SV_type = "unfiltered",
        bp_window = 100
      )
    
      manta_gr <- create_manta_granges(
        samp_name = sample_names[s],
        in_dir = manta_dir,
        SV_type = "unfiltered",
        bp_window = 100
      )
      
      # save bps that occur within 100 bp window from filtered svaba
      # and unfiltered manta:
      svaba_olaps <- findOverlaps(svaba_gr, manta_gr)
      bp <- svaba_gr[unique(queryHits(svaba_olaps))]
  
      # add filtered breakpoints:
      svaba_gr_filt <- create_svaba_granges(
        samp_name = sample_names[s],
        in_dir = svaba_dir,
        SV_type = "filtered",
        bp_window = 100
      )
    
      manta_gr_filt <- create_manta_granges(
        samp_name = sample_names[s],
        in_dir = manta_dir,
        SV_type = "filtered",
        bp_window = 100
      )
  
      bp <- c(bp, svaba_gr_filt, manta_gr_filt)
    
      # extand ranges to n kb each side of bp:
      ranges(bp) <- IRanges(
        start = start(bp) - overlap_window,
        end = end(bp) + overlap_window
      )
  
      print(
        paste0(
          "No. breakpoints for ", sample_names[s], " = ", length(bp)
        )
      )
  
      if (s==1) {
        low_conf_bp <- bp
      } else {
        low_conf_bp <- c(low_conf_bp, bp)
      }
  
      print(
        paste0(
          "Total No. breakpoints so far = ", length(low_conf_bp)
        )
      )
  
      writeLines("\n")
  
    }
  
    # remove duplicate overlap records:
    all_possible_breakpoints <- remove_overlaps(low_conf_bp)
  
    saveRDS(
      all_possible_breakpoints, 
      paste0(
        Robject_dir, 
        "non_overlapping_consensus_low_conf_breakpoints.Rdata"
      )
    )
  
  } else {
  
    all_possible_breakpoints <- readRDS(
      paste0(
        Robject_dir, 
        "non_overlapping_consensus_low_conf_breakpoints.Rdata"
      )
    )
  
  }

}

# Sanity check:
# Just to check that statistics are working as they should, I 
# would add a positive control
# Basically a subset of target break regions, or a mix of break 
# regions and random regions

# load genome breakpoints:
all_sample_bp <- readRDS(paste0(Robject_dir, "all_breakpoints.Rdata"))

# add all sample breakpoints to all possible breakpoints:
all_possible_bp <- c(
  all_possible_breakpoints , unlist(as(all_sample_bp, "GRangesList"))
)

# add 10 kb either side:
all_sample_bp_windows <- lapply(all_sample_bp, function(x) {
  start(x) <- start(x) - overlap_window
  end(x) <- end(x) + overlap_window
  return(x)
})

# find mean number of bp per genome:
mean(unlist(lapply(all_sample_bp_windows, length)))

# add total sample breakpoints and remove overlaps:
all_sample_bp_windows[["total_sample"]] <- unlist(
  as(
    all_sample_bp_windows, "GRangesList"
  )
)
all_sample_bp_windows[["total_sample"]] <- reduce(
  all_sample_bp_windows[["total_sample"]]
)

# retrieve chromosome lengths for random region control:
chr_lengths <- as.data.frame(SeqinfoForUCSCGenome("hg38"))
chr_lengths <- chr_lengths[
  grep("K|G|M|J", rownames(chr_lengths), invert = T),
]

# generate random positions (one set with length = mean no bp per genome, other with length = total no bp):
generateRandomPos <- function(n,chr,chr.sizes,width=1){
  random_chr <- sample(x=chr,size=n,prob=chr.sizes,replace=T)
  random_pos <- sapply(random_chr,function(chrTmp){sample(chr.sizes[chr==chrTmp],1)}) 
  res <- GRanges(random_chr,IRanges(random_pos,random_pos+width))
  return(res)
}

random_bp <- list()

set.seed(632)
random_bp$small_negative_control <- generateRandomPos(
  500, rownames(chr_lengths), chr_lengths$seqlengths
)
set.seed(37)
random_bp$large_negative_control <- generateRandomPos(
  8500, rownames(chr_lengths), chr_lengths$seqlengths
)

# expand random positions by 20099 bp:
random_bp <- lapply(random_bp, function(x) {
  start(x) <- start(x) - 10099
  end(x) <- end(x) + 10000
  return(x)
})

# add to all_sample_bp_windows:
all_sample_bp_windows <- append(all_sample_bp_windows, random_bp)


# create empty list for output:
output <- list()

for (i in 1:length(all_sample_bp_windows)) {

  print(
    paste0(
      "Testing for enrichment of retrotransposons at breakpoints within ", 
      names(all_sample_bp_windows)[i], "..."
    )
  )

  # add positive control coordinates to repeats list, consisting of 100
  # bp window ranges:
  repeats[["positive_control"]] <- all_sample_bp_windows[[i]][1:100]

  # set up output lists:
  overlaps <- list()
  output <- list()

  for (j in 1:length(names(repeats))) {

    print(
      paste0(
        "Testing for enrichment of ", names(repeats)[j], " at breakpoints..."
      )
    )

    # set up objects with all possible bp, and sample bp:
    breakpoints <- GRangesList()
    breakpoints[["sample"]] <- all_sample_bp_windows[[i]]
    breakpoints[["all_possible"]] <- all_possible_bp
  
    # fetch co-ordinates of retrotransposon:
    RT <- repeats[[names(repeats)[j]]]
  
    # detect overlaps between each set of breakpoints and RT coordinates:
    overlaps[[names(repeats)[j]]] <- sapply(breakpoints, function(x) {
      sum(
        countOverlaps(x, RT, ignore.strand=TRUE) > 0
      )
    })
  
    # no of overlaps of breakpoints in genome and repeat co-ordinates:
    breakR <- overlaps[[names(repeats)[j]]][["sample"]]
    # no of overlaps of repeat and all possible breakpoints - breakR:
    NbreakR <- overlaps[[names(repeats)[j]]][["all_possible"]] - breakR
    # no of breakpoints in genome - breakR:
    breakNR <- length(breakpoints[["sample"]]) - breakR
    # no of all possible breakpoints - no of breakpoints in genome - 
    # no of overlaps of repeat and all possible breakpoints
    NbreakNR <- length(breakpoints[["all_possible"]]) - 
      length(breakpoints[["sample"]]) -
      overlaps[[names(repeats)[j]]][["all_possible"]]
  
    output[[names(repeats)[j]]] <- fisher.test(
      matrix(c(breakR, breakNR, NbreakR, NbreakNR), nrow=2)
      , alternative="greater"
    )

  }

  if (length(overlaps) > 0) {

    RT_bp_overlaps <- do.call("rbind", overlaps)

    if (i==1) {

      all_outputs <- list(output)
      names(all_outputs)[i] <- names(all_sample_bp_windows)[i]

      all_pvals <- list(lapply(output, function(x) x$p.value))
      names(all_pvals)[i] <- names(all_sample_bp_windows)[i]
  
      all_overlaps <- list(RT_bp_overlaps)
      names(all_overlaps)[i] <- names(all_sample_bp_windows)[i]
  
    } else {
  
      all_outputs[[i]] <- output
      names(all_outputs)[i] <- names(all_sample_bp_windows)[i]

      all_pvals[[i]] <- lapply(output, function(x) x$p.value)
      names(all_pvals)[i] <- names(all_sample_bp_windows)[i]
  
      all_overlaps[[i]] <- RT_bp_overlaps
      names(all_overlaps)[i] <- names(all_sample_bp_windows)[i]
  
    }

  }

  print(paste0(i, " out of ", length(all_sample_bp_windows), " samples complete"))

}

all_pvals <- lapply(all_pvals, unlist)
sig_pvals <- lapply(all_pvals, function(x) x[x<0.1])

# format into df and save:
pval_df <- do.call("rbind", all_pvals)
pval_df <- subset(pval_df, select = -positive_control)
pval_df <- pval_df[
  (!rownames(pval_df) %in% c(
  "total_sample", "small_negative_control", "large_negative_control"
)),]
pval_df_rounded <- round(pval_df, 3)

FDR <- round(apply(pval_df_rounded, 2, function(x) p.adjust(x, "fdr")), 3)

FDR[FDR > 0.1] <- "nonsig"
pval_df_rounded[pval_df_rounded > 0.05] <- "nonsig"
pval_df[pval_df > 0.05] <- "nonsig"

write.table(
  pval_df_rounded,
  paste0(table_dir, "fisher_pvals.txt"),
  sep = "\t",
  quote = F,
  row.names = TRUE,
  col.names = TRUE
)

write.table(
  FDR,
  paste0(table_dir, "fisher_FDR.txt"),
  sep = "\t",
  quote = F,
  row.names = TRUE,
  col.names = TRUE
)

no_enriched_per_RT <- apply(
  pval_df_rounded, 
  2, 
  function(x) length(x) - sum(str_count(x, "nonsig"))
)

bp_assoc_RT_df <- merge(
  bp_assoc_no_df,
  data.frame(
    row.names = names(no_enriched_per_RT),
    no_enriched = no_enriched_per_RT
  ),
  by = 0
)

p_val_summaries <- apply(pval_df, 2, function(x) {

  if (all(x == "nonsig")) {
    return(NA)
  } else {

    p_only <- summary(as.numeric(x[x!="nonsig"]))
    minmax <- formatC(
      p_only[names(p_only) %in% c("Min.", "Max.")],
      format = "e",
      digits = 1
    )
    return(paste(sort(minmax), collapse = " - "))

  }
  
})

bp_assoc_RT_df <- merge(
  bp_assoc_RT_df,
  data.frame(
    Row.names = names(p_val_summaries),
    pvals = p_val_summaries
  ),
  by = "Row.names"
)

bp_assoc_RT_df <- arrange(
  bp_assoc_RT_df, 
  desc(total), 
  desc(per_genome), 
  desc(no_enriched)
)
colnames(bp_assoc_RT_df) <- c(
  "No. enriched breakpoints:",
  "Total",
  "Per genome",
  "No. samples enriched at breakpoints",
  "p-value range"
)

write.table(
  bp_assoc_RT_df,
  paste0(table_dir, "retrotransposon_assoc_breakpoint_and_enrichment_no.txt"),
  sep = "\t",
  quote = F,
  row.names = FALSE,
  col.names = TRUE
)



