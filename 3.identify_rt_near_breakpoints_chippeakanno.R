
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

sig_cutoff <- 0.05
logfc_cutoff <- 0.7
overlap_window <- 10000
active_RT <- c(
  "L1PA2", "L1HS", "AluYk2", "AluYd8", "AluYb8", "AluYh3", "AluYh9"
)
gtf_samples <- c(
  "AOCS-064", "AOCS-075", "AOCS-076", "AOCS-080", "AOCS-083",
  "AOCS-085", "AOCS-090", "AOCS-094", "AOCS-107", "AOCS-112", 
  "AOCS-114", "AOCS-116", "AOCS-122", "AOCS-128", "AOCS-130", 
  "AOCS-133", "AOCS-137"
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
    
    # save all low confidence breakpoints to retain mates:
    low_conf_combined <- c(svaba_gr, manta_gr)
    low_conf_combined <- low_conf_combined[!duplicated(low_conf_combined)]

    # combine:
    if (s==1) {

      all_bps <- list(combined)
      names(all_bps)[s] <- sample_names[s] 
      print(
        paste0(
          "Total verified breakpoints = ", length(all_bps[[s]])
        )
      )

      low_conf_bps <- list(low_conf_combined)
      names(low_conf_bps)[s] <- sample_names[s] 
      print(
        paste0(
          "Total non-verified breakpoints = ", length(low_conf_bps[[s]])
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

      low_conf_bps[[s]] <- low_conf_combined
      names(low_conf_bps)[s] <- sample_names[s] 
      print(
        paste0(
          "Total non-verified breakpoints = ", length(low_conf_bps[[s]])
        )
      )

    }
  
    writeLines("\n")
  
  }
  
  saveRDS(all_bps, paste0(Robject_dir, "all_breakpoints.Rdata"))
  saveRDS(
    low_conf_bps, 
    paste0(Robject_dir, "all_low_confidence_breakpoints.Rdata")
  )

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

# load retrotransposons DE cancer vs ctls:
top_RT_symbols <- return_top_DE(
  in_dir = cancer_vs_ctl_dir,
  up_filename = "primary_vs_FT_upregulated_retrotransposon.txt",
  down_filename = "primary_vs_FT_downregulated_retrotransposon.txt",
  max_pval = sig_cutoff,
  min_logfc = logfc_cutoff
)

top_RT_symbols <- unique(
  c(
    top_RT_symbols,
    active_RT
  )
)

# take coordinates of top RTs:
top_RT <- repeat_gtf[repeat_gtf$type %in% top_RT_symbols]

# remove duplicate entries of L1HS/L1PA2:
top_RT <- top_RT[-which(duplicated(top_RT)),]


########################################################################
### 3. Find overlaps ###
########################################################################

if (
  !file.exists(
    paste0(Robject_dir, "retrotransposon_breakpoint_assoc_genes.Rdata")
  ) & !file.exists(
    paste0(Robject_dir, "retrotransposon_assoc_breakpoint_no.Rdata")
  ) & !file.exists(
    paste0(Robject_dir, "retrotransposon_assoc_breakpoints.Rdata")
  )
) {

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
  
    # determine column numbers of different start and end co-ordinates:
    both_start_ind <- grep("start", colnames(olap_peaks))
    both_end_ind <- grep("end", colnames(olap_peaks))
    start_ind <- both_start_ind[1]
    end_ind <- both_end_ind[1]
    RT_start_ind <- both_start_ind[2]
    RT_end_ind <- both_end_ind[2]
  
    sample_RT_bp <- GRanges(
      seqnames = Rle(olap_peaks$seqnames),
      ranges = IRanges(
        start = olap_peaks[, start_ind], 
        end = olap_peaks[, end_ind]
      ),
      strand = Rle("*"),
      bp = olap_peaks$bp,
      CNV_type = olap_peaks$CNV_type,
      proximal_RT = olap_peaks$type,
      RT_start = olap_peaks[, RT_start_ind],
      RT_end = olap_peaks[, RT_end_ind],
      dist_from_bp = olap_peaks$shortestDistance,
      sample = olap_peaks$sample
    )
  
    if (i==1) {
      all_RT_bp <- list(sample_RT_bp)
    } else {
      all_RT_bp[[i]] <- sample_RT_bp
    }
  
    # find overlaps of RT-assoc bps and genes:
    library(EnsDb.Hsapiens.v86)
    gene_annot <- toGRanges(EnsDb.Hsapiens.v86, feature="gene")
  
    bp_olaps <- olaps$peaklist[["all_bps..i..///top_RT"]]
    gene_overlaps <- annotatePeakInBatch(bp_olaps, AnnotationData=gene_annot, 
      output="overlapping", maxgap=5000L)
    gene_overlaps$gene_name <- gene_annot$gene_name[
      match(gene_overlaps$feature, names(gene_annot))
    ]
  
    # remove NAs, unwanted rows and columns:
    gene_overlaps <- gene_overlaps[grep("NA", names(gene_overlaps), invert=T)]
    gene_overlaps <- gene_overlaps[
      grep("RP[0-9][0-9]-|CTD-[0-9]|AC[0-9].*\\.", gene_overlaps$gene_name, invert=T)
    ]
    mcols(gene_overlaps) <- subset(
      mcols(gene_overlaps), 
      select = c(peakNames, start_position, end_position, insideFeature, distancetoFeature,
        shortestDistance, fromOverlappingOrNearest, gene_name)
    )
    colnames(mcols(gene_overlaps)) <- c(
      "peaks", "start_pos", "end_pos", "inside", 
      "distance", "shortest_dist", "olap_or_nearest", "symbol"
    )
    gene_overlaps$peaks <- gsub(
      "__", "", 
      gsub("..i..", "",  gene_overlaps$peaks)
    )
    gene_overlaps$sample <- names(all_bps)[i]
  
    if (i==1) {
      RT_bp_genes <- list(gene_overlaps)
    } else {
      RT_bp_genes[[i]] <- gene_overlaps
    }
  
  }
  
  # save retrotransposon-associated breakpoint-associated genes:
  names(RT_bp_genes) <- names(all_bps)
  saveRDS(
    RT_bp_genes, 
    paste0(Robject_dir, "retrotransposon_breakpoint_assoc_genes.Rdata")
  )
  
  # combine retrotransposon-associated breakpoints as one:
  RT_bp <- unlist(as(all_RT_bp, "GRangesList"))
  
  # split according to RT:
  RT_bp_spl <- split(RT_bp, as.character(RT_bp$proximal_RT))
  
  # check number of bp associations for each and save:
  bp_assoc_no <- list(
    total = unlist(lapply(RT_bp_spl, length))
  )
  bp_assoc_no$per_genome <- round(bp_assoc_no$total/length(all_RT_bp))
  saveRDS(
    bp_assoc_no, 
    paste0(Robject_dir, "retrotransposon_assoc_breakpoint_no.Rdata")
  )
  
  # save retrotransposon assoc breakpoints:
  ranges(RT_bp) <- IRanges(start = RT_bp$bp, end = RT_bp$bp)
  mcols(RT_bp) <- subset(mcols(RT_bp), select = -bp)
  saveRDS(
    RT_bp, 
    paste0(Robject_dir, "retrotransposon_assoc_breakpoints.Rdata")
  )
  export(RT_bp, paste0(table_dir, "retrotransposon_assoc_breakpoints.gtf"))

} else {

  RT_bp_genes <- readRDS(
    paste0(Robject_dir, "retrotransposon_breakpoint_assoc_genes.Rdata")
  )
  bp_assoc_no <- readRDS(
    paste0(Robject_dir, "retrotransposon_assoc_breakpoint_no.Rdata")
  )
  RT_bp <- readRDS(
    paste0(Robject_dir, "retrotransposon_assoc_breakpoints.Rdata")
  )

}

# swap ranges:
bp_RT <- RT_bp
ranges(bp_RT) <- IRanges(start = RT_bp$RT_start, end = RT_bp$RT_end)
bp_RT$bp <- start(RT_bp)
bp_RT$RT_start <- NULL
bp_RT$RT_end <- NULL

# save breakpoint assoc retrotransposons:
export(bp_RT, paste0(table_dir, "breakpoint_assoc_retrotransposons.gtf"))

# save assoc nos:
bp_assoc_no_df <- data.frame(
  total = bp_assoc_no$total,
  per_genome = bp_assoc_no$per_genome
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


########################################################################
### 4. Save, gtfs for individual samples ###
########################################################################

for (i in 1:length(gtf_samples)) {

  # create sample directory:
    sample_dir <- paste0(table_dir, gtf_samples[i], "/")
    system(paste0("mkdir -p ", sample_dir))

  if (
    !file.exists(
      paste0(
        sample_dir, gtf_samples[i], "_retrotransposon_assoc_breakpoints.gtf"
      )
    )
  ) {

    # save sample e.g.:
    sample_RT_bp <- RT_bp[RT_bp$sample == gtf_samples[i]]
    # add 10000 either side of ranges:
    sample_RT_bp_windows <- sample_RT_bp
    ranges(sample_RT_bp_windows) <- IRanges(
      start = start(sample_RT_bp)-10000,
      end = end(sample_RT_bp)+10000
    )
    
    export(
      sample_RT_bp, 
      paste0(sample_dir, gtf_samples[i], "_retrotransposon_assoc_breakpoints.gtf")
    )
  
    export(
      sample_RT_bp_windows, 
      paste0(sample_dir, gtf_samples[i], "_retrotransposon_assoc_breakpoint_windows.gtf")
    )
    
    # save A0CS_063 e.g.:
    sample_bp_RT <- bp_RT[
      bp_RT$sample == gtf_samples[i]
    ]
    export(
      sample_bp_RT, 
      paste0(sample_dir, gtf_samples[i], "_breakpoint_assoc_retrotransposons.gtf")
    )
    
    # save all sample breakpoints:
    export(all_bps[[1]], paste0(sample_dir, gtf_samples[i], "_breakpoints.gtf"))
    
    # save sample RT bp-assoc gene windows:
    export(
      RT_bp_genes[[1]], 
      paste0(sample_dir, gtf_samples[i], "_retrotransposon_breakpoint_assoc_genes.gtf")
    )

  }

}


