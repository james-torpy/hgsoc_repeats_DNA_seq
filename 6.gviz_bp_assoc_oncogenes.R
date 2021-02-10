
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/jamestorpy/clusterHome/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/")
Robject_dir <- paste0(results_dir, "bp_assoc_RT/Rdata/")

in_dir <- paste0(results_dir, "bp_assoc_RT/tables/")
bam_dir <- paste0(results_dir, "bwa/")

out_path <- paste0(results_dir, "bp_assoc_RT/gviz/")
plot_dir <- paste0(out_path, "plots/")

system(paste0("mkdir -p ", plot_dir))


########################################################################
### 0.  Load packages ###
########################################################################

library(Gviz)
library(rtracklayer)
library(GenomicRanges)


########################################################################
### 1. Load non-sample specific data and tracks ###
########################################################################

# load all bps to look up joins:
all_bp <- readRDS(paste0(Robject_dir, "all_breakpoints.Rdata"))

# load pre-generated gene models of specific genes:
goi <- c(
  "EPHA3", "ANKRD11", "MECOM", "EXT2",
  "ATM", "POLD1", "CXCR2", "UBE2C", 
  "DKK2", "FANCD2", "BARD1", "CCDC6",
  "FBXW7", "C11orf65"
)

gene_annot <- import(
  paste0(genome_dir, "custom.genes.gencode.v35.gtf")
)
gene_annot <- gene_annot[
  gene_annot$gene_name %in% goi
  ]

# format for Gviz:
# remove transcript, gene and CDS entries, as well as non-essential columns:
gene_annot <- gene_annot[!(gene_annot$type %in% c("gene", "transcript", "CDS"))]
mcols(gene_annot) <- subset(
  mcols(gene_annot),
  select = c(gene_type, gene_id, exon_id, transcript_id, type, gene_name)
)

# convert into df:
gene_df <- data.frame(
  chromosome = seqnames(gene_annot),
  start = start(gene_annot),
  end = end(gene_annot),
  width = width(gene_annot),
  strand = strand(gene_annot),
  feature = gene_annot$type,
  gene = gene_annot$gene_id,
  exon = gene_annot$exon_id,
  transcript = gene_annot$transcript_id,
  symbol = gene_annot$gene_name
)

# create genomic coordinate track:
gtrack <- GenomeAxisTrack()


########################################################################
### 2. Prepare loading of sample data ###
########################################################################

# setup snapshot list:
snapshot_list <- list(
  data.frame(
    sample = "AOCS-080",
    gene_symbol = "ATM",
    chr = "chr11",
    start_coord = as.integer(108200000),
    end_coord = as.integer(108445000),
    stringsAsFactors = F,
    cov_int = 20000,
    addstring="none",
    addgene = "C11orf65",
    addbp = 108346154
  ),
  data.frame(
    sample = "AOCS-080",
    gene_symbol = "ATM",
    chr = "chr11",
    start_coord = as.integer(107750000),
    end_coord = as.integer(109000000),
    stringsAsFactors = F,     
    cov_int = 20000,
    addstring = "zoomed_out",
    addgene = "C11orf65",
    addbp = 108346154
  ),
  data.frame(
    sample = "AOCS-083",
    gene_symbol = "POLD1",
    chr = "chr19",
    start_coord = as.integer(50375000),
    end_coord = as.integer(50420000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none"
  ),
  data.frame(
    sample = "AOCS-083",
    gene_symbol = "POLD1",
    chr = "chr19",
    start_coord = as.integer(13680000),
    end_coord = as.integer(50550000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="zoomed_out"
  ),
  data.frame(
    sample = "AOCS-128",
    gene_symbol = "FANCD2",
    chr = "chr3",
    start_coord = as.integer(10000000),
    end_coord = as.integer(10110000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none",
    addbp = 10029933
  ),
  data.frame(
    sample = "AOCS-128",
    gene_symbol = "FANCD2",
    chr = "chr3",
    start_coord = as.integer(9900000),
    end_coord = as.integer(10200000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="zoomed_out",
    addbp = 10029933
  ),
  data.frame(
    sample = "AOCS-130",
    gene_symbol = "DKK2",
    chr = "chr4",
    start_coord = as.integer(106890000),
    end_coord = as.integer(107045000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none"
  ),
  data.frame(
    sample = "AOCS-130",
    gene_symbol = "DKK2",
    chr = "chr4",
    start_coord = as.integer(106800000),
    end_coord = as.integer(107350000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="zoomed_out"
  ),
  data.frame(
    sample = "AOCS-137",
    gene_symbol = "FBXW7",
    chr = "chr4",
    start_coord = as.integer(152250000),
    end_coord = as.integer(152600000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none"
  ),
  data.frame(
    sample = "AOCS-137",
    gene_symbol = "FBXW7",
    chr = "chr4",
    start_coord = as.integer(152000000),
    end_coord = as.integer(153000000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="zoomed_out"
  ),
  data.frame(
    sample = "AOCS-085",
    gene_symbol = "CXCR2",
    chr = "chr2",
    start_coord = as.integer(218120000),
    end_coord = as.integer(218160000),
    stringsAsFactors = F,     
    cov_int = 20000,
    addstring = "none"
  ),
  data.frame(
    sample = "AOCS-085",
    gene_symbol = "CXCR2",
    chr = "chr2",
    start_coord = as.integer(217000000),
    end_coord = as.integer(218750000),
    stringsAsFactors = F,     
    cov_int = 20000,
    addstring="zoomed_out",
    addbp = 217395522
  ),
  data.frame(
    sample = "AOCS-064",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = as.integer(89000000),
    end_coord = as.integer(89600000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none",
    addbp=89173179
  ),
  data.frame(
    sample = "AOCS-112",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = as.integer(89050000),
    end_coord = as.integer(89600000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none",
    start_coord = as.integer(89000000),
    end_coord = as.integer(89600000)
  ),
  data.frame(
    sample = "AOCS-122",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = as.integer(89000000),
    end_coord = as.integer(89600000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring ="none"
  ),
  data.frame(
    sample = "AOCS-130",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = as.integer(89000000),
    end_coord = as.integer(89600000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none"
  ),
  data.frame(
    sample = "AOCS-075",
    gene_symbol = "MECOM",
    chr = "chr3",
    start_coord = as.integer(168800000),
    end_coord = as.integer(169700000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="none"
  ),
  data.frame(
    sample = "AOCS-075",
    gene_symbol = "MECOM",
    chr = "chr3",
    start_coord = as.integer(168000000),
    end_coord = as.integer(170700000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="zoomed_out"
  ),
  data.frame(
    sample = "AOCS-114",
    gene_symbol = "UBE2C",
    chr = "chr20",
    start_coord = as.integer(45809000),
    end_coord = as.integer(45818500),
    stringsAsFactors = F,     
    cov_int = 500, 
    addstring="none"
  ),
  data.frame(
    sample = "AOCS-114",
    gene_symbol = "UBE2C",
    chr = "chr20",
    start_coord = as.integer(44500000),
    end_coord = as.integer(46500000),
    stringsAsFactors = F,     
    cov_int = 20000, 
    addstring="zoomed_out"
  )
)

for (i in 1:length(snapshot_list)) {
  
  print(
    paste0(
      "Plotting ", snapshot_list[[i]]$gene_symbol, 
      " snapshot within ", snapshot_list[[i]]$sample
    )
  )
  
  # create chromosome ideogram track:
  itrack <- IdeogramTrack(
    genome = "hg38", 
    chromosome = snapshot_list[[i]]$chr
  )
  
  # load breakpoint and retrotransposon gtfs:
  gtfs <- list()
  
  gtfs$RT <- import(
    paste0(
      in_dir, snapshot_list[[i]]$sample, "/", snapshot_list[[i]]$sample, 
      "_breakpoint_assoc_retrotransposons.sorted.gtf"
    )
  )
  gtfs$bp <- import(
    paste0(
      in_dir, snapshot_list[[i]]$sample, "/", snapshot_list[[i]]$sample, 
      "_retrotransposon_assoc_breakpoints.sorted.gtf"
    )
  )
  
  
  # remove duplicates:
  gtfs <- lapply(gtfs, function(x) x[!duplicated(x)])
  
  if ("addbp" %in% colnames(snapshot_list[[i]])) {
    # add additional breakpoint:
    gtfs$bp <- c(
      gtfs$bp,
      GRanges(
        seqnames = snapshot_list[[i]]$chr,
        ranges = IRanges(
          start = snapshot_list[[i]]$addbp,
          end = snapshot_list[[i]]$addbp
        ),
        strand = "*"
      )
    )
  }
  
  if ("addbp2" %in% colnames(snapshot_list[[i]])) {
    # add additional breakpoint:
    gtfs$bp <- c(
      gtfs$bp,
      GRanges(
        seqnames = snapshot_list[[i]]$chr,
        ranges = IRanges(
          start = snapshot_list[[i]]$addbp2,
          end = snapshot_list[[i]]$addbp2
        ),
        strand = "*"
      )
    )
  }
  
  
  ########################################################################
  ### 3. Load additional tracks ###
  ########################################################################
  
  # create breakpoint track:
  bp_track = AnnotationTrack(
    gtfs$bp, 
    name = "Bp",
    col = "#4A70D1",
    fill = "#4A70D1",
    background.title = "#4A70D1"
  )
  
  # set genome and chromosome for bp_track:
  genome(bp_track) <- "hg38"
  chromosome(bp_track) <- snapshot_list[[i]]$chr
  
  # format retrotransposon gene data into df:
  RT_df <- data.frame(
    chromosome = seqnames(gtfs$RT),
    start = start(gtfs$RT),
    end = end(gtfs$RT),
    width = width(gtfs$RT),
    strand = strand(gtfs$RT),
    feature = "exon",
    symbol = gtfs$RT$proximal_RT,
    stringsAsFactors = F
  )
  
  # create retrotranposon gene model:
  RT_track <- GeneRegionTrack(
    RT_df, 
    genome = "hg38", 
    chromosome = snapshot_list[[i]]$chr,
    name = "Rt",
    transcriptAnnotation = "symbol",
    col = "#1B9E77",
    fill = "#1B9E77",
    background.title = "#1B9E77"
  )
  
  # create protein coding gene model:
  grtrack <- GeneRegionTrack(
    gene_df, 
    genome = "hg38", 
    chromosome = snapshot_list[[i]]$chr,
    name = "Breakpoint assoc. gene",
    transcriptAnnotation = "symbol",
    background.title = "#F9D480"
  )
  
  
  ########################################################################
  ### 4. Load sample coverage and isolate required segment ###
  ########################################################################
  
  # fetch mean coverage:
  mean_cov <- round(
    read.table(
      paste0(bam_dir, snapshot_list[[i]]$sample, "-1.mean.coverage.txt"),
      header = F
    )[1,1]
  )
  
  if (snapshot_list[[i]]$addstring == "none") {
    bedfile <- paste0(
      bam_dir, snapshot_list[[i]]$sample,
      "-1.", snapshot_list[[i]]$gene_symbol,
      ".coverage.bed"
    )
  } else {
    bedfile <- paste0(
      bam_dir, snapshot_list[[i]]$sample,
      "-1.", snapshot_list[[i]]$gene_symbol,
      "_", snapshot_list[[i]]$addstring,
      ".coverage.bed"
    )
  }
  
  if (!file.exists(bedfile)) {
    
    # create temp bed file containing region of interest:
    write.table(
      subset(snapshot_list[[i]], select = -c(sample, gene_symbol)),
      paste0(bam_dir, snapshot_list[[i]]$sample, "-1.temp.bed"),
      sep = "\t",
      col.names = F,
      row.names = F,
      quote = F
    )
    
    # create bed file with specific coordinates:
    print(paste0("Creating segment coverage bed file..."))
    system(
      paste0(
        "bedtools intersect -a ", bam_dir, snapshot_list[[i]]$sample,
        "-1.coverage.bed -b ", bam_dir, snapshot_list[[i]]$sample,
        "-1.temp.bed > ", bedfile
      )
    )
    
    # remove temp bedfile:
    system(
      paste0(
        "rm ", bam_dir, snapshot_list[[i]]$sample, "-1.temp.bed")
    )
    
  }
  
  # load bed file:
  cov <- import(bedfile)
  colnames(mcols(cov)) <- "coverage"
  
  # centre values on mean coverage:
  cov$centered_coverage <- as.numeric(cov$coverage) - mean_cov
  mcols(cov) <- subset(mcols(cov), select = centered_coverage)
  
  # bin coverage values and calculate means:
  # create bins:
  bin_vec <- snapshot_list[[i]]$start_coord:snapshot_list[[i]]$end_coord
  bin_list <- split(bin_vec, ceiling(seq_along(bin_vec)/snapshot_list[[i]]$cov_int))
  bin_coords <- do.call("rbind", lapply(bin_list, function(x) c(min(x), max(x))))
  bin_gr <- GRanges(
    seqnames <- Rle(snapshot_list[[i]]$chr),
    ranges = IRanges(start = bin_coords[,1], end = bin_coords[,2]),
    strand = "*"
  )
  
  # create Rle object from coverage values:
  cov_rle <- Rle(values = cov$centered_coverage, lengths = width(cov))
  cov_rle <- mcolAsRleList(cov, "centered_coverage")
  binned_cov <- binnedAverage(bin_gr, cov_rle,  "binned_centered_coverage")
  binned_cov$binned_centered_coverage[1] <- mean(
    cov_rle[[1]]@values[
      !is.na(cov_rle[[1]]@values)][1:width(binned_cov)[1]
                                   ]
  )
  binned_cov$binned_centered_coverage <- round(
    binned_cov$binned_centered_coverage,
    1
  )
  
  # add extra range to start and end for plotting:
  if (width(binned_cov)[1] > 2) {
    
    binned_cov <- c(
      GRanges(
        seqnames = Rle(snapshot_list[[i]]$chr),
        ranges = IRanges(
          start = start(binned_cov)[1],
          end = start(binned_cov)[1] + 1
        ),
        strand = "*",
        binned_centered_coverage = binned_cov$binned_centered_coverage[1]
      ),
      binned_cov
    )
    start(binned_cov)[2] <- start(binned_cov)[2]+2
    
  }
  
  if (width(binned_cov)[length(binned_cov)] > 2) {
    
    binned_cov <- c(
      binned_cov,
      GRanges(
        seqnames = Rle(snapshot_list[[i]]$chr),
        ranges = IRanges(
          start = end(binned_cov)[length(binned_cov)] - 1,
          end = end(binned_cov)[length(binned_cov)]
        ),
        strand = "*",
        binned_centered_coverage = binned_cov$binned_centered_coverage[
          length(binned_cov)
          ]
      )
    )
    end(binned_cov)[length(binned_cov)-1] <- end(binned_cov)[
      length(binned_cov)-1
      ]-2
  }
  
  # ensure coverage coordinates are within xlims:
  gr_sub <- GRanges(
    seqnames = seqnames(binned_cov)[1],
    ranges = IRanges(
      start = snapshot_list[[i]]$start_coord, 
      end = snapshot_list[[i]]$end_coord
    ),
    strand="*"
  )
  
  binned_cov <- subsetByOverlaps(binned_cov, gr_sub)
  
  cov_track <- DataTrack(
    binned_cov,
    name = "Mean-centered coverage",
    col = "#B066B2",
    background.title = "#B066B2",
    ylim = c(-80, 80)
  )
  
  ########################################################################
  ### 5. Save track/granges for each sample/gene combo ###
  ########################################################################
  
  if (i==1) {
    
    all_track_gr <- list(
      list(
        cov_track = cov_track,
        RT_df = RT_df,
        bp_track = bp_track
      )
    )
    
  } else {
    
    all_track_gr[[i]] <- list(
      cov_track = cov_track,
      RT_df = RT_df,
      bp_track = bp_track
    )
    
  }
 
  if (snapshot_list[[i]]$addstring == "zoomed_out") {
    names(all_track_gr)[i] <- paste0(
      snapshot_list[[i]]$gene_symbol, "_", 
      snapshot_list[[i]]$sample, "_",
      snapshot_list[[i]]$addstring
    )
  } else {
    names(all_track_gr)[i] <- paste0(
      snapshot_list[[i]]$gene_symbol, "_", 
      snapshot_list[[i]]$sample
    )
  }
  
  
  ########################################################################
  ### 6. Plot snapshots ###
  ########################################################################
  
  plotTracks(
    list(
      itrack,
      gtrack,
      cov_track,
      RT_track,
      bp_track,
      grtrack
    ),
    from = snapshot_list[[i]]$start_coord, to = snapshot_list[[i]]$end_coord,
    type = "s"
  )
  
  if (snapshot_list[[i]]$addstring != "none") {
    
    filepath <- paste0(
      plot_dir, snapshot_list[[i]]$sample, "_", 
      snapshot_list[[i]]$gene_symbol, "_snapshot_", 
      snapshot_list[[i]]$addstring
    )
    
  } else {
    
    filepath <- paste0(
      plot_dir, snapshot_list[[i]]$sample, "_", 
      snapshot_list[[i]]$gene_symbol, "_snapshot"
    )
    
  }
  
  pdf(
    paste0(filepath, ".pdf"),
    height = 4,
    width = 8
  )
  plotTracks(
    list(
      itrack, 
      gtrack, 
      cov_track,
      RT_track,
      bp_track,
      grtrack
    ),
    from = snapshot_list[[i]]$start_coord, to = snapshot_list[[i]]$end_coord,
    type = "s"
  )
  dev.off()
  
  png(
    paste0(filepath, ".png"),
    height = 4,
    width = 8,
    res = 300,
    units = "in"
  )
  plotTracks(
    list(
      itrack, 
      gtrack, 
      cov_track,
      RT_track,
      bp_track,
      grtrack
    ),
    from = snapshot_list[[i]]$start_coord, to = snapshot_list[[i]]$end_coord,
    type = "s"
  )
  dev.off()
  
}


########################################################################
### 7. Combine EPHA3 tracks ###
########################################################################

# select only elements with EPHA3_AOCS in name without zoomed out plots:
EPHA3_track_gr <- all_track_gr[
  grep("EPHA3_AOCS", names(all_track_gr))
]
EPHA3_track_gr <- EPHA3_track_gr[
  grep("zoomed_out", names(EPHA3_track_gr), invert = T)
]

# combine RT tracks:
RT_combined <- lapply(EPHA3_track_gr, function(x) x$RT_df)
RT_combined <- do.call("rbind", RT_combined)
RT_combined <- GRanges(
  seqnames = Rle(RT_combined$chromosome),
  ranges = IRanges(
    start = RT_combined$start,
    end = RT_combined$end
  ),
  strand = "*",
  symbol = RT_combined$symbol
)
RT_combined <- RT_combined[!duplicated(RT_combined)]

combined_RT_track <- GeneRegionTrack(
  RT_combined, 
  genome = "hg38", 
  chromosome = "chr3",
  name = "Rt",
  transcriptAnnotation = "symbol",
  col = "#1B9E77",
  fill = "#1B9E77",
  background.title = "#1B9E77"
)

# rename coverage tracks:
for (j in 1:length(EPHA3_track_gr)) {
  EPHA3_track_gr[[j]]$cov_track@name <- paste0("Tumour ", j, " mean-centered coverage")
}

# rename breakpoint tracks:
for (j in 1:length(EPHA3_track_gr)) {
  EPHA3_track_gr[[j]]$bp_track@name <- paste0("Tumour ", j, " bp")
}

plotTracks(
  list(
    itrack,
    gtrack,
    combined_RT_track,
    EPHA3_track_gr[[1]]$cov_track,
    EPHA3_track_gr[[1]]$bp_track,
    EPHA3_track_gr[[2]]$cov_track,
    EPHA3_track_gr[[2]]$bp_track,
    EPHA3_track_gr[[3]]$cov_track,
    EPHA3_track_gr[[3]]$bp_track,
    EPHA3_track_gr[[4]]$cov_track,
    EPHA3_track_gr[[4]]$bp_track,
    grtrack
  ),
  from = snapshot_list[[i]]$start_coord, to = snapshot_list[[i]]$end_coord,
  type = "s"
)

# write plots:
filepath <- paste0(
  plot_dir, snapshot_list[[i]]$gene_symbol, "_combined_snapshot"
)

pdf(
  paste0(filepath, ".pdf"),
  height = 9,
  width = 8
)
plotTracks(
  list(
    itrack,
    gtrack,
    combined_RT_track,
    EPHA3_track_gr[[1]]$cov_track,
    EPHA3_track_gr[[1]]$bp_track,
    EPHA3_track_gr[[2]]$cov_track,
    EPHA3_track_gr[[2]]$bp_track,
    EPHA3_track_gr[[3]]$cov_track,
    EPHA3_track_gr[[3]]$bp_track,
    EPHA3_track_gr[[4]]$cov_track,
    EPHA3_track_gr[[4]]$bp_track,
    grtrack
  ),
  from = 89000000, to = 89600000,
  type = "s"
)
dev.off()

png(
  paste0(filepath, ".png"),
  height = 9,
  width = 8,
  res = 300,
  units = "in"
)
plotTracks(
  list(
    itrack,
    gtrack,
    combined_RT_track,
    EPHA3_track_gr[[1]]$cov_track,
    EPHA3_track_gr[[1]]$bp_track,
    EPHA3_track_gr[[2]]$cov_track,
    EPHA3_track_gr[[2]]$bp_track,
    EPHA3_track_gr[[3]]$cov_track,
    EPHA3_track_gr[[3]]$bp_track,
    EPHA3_track_gr[[4]]$cov_track,
    EPHA3_track_gr[[4]]$bp_track,
    grtrack
  ),
  from = 89000000, to = 89600000,
  type = "s"
)
dev.off()


