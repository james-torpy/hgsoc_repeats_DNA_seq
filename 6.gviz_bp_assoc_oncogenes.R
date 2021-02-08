
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
#home_dir <- "/Users/jamestorpy/clusterHome/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/")

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

# load pre-generated gene models of specific genes:
goi <- c(
  "EPHA3", "ANKRD11", "MECOM", "EXT2",
  "ATM", "POLD1", "CXCR2", "UBE2C", 
  "DKK2", "FANCD2", "BARD1", "CCDC6",
  "FBXW7"
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
    sample = "AOCS-064",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = 89050000,
    end_coord = 89500000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-112",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = 89050000,
    end_coord = 89500000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-122",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = 89050000,
    end_coord = 89500000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-130",
    gene_symbol = "EPHA3",
    chr = "chr3",
    start_coord = 89050000,
    end_coord = 89500000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-075",
    gene_symbol = "ANKRD11",
    chr = "chr16",
    start_coord = 89240000,
    end_coord = 89500000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-075",
    gene_symbol = "MECOM",
    chr = "chr3",
    start_coord = 168900000,
    end_coord = 169700000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-076",
    gene_symbol = "EXT2",
    chr = "chr11",
    start_coord = 44050000,
    end_coord = 44250000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-080",
    gene_symbol = "ATM",
    chr = "chr11",
    start_coord = 108270000,
    end_coord = 108440000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-083",
    gene_symbol = "POLD1",
    chr = "chr19",
    start_coord = 50380000,
    end_coord = 50420000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-085",
    gene_symbol = "CXCR2",
    chr = "chr2",
    start_coord = 218120000,
    end_coord = 218155000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-114",
    gene_symbol = "UBE2C",
    chr = "chr20",
    start_coord = 45810000,
    end_coord = 45818000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-130",
    gene_symbol = "DKK2",
    chr = "chr4",
    start_coord = 106910000,
    end_coord = 107040000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-128",
    gene_symbol = "FANCD2",
    chr = "chr3",
    start_coord = 10020000,
    end_coord = 10105000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-130",
    gene_symbol = "BARD1",
    chr = "chr2",
    start_coord = 214710000,
    end_coord = 214820000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-130",
    gene_symbol = "CCDC6",
    chr = "chr10",
    start_coord = 59770000,
    end_coord = 59910000,
    stringsAsFactors = F
  ),
  data.frame(
    sample = "AOCS-137",
    gene_symbol = "FBXW7",
    chr = "chr4",
    start_coord = 152300000,
    end_coord = 152550000,
    stringsAsFactors = F
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
  gtfs$bp_windows <- import(
    paste0(
      in_dir, snapshot_list[[i]]$sample, "/", snapshot_list[[i]]$sample, 
      "_retrotransposon_assoc_breakpoint_windows.sorted.gtf"
    )
  )
  
  # remove duplicates:
  gtfs <- lapply(gtfs, function(x) x[!duplicated(x)])
  
  # reduce breakpoint windows:
  gtfs$bp_windows <- reduce(gtfs$bp_windows)
  
  # adjust coords for visualisation if needed:
  if (snapshot_list[[i]]$sample == "AOCS-080" & 
    snapshot_list[[i]]$gene_symbol == "ATM") {
    
    ranges(gtfs$bp)[30] <- IRanges(start = 108348000, end = 108348000)
    ranges(gtfs$bp)[31] <- IRanges(start = 108349500, end = 108349500)
    
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
  
  if (
    !file.exists(
      paste0(
        bam_dir, snapshot_list[[i]]$sample, 
        "-1.", snapshot_list[[i]]$gene_symbol, 
        ".coverage.bed"
      )
    )
  ) {
    
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
    system(
      paste0(
        "bedtools intersect -a ", snapshot_list[[i]]$sample, 
        "-1.coverage.bed -b ", snapshot_list[[i]]$sample, 
        "-1.temp.bed > ", snapshot_list[[i]]$sample, "-1.", 
        snapshot_list[[i]]$gene_symbol, ".coverage.bed"
      )
    )
    
    
    # remove temp bedfile:
    system(
      paste0(
        "rm ", bam_dir, snapshot_list[[i]]$sample, "-1.temp.bed")
    )
    
  }
  
  # load bed file:
  cov <- import(
    paste0(
      bam_dir, snapshot_list[[i]]$sample, "-1.", 
      snapshot_list[[i]]$gene_symbol, ".coverage.bed"
    )
  )
  colnames(mcols(cov)) <- "coverage"

  # centre values on mean coverage:
  cov$centered_coverage <- as.numeric(cov$coverage) - mean_cov
  mcols(cov) <- subset(mcols(cov), select = centered_coverage)
  
  # bin coverage values and calculate means:
  # create bins:
  bin_vec <- snapshot_list[[i]]$start_coord:snapshot_list[[i]]$end_coord
  bin_list <- split(bin_vec, ceiling(seq_along(bin_vec)/20000))
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
  
  covTrack <- DataTrack(binned_cov, name = "Mean-centered coverage")
  plotTracks(covTrack, type = "s")
  
  
  ########################################################################
  ### 5. Plot snapshots ###
  ########################################################################
  
  plotTracks(
    list(
      itrack, 
      gtrack, 
      covTrack,
      RT_track,
      bp_track,
      grtrack
    ),
    from = snapshot_list[[i]]$start_coord, to = snapshot_list[[i]]$end_coord,
    type = "s"
  )
  
  pdf(
    paste0(
      plot_dir, snapshot_list[[i]]$sample, "_", 
      snapshot_list[[i]]$gene_symbol, "_snapshot.pdf"),
    height = 3,
    width = 8
  )
    plotTracks(
      list(
        itrack, 
        gtrack, 
        RT_track,
        bp_track,
        grtrack
      ),
      from = snapshot_list[[i]]$start_coord, to = snapshot_list[[i]]$end_coord
    )
  dev.off()
  
  png(
    paste0(plot_dir, snapshot_list[[i]]$sample, "_", 
      snapshot_list[[i]]$gene_symbol, "_snapshot.png"),
    height = 3,
    width = 8,
    res = 300,
    units = "in"
  )
    plotTracks(
      list(
        itrack, 
        gtrack, 
        RT_track,
        bp_track,
        grtrack
      ),
      from = snapshot_list[[i]]$start_coord, to = snapshot_list[[i]]$end_coord
    )
  dev.off()
  
}






