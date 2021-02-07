
project_name <- "hgsoc_repeats/DNA-seq"

#home_dir <- "/share/ScratchGeneral/jamtor/"
home_dir <- "/Users/jamestorpy/clusterHome/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/")

in_dir <- paste0(results_dir, "bp_assoc_RT/tables/")

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
### 2. Load sample data ###
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
    start_coord = 108210000,
    end_coord = 108370000,
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
  
  
  ########################################################################
  ### 3. Load remaining tracks ###
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
  ### 4. Plot snapshots ###
  ########################################################################
  
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






