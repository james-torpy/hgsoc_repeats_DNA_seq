
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/map_and_count/functions/")

svaba_dir <- paste0(results_dir, "svaba/")
manta_dir <- paste0(results_dir, "manta/")
col_dir <- paste0(home_dir, "R/colour_palettes/")


########################################################################
### 0.  Load packages and functions ###
########################################################################

library(GenomicRanges)


########################################################################
### 1.  Load and format VCFs ###
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

create_svaba_granges <- function(
  samp_name,
  in_dir,
  SV_type = "filtered",
  bp_window = 50
) {

  # load svaba vcf and format:
  if (SV_type == "filtered") {
    fname <- ".svaba.somatic.sv.vcf"
  } else {
    fname <- ".svaba.unfiltered.somatic.sv.vcf"
  }

  svaba_bp <- read.table(
    paste0(in_dir, sample_names[s], "/", sample_names[s], fname)
  )[,c(1:2,8)]
  colnames(svaba_bp) <- c("chr", "pos", "type")
  svaba_bp$type <- gsub("^.*;SVTYPE=", "", svaba_bp$type)
  
  # isolate breakpoints
  print(
    paste0(
      "No. svaba entries before filtering for breakpoints only and ",
      "removing non-standard chrs = ", 
      nrow(svaba_bp)
    )
  )
  
  svaba_bp <- svaba_bp[svaba_bp$type == "BND",]
  
  # remove non-standard chromosomes:
  svaba_bp <- svaba_bp[grep("K|G|M", svaba_bp$chr, invert = T),]
  
  print(
    paste0(
      "No. svaba breakpoints = ", 
      nrow(svaba_bp)
    )
  )
  
  # create GRanges of bps with 50 base pair window:
  return(
    GRanges(
      seqnames = Rle(svaba_bp$chr),
      ranges = IRanges(
        start = svaba_bp$pos-(bp_window/100), 
        end = svaba_bp$pos+(bp_window/100)
      ),
      strand = Rle("*"),
      bp = svaba_bp$pos
    )
  )

}

create_manta_granges <- function(
  samp_name,
  in_dir,
  SV_type = "filtered",
  bp_window = 50
) {

  # load manta vcf, format and isolate breakpoints:
  if (SV_type == "filtered") {
    fname <- "somaticSV.vcf.gz"
  } else {
    fname <- "candidateSV.vcf.gz"
  }

  manta_CNV <- read.table(
    paste0(in_dir, samp_name, "/results/variants/", fname)
  )[,c(1:2,8)]
  colnames(manta_CNV) <- c("chr", "pos", "type")
  manta_CNV$type <- gsub(";.*", "", manta_CNV$type)
  
  # fetch breakpoints at the ends of CNVs:
  manta_CNV_only <- manta_CNV[grep("BND", manta_CNV$type, invert = T),]
  manta_CNV_bp1 <- subset(manta_CNV_only, select = c(chr, pos))
  manta_CNV_bp2 <- subset(manta_CNV_only, select = c(chr, type))
  colnames(manta_CNV_bp2) <- colnames(manta_CNV_bp1)
  manta_CNV_bp <- rbind(manta_CNV_bp1, manta_CNV_bp2)
  manta_CNV_bp$pos <- gsub("END=", "", manta_CNV_bp$pos)
  manta_CNV_bp$type <- "BND"

  # add to other breakpoints:
  manta_bp <- manta_CNV[grep("BND", manta_CNV$type),]
  manta_bp$type <- "BND"

  # joint to other breakpoints:
  manta_bp <- rbind(
    manta_bp,
    manta_CNV_bp
  )

  print(
    paste0(
      "No. manta entries before filtering for breakpoints only = ", 
      nrow(manta_bp)
    )
  )
  manta_bp <- manta_bp[grep("BND", manta_bp$type),]
  print(
    paste0(
      "No. manta breakpoints = ", 
      nrow(manta_bp)
    )
  )
  
#  # create column with both chr and pos info:
#  manta_bp$joint <- paste0(manta_bp$chr, "_", manta_bp$pos)

  # create GRanges of bps with 50 base pair window:
  manta_bp$pos <- as.numeric(manta_bp$pos)
  
  return(
    GRanges(
      seqnames = Rle(manta_bp$chr),
      ranges = IRanges(
        start = manta_bp$pos-(bp_window/100), 
        end = manta_bp$pos+(bp_window/100)
      ),
      strand = Rle("*"),
      bp = manta_bp$pos
    )
  )

}

for (s in 1:length(sample_names)) {

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
  olaps1 <- findOverlaps(svaba_gr_filt, manta_gr)
  print(length(olaps1))

  # save bps that occur within 50 bp window from filtered manta
  # and unfiltered svaba:
  olaps2 <- findOverlaps(manta_gr_filt, svaba_gr)
  print(length(olaps2))

}







# load VCFs:
lapply(sample_names, function(x) {

  svaba_CNV <- read.table(
      paste0(svaba_dir, x, "/", x, ".svaba.unfiltered.somatic.sv.vcf")
    )
  manta_CNV <- read.table(
    paste0(manta_dir, x, "/results/variants/somaticSV.vcf.gz")
  )

})





