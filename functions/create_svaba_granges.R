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
  )[,c(1:5,8)]
  colnames(svaba_bp) <- c("chr", "pos", "id", "joining_nt", "join", "info")
  svaba_bp$type <- gsub("^.*;SVTYPE=", "", svaba_bp$info)
  
  # isolate breakpoints
  print(
    paste0(
      "No. svaba entries before filtering for breakpoints only and ",
      "removing non-standard chrs = ", 
      nrow(svaba_bp)
    )
  )
  
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
        start = svaba_bp$pos-(bp_window/2), 
        end = svaba_bp$pos+(bp_window/2)
      ),
      strand = Rle("*"),
      bp = svaba_bp$pos,
      id = svaba_bp$id,
      CNV_type = svaba_bp$type,
      sample = samp_name,
      joining_nt = svaba_bp$joining_nt,
      join = svaba_bp$join,
      bp_info = svaba_bp$info
    )
  )

}