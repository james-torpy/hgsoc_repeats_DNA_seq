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

  # fetch breakpoints at the start of CNVs:
  manta_CNV_only <- manta_CNV[grep("SVTYPE=BND", manta_CNV$type, invert = T),]
  manta_CNV_bp1 <- subset(manta_CNV_only, select = c(chr, pos, type))
  manta_CNV_bp1$type <- gsub(
    ";.*$", "",
    gsub("^.*SVTYPE=", "", manta_CNV_bp1$type)
  )

  # fetch breakpoints at the end of CNVs:
  manta_CNV_bp2 <- subset(manta_CNV_only, select = c(chr, type))
  manta_CNV_bp2 <- data.frame(
    chr = manta_CNV_bp2$chr,
    pos = gsub(
      ";.*$", "",
      gsub("END=", "", manta_CNV_bp2$type)
    ),
    type = gsub(
      ";.*$", "",
      gsub("^.*SVTYPE=", "", manta_CNV_bp2$type)
    )
  )

  # merge:
  manta_CNV_bp <- rbind(manta_CNV_bp1, manta_CNV_bp2)

  # fetch other breakpoints:
  manta_bp <- manta_CNV[grep("BND", manta_CNV$type),]
  manta_bp$type <- "BND"

  # join to other breakpoints:
  manta_bp <- rbind(
    manta_bp,
    manta_CNV_bp
  )

  print(
    paste0(
      "No. manta entries before splitting CNVs to breakpoints = ", 
      nrow(manta_CNV)
    )
  )
  
  print(
    paste0(
      "No. manta breakpoints = ", 
      nrow(manta_bp)
    )
  )
  
#  # create column with both chr and pos info:
#  manta_bp$joint <- paste0(manta_bp$chr, "_", manta_bp$pos)

  # create GRanges of bps with 100 base pair window:
  manta_bp$pos <- as.numeric(manta_bp$pos)
  
  return(
    GRanges(
      seqnames = Rle(manta_bp$chr),
      ranges = IRanges(
        start = manta_bp$pos-(bp_window/2), 
        end = manta_bp$pos+(bp_window/2)
      ),
      strand = Rle("*"),
      bp = manta_bp$pos,
      CNV_type = manta_bp$type,
      sample = samp_name
    )
  )

}