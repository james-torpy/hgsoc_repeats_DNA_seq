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
  )[,c(1:5,8)]
  colnames(manta_CNV) <- c("chr", "pos", "id", "joining_nt", "join", "info")

  # fetch breakpoints at the start of CNVs:
  manta_CNV_only <- manta_CNV[grep("SVTYPE=BND", manta_CNV$info, invert = T),]
  manta_CNV_bp1 <- subset(manta_CNV_only, select = c(chr, pos, id, joining_nt, join, info))
  manta_CNV_bp1$type <- gsub(
    ";.*$", "",
    gsub("^.*SVTYPE=", "", manta_CNV_bp1$info)
  )

  # fetch breakpoints at the end of CNVs:
  manta_CNV_bp2 <- subset(manta_CNV_only, select = c(chr, id, joining_nt, join, info))
  manta_CNV_bp2 <- data.frame(
    chr = manta_CNV_bp2$chr,
    pos = gsub(
      ";.*$", "",
      gsub("END=", "", manta_CNV_bp2$info)
    ),
    id = manta_CNV_bp2$id,
    joining_nt = manta_CNV_bp2$joining_nt,
    join = manta_CNV_bp2$join,
    info = manta_CNV_bp2$info,
    type = gsub(
      ";.*$", "",
      gsub("^.*SVTYPE=", "", manta_CNV_bp2$info)
    )
  )

  # merge:
  manta_CNV_bp <- rbind(manta_CNV_bp1, manta_CNV_bp2)

  # fetch other breakpoints:
  manta_bp <- manta_CNV[grep("BND", manta_CNV$info),]
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
      id = manta_bp$id,
      CNV_type = manta_bp$type,
      sample = samp_name,
      joining_nt = manta_bp$joining_nt,
      join_nt = manta_bp$join,
      bp_info = manta_bp$info
    )
  )

}