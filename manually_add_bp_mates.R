AOCS_075_bp <- import(
  paste0(table_dir, "AOCS-075/AOCS-075_retrotransposon_assoc_breakpoints.sorted.gtf")
)

cancer_assoc[cancer_assoc$symbol == "ANKRD11",]

to_add <- total_bp[total_bp$id == "175056237:1"]
to_add_gr <- GRanges(
  seqnames = seqnames(to_add),
  ranges = IRanges(start = to_add$bp, end = to_add$bp),
  strand = "*",
  bp_id = to_add$id,
  CNV_type = to_add$CNV_type,
  sample = to_add$sample
)

res <- c(AOCS_075_bp, to_add_gr)

export(
  res,
  paste0(table_dir, "AOCS-075/AOCS-075_retrotransposon_assoc_breakpoints.sorted.gtf")
)