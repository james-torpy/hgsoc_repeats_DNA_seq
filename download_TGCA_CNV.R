
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
ref_dir <- paste0(project_dir, "refs/")
out_path <- paste0(ref_dir, "TGCA/OV_CNV/")
raw_dir <- paste0(out_path, "/raw_files/")
func_dir <- paste0(project_dir, "scripts/functions/")

incl_non_hg38 <- TRUE

if (incl_non_hg38) {
  out_dir <- paste0(out_path, "/CNV_data/all_genome_builds/")
} else {
  out_dir <- paste0(out_path, "/CNV_data/hg38_only/")
}

dir.create(out_dir, recursive = T)

library(dplyr)
library(RTCGA)
library(GenomicRanges)

if (incl_non_hg38) {

  # load functions:
  do_liftover <- dget(paste0(func_dir, "do_liftover.R"))
  collapse_lifted_CNA <- dget(paste0(func_dir, "collapse_lifted_CNA.R"))

}

# check available datasets:
checkTCGA('Dates')

# check available cohorts:
(cohorts <- infoTCGA() %>% 
   rownames() %>% 
   sub("-counts", "", x=.))

# check available datasets:
avail_ds <- checkTCGA("DataSets", "OV", date = "2016-01-28")
CNV_ds <- avail_ds[grep("CNA", avail_ds$Name, ignore.case = T),]

# specify dataset details:
releaseDate <- "2016-01-28"
cohort <- "OV"
datasets <- CNV_ds$Name

# download data:
all_cna_files <- list.files(
  raw_dir, pattern = "cna", recursive = TRUE, full.names = TRUE
)

if (length(all_cna_files) == 0) {
  for (dataset in datasets) {
    try(downloadTCGA(
      cancerTypes = cohort, 
      destDir = raw_dir, 
      date = releaseDate, 
      dataSet = dataset
      ),
      silent=TRUE
    )
  }
}

# read and format data:
all_cna_files <- list.files(
  raw_dir, pattern = "cna.*txt", recursive = TRUE, full.names = TRUE
)

if (!incl_non_hg38) {
  all_cna_files <- grep("hg18|hg19", all_cna_files, invert = TRUE, value = TRUE)
}

for (c in 1:length(all_cna_files)) {

  print(paste0("Adding CNA co-ordinates from ", all_cna_files[c]))

  # load CNAs:
  CNA <- read.table(all_cna_files[c], header = T)

  # change chromosome 23 and 24 to X and Y:
  CNA$Chromosome[CNA$Chromosome == "23"] <- "X"

  # convert to granges object:
  CNA_gr <- GRanges(
    seqnames = Rle(paste0("chr", CNA$Chromosome)),
    ranges = IRanges(
      start = CNA$Start, 
      end = CNA$End
    ),
    strand = Rle("*")
  )

  print(paste0("Length of granges object is: ", length(CNA_gr)))

  if ( length(grep("hg18", all_cna_files[c])) > 0 ) {

    new_filename <- gsub("(.*)hg18.*", "\\1hg38.Rdata", all_cna_files[c])
    
    if (!file.exists(new_filename)) {

      print("Co-ordinates mapped to hg18, converting to hg38...")

      # liftover to hg38:
      CNA_gr_hg38 <- do_liftover(CNA_gr, "hg18", "hg38")
  
      # for each element, keep only min and max coordinates of main chromosome:
      CNA_gr_hg38 <- lapply(CNA_gr_hg38, collapse_lifted_CNA)
      CNA_gr_hg38 <- do.call("c", CNA_gr_hg38)
  
      saveRDS(CNA_gr_hg38, new_filename)

    } else {

      print("Loading hg38 version and adding...")
      CNA_gr_hg38 <- readRDS(new_filename)

    }

    CNA_gr <- CNA_gr_hg38

  } else if ( length(grep("hg19", all_cna_files[c])) > 0 ) {

    new_filename <- gsub("(.*)hg19.*", "\\1hg38.Rdata", all_cna_files[c])
    
    if (!file.exists(new_filename)) {

      print("Co-ordinates mapped to hg18, converting to hg38...")
  
      # liftover to hg38:
      CNA_gr_hg38 <- do_liftover(CNA_gr, "hg19", "hg38")
  
      # for each element, keep only min and max coordinates of main chromosome:
      CNA_gr_hg38 <- lapply(CNA_gr_hg38, collapse_lifted_CNA)
      CNA_gr_hg38 <- do.call("c", CNA_gr_hg38)
  
      saveRDS(CNA_gr_hg38, new_filename)

    } else {

      print("Loading hg38 version...")
      CNA_gr_hg38 <- readRDS(new_filename)

    }

    CNA_gr <- CNA_gr_hg38

  }

  if (c==1) {
    all_CNA <- CNA_gr
  } else {
    all_CNA <- c(all_CNA, CNA_gr)
  }

}

# convert to breakpoints:
bp_temp <- GRanges(
  seqnames = seqnames(all_CNA),
  ranges = IRanges(
    start = start(all_CNA), 
    end = start(all_CNA)
  ),
  strand = Rle("*")
)

all_bp <- c(
  bp_temp,
  GRanges(
    seqnames = seqnames(all_CNA),
    ranges = IRanges(
      start = end(all_CNA), 
      end = end(all_CNA)
    ),
    strand = Rle("*")
  )
)

# reduce to non-overlapping co-ordinates:
all_bp_reduced <- reduce(all_bp)

# save as RDS:
saveRDS(all_bp_reduced, paste0(out_dir, "TGCA_OV_breakpoints.Rdata"))

