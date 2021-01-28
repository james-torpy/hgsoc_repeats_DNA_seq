
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
ref_dir <- paste0(project_dir, "refs/")
out_path <- paste0(ref_dir, "TGCA/OV_CNV/")
raw_dir <- paste0(out_path, "/raw_files/")
out_dir <- paste0(out_path, "/CNV_data/")

dir.create(out_dir, recursive = T)

library(dplyr)
library(RTCGA)

# create liftover function:
do_liftover <- function(gr, original, output) {

  library(rtracklayer)

  output <- gsub("h", "H", output)
  chain <- import.chain(paste0(ref_dir, "liftover/", original, 
    "To", output, ".over.chain"))

  return(liftOver(gr, chain))

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

# read and format data:
all_cna_files <- list.files(
  raw_dir, pattern = "cna", recursive = TRUE, full.names = TRUE
)


for (c in 1:length(all_cna_files)) {

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

  if ( length(grep("hg18", all_cna_files[c])) > 0 ) {
    CNA_gr_test <- do_liftover(CNA_gr, "hg18", "hg38")
  }

  if ( length(grep("hg19", all_cna_files[c])) > 0 ) {
    CNA_gr_test <- do_liftover(CNA_gr, "hg19", "hg38")
  }

  if (c==1) {
    all_CNA <- CNA
  } else {
    all_CNA <- rbind(all_CNA, CNA)
  }

}
   


