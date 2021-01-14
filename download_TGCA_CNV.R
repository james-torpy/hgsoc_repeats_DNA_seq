
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
ref_dir <- paste0(project_dir, "refs/")
out_path <- paste0(ref_dir, "TGCA/OV_CNV/")
raw_dir <- paste0(out_path, "/raw_files/")
out_dir <- paste0(out_path, "/CNV_data/")

dir.create(out_dir, recursive = T)

library(RTCGA)

# check available datasets:
checkTCGA('Dates')

# check available cohorts:
(cohorts <- infoTCGA() %>% 
   rownames() %>% 
   sub("-counts", "", x=.))

# check available datasets:
avail_ds <- checkTCGA("DataSets", "OV", date = "2016-01-28")
CNV_ds <- avail_ds[grep("CNV|CNA|gistic", avail_ds$Name, ignore.case = T),]


# specify dataset details:
releaseDate <- "2016-01-28"
cohort <- "OV"
datasets <- c("")

# download data:
for (cohort in cohorts) {
  try(downloadTCGA( cancerTypes = cohort, destDir = raw_dir, 
  	date = releaseDate, 
  	dataSet = "" ),
    silent=TRUE
  )
}

# read and format data:
allCNVFiles <- list.files(
  raw_dir, pattern = "cnv", recursive = TRUE, full.names = TRUE
)

for (CNVFile in allCNVFiles) {

  CNV <- read.table(CNVFile,h=T) 
   
  cohortName <- strsplit(
  	strsplit(CNVFile, split = "/"
  )[[1]][4], "\\.")[[1]][1]
  name = paste0(cohortName, ".CNV")
  assign(name, CNV)
  save(
    list = name, 
    file=paste0(out_dir, name, ".rda"), 
    compression_level = 9, 
    compress = "xz"
  )

}


