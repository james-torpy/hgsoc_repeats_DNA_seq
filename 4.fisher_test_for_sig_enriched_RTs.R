
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")

in_dir <- paste0(results_dir, "bp_assoc_RT/")
Robject_dir <- paste0(in_dir, "/Rdata/")

svaba_dir <- paste0(results_dir, "svaba/")
manta_dir <- paste0(results_dir, "manta/")

create_svaba_granges <- dget(paste0(func_dir, "create_svaba_granges.R"))
create_manta_granges <- dget(paste0(func_dir, "create_manta_granges.R"))

overlap_window <- 10000


########################################################################
### 0. Load packages and functions ###
########################################################################

library(GenomicRanges)

create_svaba_granges <- dget(paste0(func_dir, "create_svaba_granges.R"))


########################################################################
### 1. Load repetitive elements ###
########################################################################

# This is a list, with each memeber representing a GRanges object of 
# location:
if (!file.exists(paste0(Robject_dir, "repeat_gtf.Rdata"))) {

  repeat_gtf <- rtracklayer::import(
    paste0(genome_dir, "repeats.hg38.gtf")
  )
  
  # remove alternate chrs:
  repeat_gtf <- repeat_gtf[
    grep("K|G|M", seqnames(repeat_gtf), invert = T)
  ]

  saveRDS(repeat_gtf, paste0(Robject_dir, "repeat_gtf.Rdata"))

} else {
  repeat_gtf <- readRDS(paste0(Robject_dir, "repeat_gtf.Rdata"))
}

repeats <- GRangesList(split(repeat_gtf, repeat_gtf$type))
repeats <- repeats[c("L1MD3","L1PA2","L1HS", "AluYb8", "AluYh3")]


########################################################################
### 2. Load breakpoints ###
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

# fetch low confidence breakpoints present in both svaba and manta
# outputs:
if (
  !file.exists(
    paste0(Robject_dir, "consensus_low_conf_breakpoints.Rdata")
  )
) {

  for (s in 1:length(sample_names)) {
    
    print(
      paste0(
        "Fetching and verifying breakpoints for ", sample_names[s]
      )
    )
  
    svaba_gr <- create_svaba_granges(
      samp_name = sample_names[s],
      in_dir = svaba_dir,
      SV_type = "unfiltered",
      bp_window = 100
    )
  
    manta_gr <- create_manta_granges(
      samp_name = sample_names[s],
      in_dir = manta_dir,
      SV_type = "unfiltered",
      bp_window = 100
    )
    
    # save bps that occur within 100 bp window from filtered svaba
    # and unfiltered manta:
    svaba_olaps <- findOverlaps(svaba_gr, manta_gr)
    bp <- svaba_gr[unique(queryHits(svaba_olaps))]
  
    # extand ranges to n kb each side of bp:
    ranges(bp) <- IRanges(
      start = bp$bp - overlap_window,
      end = bp$bp + overlap_window
    )

    print(
      paste0(
        "No. breakpoints for ", sample_names[s], " = ", length(bp)
      )
    )

    if (s==1) {
      low_conf_bp <- bp
    } else {
      low_conf_bp <- c(low_conf_bp, bp)
    }

    print(
      paste0(
        "Total No. breakpoints so far = ", length(low_conf_bp)
      )
    )

    writeLines("\n")

  }

  # remove those coordinates that overlap:



  

  saveRDS(
    low_conf_bp, 
    paste0(Robject_dir, "consensus_low_conf_breakpoints.Rdata")
  )

} else {

  low_conf_bp <- readRDS(
    paste0(Robject_dir, "consensus_low_conf_breakpoints.Rdata")
  )

}


# Sanity check:
# Just to check that statistics are working as they should, I 
# would add a positive control
# Basically a subset of target break regions, or a mix of break 
# regions and random regions

# load genome breakpoints:
all_bp <- readRDS(paste0(Robject_dir, "all_breakpoints.Rdata"))



for (i in 1:length(all_bp)) {

  target <- all_bp[[i]]

  repeats[["control"]] <- target[1:100]


  ########################################################################
  ### 2. Fetch a list of regions that you want to look into overlaps 
  # with ###
  ########################################################################
  
  breakPoints <- GRangesList()
  breakPoints[["all"]] <- 
  breakPoints[["target"]] <- target
  breakPoints[["random"]] <- grRandom

  #make sure you have a) all regions b) target regions per sample c) control regions with ~1000 random regions (let me know if you need script for this) in that object
  breakPoints<-GRangesList()
  breakPoints[["all"]]<-grAll
  breakPoints[["target"]]<-target
  breakPoints[["random"]]<-grRandom
  #the analysis is performed using fisher test
  results<-list()
  counts<-list()
  output<-list()
  for(repeatName in names(repeats)){
    results[[repeatName]]<-sapply(breakPoints,function(x){sum(countOverlaps(x,repeats[[repeatName]],ignore.strand=T)>0)/length(x)})
    counts[[repeatName]]<-sapply(breakPoints,function(x){sum(countOverlaps(x,repeats[[repeatName]],ignore.strand=T)>0)})
    n=length(breakPoints[["all"]])
    nTarget=length(breakPoints[["target"]])
    upR<-counts[[repeatName]][3]
    upNR<-nTarget-counts[[repeatName]][3]
    NupR<-counts[[repeatName]][5]-counts[[repeatName]][3]
    NupNR<-n-nTarget-counts[[repeatName]][5]
    output[[repeatName]]<-fisher.test(matrix(c(upR,upNR,NupR,NupNR),nrow=2),alternative="greater")
  }
  dfC<-do.call("rbind",counts)
  #check the results in the output object

}
