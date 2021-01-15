
project_name <- "hgsoc_repeats/DNA-seq"

home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "projects/", project_name, "/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "refs/")
func_dir <- paste0(project_dir, "scripts/functions/")
genome_dir <- paste0(project_dir, "genome/") 

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
library(rtracklayer)

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
repeats <- repeats[
  c("L1MD3","L1PA2","L1HS", "AluYb8", "AluYh3", "MER66B")
]


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

    # add filtered breakpoints:
    svaba_gr_filt <- create_svaba_granges(
      samp_name = sample_names[s],
      in_dir = svaba_dir,
      SV_type = "filtered",
      bp_window = 100
    )
  
    manta_gr_filt <- create_manta_granges(
      samp_name = sample_names[s],
      in_dir = manta_dir,
      SV_type = "filtered",
      bp_window = 100
    )

    bp <- c(bp, svaba_gr_filt, manta_gr_filt)
  
    # extand ranges to n kb each side of bp:
    ranges(bp) <- IRanges(
      start = start(bp) - overlap_window,
      end = end(bp) + overlap_window
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

  # reduce to remove those coordinates that overlap:
  # remove duplicate overlap records:
  remove_overlaps <- function(gr) {
  
    print("Finding granges object overlaps with itself...")
    ol <- findOverlaps(gr, gr)
  
    # remove matches to themselves:
    print("Removing matches of each range to itself...")
    ol <- as.data.frame(ol[queryHits(ol) != subjectHits(ol)])
  
    # remove duplicates:
    print("Marking duplicates...")
    ol$keep = T
    for (o in 1:nrow(ol)) {

      print(paste0("Checked ", o, " out of ", nrow(ol)))
    
      if (ol$keep[o]) {
        ol$keep[
          ol$queryHits == ol$subjectHits[o] & 
          ol$subjectHits == ol$queryHits[o]
        ] <- F
      }
      
    }
    print("Removing duplicates...")
    ol <- ol[ol$keep,]
  
    # remove subjectHits duplicates:
    print("Removing duplicates of co-ordinates to remove from granges object...")
    ol <- ol[!duplicated(ol$subjectHits),]
  
    # remove remaining subjectHits from gr:
    print("Removing second range of each overlapping pair...")
    gr <- gr[-ol$subjectHits]
  
    return(ol)
  
  }

  low_conf_bp <- remove_overlaps(low_conf_bp)


  saveRDS(
    low_conf_bp, 
    paste0(Robject_dir, "non_overlapping_consensus_low_conf_breakpoints.Rdata")
  )

} else {

  low_conf_bp <- readRDS(
    paste0(Robject_dir, "non_overlapping_consensus_low_conf_breakpoints.Rdata")
  )

}

# Sanity check:
# Just to check that statistics are working as they should, I 
# would add a positive control
# Basically a subset of target break regions, or a mix of break 
# regions and random regions

# load genome breakpoints:
all_bp <- readRDS(paste0(Robject_dir, "all_breakpoints.Rdata"))

# add 10 kb either side:
all_bp_windows <- lapply(all_bp, function(x) {
  start(x) <- start(x) - overlap_window
  end(x) <- end(x) + overlap_window
  return(x)
})

# retrieve chromosome lengths for random region control:
chr_lengths <- as.data.frame(SeqinfoForUCSCGenome("hg38"))
chr_lengths <- chr_lengths[
  grep("K|G|M|J", rownames(chr_lengths), invert = T),
]

# generate random positions:
generateRandomPos <- function(n,chr,chr.sizes,width=1){
  random_chr <- sample(x=chr,size=n,prob=chr.sizes,replace=T)
  random_pos <- sapply(random_chr,function(chrTmp){sample(chr.sizes[chr==chrTmp],1)}) 
  res <- GRanges(random_chr,IRanges(random_pos,random_pos+width))
  return(res)
}

set.seed(632)
random_bp <- generateRandomPos(
  1000, rownames(chr_lengths), chr_lengths$seqlengths
)

#chr_lengths$cum_sum <- cumsum(as.numeric(chr_lengths$seqlengths))
#chr_lengths$cum_start <- c(1, chr_lengths$cum_sum[1:(nrow(chr_lengths)-1)])
#chr_lengths$cum_end <- chr_lengths$cum_sum
#
## create list of cumulative chromosome intervals:
#for (r in 1:nrow(chr_lengths)) {
#  if (r==1) {
#    cum_chr_int <- list(chr_lengths$cum_start[r]:chr_lengths$cum_end[r])
#    names(cum_chr_int) <- rownames(chr_lengths)[r]
#  } else {
#    cum_chr_int[[r]] <- chr_lengths$cum_start[r]:chr_lengths$cum_end[r]
#    names(cum_chr_int) <- rownames(chr_lengths)[r]
#  }
#}
#
## pick 1000 random start positions:
#random_starts <- sample(1:max(chr_lengths$cum_sum), 1000)
#
## assign to chromosomes and chromosomal positions:
#starts <- lapply(cum_chr_int, function(x) {
#  which(random_starts)
#})


for (i in 1:length(all_bp)) {

  repeats[["control"]] <- all_bp[[i]][1:100]

  ########################################################################
  ### 2. Fetch a list of regions that you want to look into overlaps 
  # with ###
  ########################################################################

  # add to granges list: a) all regions b) target regions per sample 
  # c) control regions with ~1000 random regions
  breakPoints <- GRangesList()
  breakPoints[["all"]] <- low_conf_bp
  breakPoints[["target"]] <- all_bp_windows[[i]]
  breakPoints[["random"]] <- random_bp

  #the analysis is performed using fisher test
  results <- list()
  counts <- list()
  output <- list()

  for (repeatName in names(repeats)) {

    # fraction of number of overlaps between set of breakpoints
    # and retrotransposon ranges vs number of breakpoints in set: 
    results[[repeatName]] <- sapply(breakPoints, function(x) {
      sum(
        countOverlaps(
          x, repeats[[repeatName]], ignore.strand=T
        ) > 0
      )/length(x)
    })

    # number of overlaps between set of breakpoints
    # and retrotransposon ranges:
    counts[[repeatName]] <- sapply(breakPoints, function(x) {
      sum(
        countOverlaps(
          x,repeats[[repeatName]], ignore.strand=T
        ) > 0
      )
    })

    # number of total possible breakpoints:
    nAll = length(breakPoints[["all"]])

    # number of detected breakpoints:
    nTarget = length(breakPoints[["target"]])

    bp_RT <- counts[[repeatName]][3]
    bp_no_RT < -nTarget-counts[[repeatName]][3]
    no_bp_RT <- counts[[repeatName]][5]-counts[[repeatName]][3]
    no_bp_no_RT <- n-nTarget-counts[[repeatName]][5]
    output[[repeatName]]<-fisher.test(matrix(c(upR,upNR,NupR,NupNR),nrow=2),alternative=“greater”)
  }
  
  dfC<-do.call("rbind",counts)
  #check the results in the output object

}
