### 4.find_repeat_insertions.R ###

# This script takes discordant reads from svaba output overlaps with a 
# custom DE repeats annotation to find retrotransposon insertions, 
# and filters those lying within breakpoints from CNVnator output within
# following distances: 0, 2, 5, 10, 100, 1000 kb

# module load R/3.5.1 for NCI

library(rtracklayer)
library(GenomicRanges)
library(org.Hs.eg.db)
library(GEOquery)
library("BSgenome.Hsapiens.UCSC.hg19")
library(ggplot2)
library("RColorBrewer")
library(parallel)


### 0. Define variables/paths ###

project <- "hgsoc_repeats"

#home_dir="/g/data1a/ku3/jt3341/"
home_dir <- "/Users/jamestorpy/nciShare/"
project_dir=paste0(home_dir, "/projects/", project, "/genome/")
results_dir=paste0(project_dir, "/results/")
Robject_dir=paste0(project_dir, "/Robjects/")
system(paste0("mkdir -p ", results_dir))
system(paste0("mkdir -p ", Robject_dir))

ref_dir="/Users/jamestorpy/nciShare/genomes/"
#ref_dir="/g/data1a/ku3/jt3341/genomes/"
svaba_dir=paste0(results_dir, "/svaba/")
cnvnator_dir=paste0(results_dir, "/cnvnator/")
table_dir=paste0(results_dir, "/tables/")
plot_dir=paste0(results_dir, "/plots/")
system(paste0("mkdir -p ", table_dir))
system(paste0("mkdir -p ", plot_dir))

#sample_id <- "AOCS-153-1-2_subset"

# define distances from breakpoints in which to look for repeat insertions:
dist <- list(1000, 2000, 5000, 10000, 100000, 1000000)

c("L1MD3", "L1M8", "L1P2", "L1P3", "L1M3f", "L1ME5", "L1MB4", "L1MEg2", 
  "L1M4c", "L1M4b", "L1MA5", "L1MEg", "L1MEd", "L1MEj", "L1MA5A", "L1PB2", 
  "L1PA10", "L1M4a2", "L1P4", "L1M4a1", "L1M3c", "L1P3b", "AluYi6", 
  "AluYd8")

##############################################################################
### 1. Load in svaba discrepencies and repeat annotation ###
##############################################################################

if ( file.exists(paste0(Robject_dir, "/rep_hg19_gr.rds")) ) {
  
  rep_gr <- readRDS(paste0(Robject_dir, "/rep_hg19_gr.rds"))

} else {
  
  rep_annot <- read.table(paste0(ref_dir, "/repeats/repeats.hg19.gff"))
  
  # convert repeat annotation to GRanges object:
  rep_gr <- GRanges(
    seqnames = Rle(rep_annot$V1),
    ranges = IRanges(start=as.integer(rep_annot$V4), 
                     end = as.integer(rep_annot$V5)),
    strand = Rle("*"),
    gene_id = gsub("ID=", "", rep_annot$V9)
  )
  # remove unwanted chromosomes:
  rep_gr <- rep_gr[grep("M|G|K", seqnames(rep_gr), 
                        invert = T)]
  # remove dust, trf, rRNA:
  rep_gr <- rep_gr[!(rep_gr$gene_id %in% c("dust", "trf", "rRNA"))]
  # save rep_gr as RDS file:
  saveRDS(rep_gr, paste0(Robject_dir, "/rep_hg19_gr.rds"))
  
}

if ( file.exists(paste0(Robject_dir, "/discord_list_hg19.rds")) ) {
  
  discord <- readRDS(paste0(Robject_dir, "/discord_list_hg19.rds"))

} else {
  
  # fetch all svaba discordant reads files from svaba_dir:
  discord_files <- grep("subset", list.files(svaba_dir, pattern = "discordant", 
    full.names = T), value = T, invert = T)
  
  # keep only discordant reads matched for which CNVnator output exists:
  cnv_samples <- basename(
    gsub(
      "_CNV_calls.txt", "", grep(
        "subset", list.files(cnvnator_dir, pattern = "CNV_calls", 
        full.names = T), value = T, invert = T
      )
    )
  )
  
  discord_files <- discord_files[
    basename(
      gsub(
        ".gz", "",
        gsub(
          ".discordant.txt", "", discord_files
        )
      )
    ) %in% cnv_samples
  ]
  
  # split into list:
  discord_files <- split(discord_files, 1:length(discord_files))
  sample_names <- unlist(lapply(discord_files, function(x) {
    return(gsub(".discordant.*$", "", basename(x)))
  }))
  
  discord_list <- lapply(discord_files, function(x) {
    sample_id <- gsub(".discordant.*$", "", basename(x))
    # gunzip svaba discordant reads if needed and load table:
    if ( grepl("gz", x) ) {
      gunzip(x)
    }
    
    discord <- read.table(paste0(svaba_dir, "/", sample_id, ".discordant.txt"), 
                          sep = "\t", header = T, fill = T)
    
    # convert discord into GRanges object with both donor and acceptor data as main ranges:
    both_seqnames <- c(discord$chr1, discord$chr2)
    both_start <- c(discord$pos1, discord$pos2)
    both_actual_strand <- c(as.character(discord$strand1), as.character(discord$strand2))
    both_match_seqnames <- c(discord$chr2, discord$chr1)
    both_match_pos <- c(discord$pos2, discord$pos1)
    both_match_strand <- c(as.character(discord$strand2), as.character(discord$strand1))
    
    discord_gr <- GRanges(
      seqnames = Rle(paste0("chr",both_seqnames)),
      ranges = IRanges(start=as.integer(both_start), 
                       width=1),
      strand = Rle("*"),
      actual_strand = Rle(both_actual_strand),
      match_seqnames <- Rle(paste0("chr", both_match_seqnames)),
      match_pos <- as.integer(both_match_pos),
      match_strand <- both_match_strand
    )
    colnames(values(discord_gr)) <- c("actual_strand", "match_seqnames", 
                                      "match_pos", "match_strand")
    # remove unwanted chromosomes:
    discord_gr <- discord_gr[grep("M|G|K", as.character(seqnames(discord_gr)), 
                                  invert = T)]
    return(discord_gr)
    
  })
}
names(discord_list) <- sample_names

discord_list <- lapply(discord_list, function(x) {
  x$supporting_reads <- NA
  uniq <- unique(x)
  for ( j in 1:length(uniq) ) {
    print(j)
    uniq$supporting_reads[j] <- length(x[start(ranges(x))==start(ranges(uniq[j]))])
  }
  return(uniq)
})
saveRDS(discord_list, paste0(Robject_dir, "/discord_list_hg19.rds"))

verified_discords <- lapply(discord_list, function(x) return(x[x$supporting_reads > 2]))
saveRDS(verified_discords, paste0(Robject_dir, "/verified_discords_hg19.rds"))


##############################################################################
### 2. Find retrotransposon insertions ###
##############################################################################

if ( file.exists(paste0(Robject_dir, "/rep_insertions_hg19.rds"))) {
  
  rep_insertions <- readRDS(paste0(Robject_dir, "/rep_insertions_hg19.rds"))
  
} else {
  
  rep_insertions <- lapply(verified_discords, function(x) {
    # find overlaps of discordant reads with repeat annotation:
    olaps <- findOverlaps(x, rep_gr)
    
    incl_discord <- x[queryHits(olaps)]
    incl_rep <- rep_gr[subjectHits(olaps)]
    
    rep_insert_gr <- GRanges(
      seqnames = seqnames(incl_rep),
      ranges = ranges(incl_rep),
      strand = strand(incl_rep),
      gene_id = incl_rep$gene_id,
      fusion_half_1_seqnames = as.character(seqnames(incl_discord)),
      fusion_half_1_pos = as.integer(start(ranges(incl_discord))),
      fusion_half_1_strand = incl_discord$actual_strand,
      fusion_half_2_seqnames = incl_discord$match_seqnames,
      fusion_half_2_pos = incl_discord$match_pos,
      fusion_half_2_strand = incl_discord$match_strand
    )
    return(rep_insert_gr)
  })
  
  saveRDS(rep_insertions, paste0(Robject_dir, "/rep_insertions_hg19.rds"))
  
}

rep_inserts_pooled <- do.call(getMethod(c, "GenomicRanges"), rep_insertions)

saveRDS(rep_inserts_pooled, paste0(Robject_dir, "/repeat_insertions_hg19_pooled.rds"))

# make rep insertions gtf with co-ordinate of where repeat inserted:
rep_inserts_pooled_gtf <- GRanges(
  seqnames = Rle(rep_inserts_pooled$fusion_half_2_seqnames),
  ranges = IRanges(start=rep_inserts_pooled$fusion_half_2_pos, width=1),
  strand = "*",
  source = "svaba",
  feature = "repeat",
  score = ".",
  frame = ".",
  repeat_id = rep_inserts_pooled$gene_id
)
rep_inserts_pooled_gtf <- unique(rep_inserts_pooled_gtf)
export(rep_inserts_pooled_gtf, paste0(table_dir, "/repeat_inserts_pooled_hg19.gtf"))

rt_inserts_gtf <- rep_inserts_pooled_gtf[grep("^L1|Alu", rep_inserts_pooled_gtf$repeat_id)]
export(rt_inserts_gtf, paste0(table_dir, "/rt_inserts_pooled_hg19.gtf"))


##############################################################################
### 4. Load in CNVnator results and filter out CNVnator breakpoints not 
# supported by at least 4 additional svaba discordant read ###
##############################################################################

CNV_files <- split(
  grep(
    "subset", list.files(cnvnator_dir, pattern = "CNV_calls", 
                         full.names = T), value = T, invert = T
  ), 1:5
)
# CNV_files <- split(
#   CNV_files[gsub("_CNV_calls.txt", "", basename(CNV_files)) %in% samplenames], 1:length(CNV_files)
# )

# read in CNVnator results:
CNV_list <- lapply(CNV_files, function(x) {
  CNV_breakpoints <- read.table(x, header = F, sep = "\t")
  
  # format data frame:
  CNV_df <- CNV_breakpoints[grep("M|G|K", CNV_breakpoints$V2, invert = T),]
  CNV_df$type <- CNV_df$V1
  split_df <- do.call("rbind", strsplit(as.character(CNV_df$V2), ":"))
  CNV_df$seqnames <- split_df[,1]
  CNV_df$start <- do.call("rbind", strsplit(as.character(split_df[,2]), "-"))[,1]
  CNV_df$end <- do.call("rbind", strsplit(as.character(split_df[,2]), "-"))[,2]
  CNV_df <- subset(CNV_df, select = c("type", "seqnames", "start", "end"))
  
  # convert to GRanges object with all breakpoints as both start and end ranges:
  CNV_gr <- GRanges(
    seqnames = Rle(c(CNV_df$seqnames, CNV_df$seqnames)),
    ranges = IRanges(start=c(as.integer(CNV_df$start), as.integer(CNV_df$end)), 
                     end=c(as.integer(CNV_df$start), as.integer(CNV_df$end))),
    strand = Rle("*"),
    type = c(as.character(CNV_df$type), as.character(CNV_df$type)),
    actual_start = c(as.integer(CNV_df$start), as.integer(CNV_df$start)),
    actual_end = c(as.integer(CNV_df$end), as.integer(CNV_df$end))
  )
  return(CNV_gr)
})

i=1
verified_breakpoints <- lapply(CNV_list, function(x) {
  bp_olaps <- findOverlaps(x, discord_list[[i]])
  i <<- i+1
  return(unique(x[queryHits(bp_olaps)]))
})
pooled_verified_breakpoints <- do.call(getMethod(c, "GenomicRanges"), verified_breakpoints)
  

# create deletions/duplications gtf file:
pv_deletions <- pooled_verified_breakpoints[pooled_verified_breakpoints$type == "deletion"]
pooled_deletion_gtf <- GRanges(
  seqnames = seqnames(pv_deletions),
  ranges = ranges(pv_deletions),
  strand = "*",
  source = "cnvnator",
  feature = "breakpoint",
  score = ".",
  frame = ".",
  type = pv_deletions$type
)
export(pooled_deletion_gtf, paste0(table_dir, "/deletion_breakpoints_pooled.gtf"))

pv_duplications <- pooled_verified_breakpoints[pooled_verified_breakpoints$type == "duplication"]
pooled_duplication_gtf <- GRanges(
  seqnames = seqnames(pv_duplications),
  ranges = ranges(pv_duplications),
  strand = "*",
  source = "cnvnator",
  feature = "breakpoint",
  score = ".",
  frame = ".",
  type = pv_duplications$type
)
export(pooled_duplication_gtf, paste0(table_dir, "/duplication_breakpoints_pooled.gtf"))


##############################################################################
### 5. Find odds of any repeat insertions at varying distances from verified 
# breakpoints ###
##############################################################################

# convert to list with first element deletions, second duplications:
deletion_gr <- pooled_verified_breakpoints[pooled_verified_breakpoints$type == "deletion",]
duplication_gr <- pooled_verified_breakpoints[pooled_verified_breakpoints$type == "duplication",]

# create function to expand genomic ranges either direction from start/end
# and optionally remove original ranges:
n=1
exp_gr <- function(Length, gr, orig) {
  
  print(paste0("Processing range ", n))
  
  # define lengths of chromosomes:
  seq_lengths <- seqlengths(Hsapiens)[!grepl("_",
                                             names(seqlengths(Hsapiens)))]
  # remove unwanted chromosome constructs:
  gr <- gr[grep("[0-9].[0-9]|MT|G|K", seqnames(gr), 
                      invert = T)]
  
  if ( length(gr) > 0 ) {
    # reduce ranges to ensure no overlaps/double annotations:
    gr <- unique(reduce(gr))

    # assign length of all chromosomes to the variable name 'seq_lengths':
    gr$seq_lengths <- rep(NA, length(gr))
    
    for ( v in names(seq_lengths) ) {
      gr$seq_lengths[as.character(seqnames(gr)) == v] <- seq_lengths[v]
    }
    
    # add length upstream of each range to another gr:
    start_gr <- gr
    # make end position equal to original start position:
    end(ranges(start_gr)) <- start(ranges(start_gr))
    # make start position the original start position - Length if the original start is at least Length
    # away from start of chromosome:
    start(start_gr)[start(ranges(start_gr)) >= Length] <- 
      start(ranges(start_gr))[start(ranges(start_gr)) >= Length] - Length
    
    # add length downstream of each range to another gr:
    end_gr <- gr
    
    # make start position equal to original end position:
    start(ranges(end_gr)) <- end(ranges(end_gr))
    # make end position the original end position + Length if the end is at least Length
    # away from end of chromosome:
    far_enuff <- end(ranges(end_gr)) <= (end_gr$seq_lengths - Length)
    end(ranges(end_gr))[far_enuff] <- end(ranges(end_gr))[far_enuff] + Length
    
    if (orig) {
      # include original sequences in GRanges object:
      start(ranges(gr)) <- start(ranges(start_gr))
      end(ranges(gr)) <- end(ranges(end_gr))
      # remove ranges with end=1:
      gr <- gr[!(end(ranges(gr))==1)]
      n <<- n+1
      return(gr)
    }
    
    # remove co-ordinates that fall after the end of seq_lengths values:
    end_gr <- end_gr[!(end(ranges(end_gr)) > end_gr$seq_lengths)]
    
    # make end position the end of the chromosome if the original end is not at least Length
    # away from end of chromosome:
    not_far_enuff <- end(ranges(end_gr)) > (end_gr$seq_lengths - Length)
    if ( length(not_far_enuff) > 0 ) {
      end(ranges(end_gr))[not_far_enuff] <- end_gr$seq_lengths[not_far_enuff]
    }
    gr <- c(start_gr, end_gr)
    # remove ranges with end=1:
    gr <- gr[!(end(ranges(gr))==1)]
    # return combined start and end annotations:
    n <<- n+1
    return(gr)
  }
}

# Create list of deletions with ranges extending different distances outwards 
# from segments, without original (deleted) segment:
n=1
deletion_borders <- lapply(dist, exp_gr, gr=deletion_gr, orig=FALSE)
non_deletion <- lapply(deletion_borders, gaps)

# Create list of duplications with ranges extending different distances outwards 
# from segments, keeping original segment:
n=1
duplication_borders <- lapply(dist, exp_gr, gr=duplication_gr, orig=FALSE)
non_duplication <- lapply(duplication_borders, gaps)


# create function to calculate odds of breakpoint region enrichment for 
# repeat insertions compared with non-breakpoint regions:
k=1
oddsRatio <- function(region, non_region, annot){
  
  # determine how many repeats in region:
  r_olaps <- findOverlaps(annot, region)
  if ( length(r_olaps) > 0 ) {
    region_no <- length(unique(queryHits(r_olaps)))
    # determine length of region:
    region_length <- sum(width(reduce(region)))
    
    # determine how many repeats in non-region:
    nr_olaps <- findOverlaps(annot, non_region[[k]])
    non_region_no <- length(unique(queryHits(nr_olaps)))
    # determine length of non-region:
    non_region_length <- sum(width(reduce(non_region[[k]])))
    
    # calculate odds of region containing a repeat insertion:
    first_ratio <- region_no/region_length
    
    # calculate odds of non-region containing a repeat insertion:
    second_ratio <- non_region_no/non_region_length
    
    # calculate odds ratio:
    odds <- first_ratio/second_ratio
    
    # calculate standard error for odds:
    se <- sqrt(1/region_length + 1/region_no + 1/non_region_length + 
                 1/non_region_no)
    print("Odds are:")
    cat(odds)
    cat("\n")
    print("Std error is:")
    cat(se)
    cat("\n")
    
    # calculate p-value for odds:
    database_enrichment <-
      matrix( c(as.numeric(region_no), as.numeric(non_region_no), 
                as.numeric(region_length), as.numeric(non_region_length)), nrow = 2,
              ncol = 2 )
    colnames(database_enrichment) = c("insertions", "length")
    rownames(database_enrichment) = c("region", "non_region")
    
    pval <- chisq.test(database_enrichment)$p.value
    print("P-value is:")
    cat(pval)
    cat("\n")
    
    result <- data.frame(odds, se, pval)
    if ( result[1,3] == 0 ) {
      result[1,3] <- "< 0.1e15"
    }
    colnames(result) <- c("odds", "std_error", "pval")
    
    k <<- k+1
    return(result)
  } else {
    repeat_id <- unique(annot$repeat_id)
    print(paste0("No ", repeat_id, " in region"))
    return(NA)
  }
}

k=1
deletion_repeat_enrichment_odds <- lapply(deletion_borders, oddsRatio, 
  non_deletion, rt_inserts_gtf)
names(deletion_repeat_enrichment_odds) <- as.integer(unlist(dist))

k=1
duplication_repeat_enrichment_odds <- lapply(duplication_borders, oddsRatio, 
  non_duplication, rt_inserts_gtf)
names(duplication_repeat_enrichment_odds) <- as.integer(unlist(dist))

# plot odds for each length of breakpoint regions:
deletion_odds <- do.call("rbind", deletion_repeat_enrichment_odds)
deletion_odds$distance_interval_from_breakpoint <- 
  factor(rownames(deletion_odds), levels = as.integer(unlist(dist)))

deletion_log_odds <- deletion_odds
deletion_log_odds$log_odds <- log(deletion_log_odds$odds)

delcols <- rev(brewer.pal(6, "Greens"))

p <- ggplot(data = deletion_log_odds, aes(x=distance_interval_from_breakpoint, 
  y=log_odds))
p <- p + geom_bar(aes(fill = distance_interval_from_breakpoint), 
stat = "identity")
p <- p + scale_fill_manual(values=delcols)
p <- p + geom_errorbar(aes(ymin=log_odds-std_error, ymax=log_odds+std_error), 
  width=.2, position=position_dodge(.9))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.55), 
               axis.text=element_text(size=14), axis.title = element_text(size=14))
p <- p + ylab("Log odds")
p <- p + xlab("Distance interval from breakpoint")
p <- p + guides(fill=FALSE)
pdf(paste0(plot_dir, "/repeat_enrichment_at_different_intervals_from_deletions.pdf"))
p
dev.off()

dupcols <- rev(brewer.pal(6, "Purples"))

duplication_odds <- do.call("rbind", duplication_repeat_enrichment_odds)
duplication_odds$distance_interval_from_breakpoint <- 
  factor(rownames(duplication_odds), levels = as.integer(unlist(dist)))

duplication_log_odds <- duplication_odds
duplication_log_odds$log_odds <- log(duplication_log_odds$odds)

p <- ggplot(data = duplication_log_odds, aes(x=distance_interval_from_breakpoint, 
                                          y=log_odds))
p <- p + geom_bar(aes(fill = distance_interval_from_breakpoint), 
                  stat = "identity")
p <- p + scale_fill_manual(values=dupcols)
p <- p + geom_errorbar(aes(ymin=log_odds-std_error, ymax=log_odds+std_error), 
                       width=.2, position=position_dodge(.9))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.55), 
  axis.text=element_text(size=14), axis.title = element_text(size=14))
p <- p + ylab("Log odds")
p <- p + xlab("Distance interval from breakpoint")
p <- p + guides(fill=FALSE)
pdf(paste0(plot_dir, "/repeat_enrichment_at_different_intervals_from_duplications.pdf"))
p
dev.off()


##############################################################################
### 6. Find odds of specific repeat insertions at varying distances from 
# verified breakpoints ###
##############################################################################

# create function to remove NA from lists:
na.omit.list <- function(y) { return(y[!sapply(y, function(x) all(is.na(x)))]) }

# split DE repeat insert GRanges by repeat ID:
DE_inserts_gtf <- rt_inserts_gtf[rt_inserts_gtf$repeat_id %in% DE]
DE_inserts_split <- split(DE_inserts_gtf, DE_inserts_gtf$repeat_id)

DE_rt_odds_dup <- lapply(DE_inserts_split, function(y) {
  print(unique(y$repeat_id))
  res <- lapply(duplication_borders, oddsRatio, non_duplication, y)
  names(res) <- as.character(as.integer(unlist(dist)))
  # remove NAs:
  res <- na.omit.list(res)
  # remove empty elements:
  res <- res[lapply(res, length) > 0]
  # remove elements with p > 0.05:
  res <- res[unlist(lapply(res, function(a) a$pval < 0.05))]
  res <- res[unlist(lapply(res, function(a) a$odds != "Inf"))]
  
  k <<- 1
  return(res)
})
DE_rt_odds_dup <- DE_rt_odds_dup[lapply(DE_rt_odds_dup, length) > 0]

dup_plots <- lapply(as.character(as.integer(dist)), function(y) {
  dup <- lapply(DE_rt_odds_dup, function(x) {
    return(x[grep(paste0(y, "$"), names(x))])
  })
  dup <- dup[lapply(dup, length) > 0]
  if ( length(dup) > 0 ) {
    for ( d in 1:length(dup) ) {
      if (d==1) {
        res <- data.frame(dup[[d]])
      } else {
        res <- rbind(res, unlist(dup[[d]]))
      }
    }
    rownames(res) <- names(dup)
    colnames(res) <- c("odds", "std_error", "pval")
    return(res)
  }
})
names(dup_plots) <- as.character(as.integer(dist))

######
dup_plot <- do.call("rbind", dup_plots)
dup_plot$distance <- gsub("\\..*$", "", rownames(dup_plot))
dup_plot <- dup_plot[c(1, 4:6, 8, 11:12),]
rownames(dup_plot)[1] <- "L1M4b"
rownames(dup_plot) <- gsub("^.*\\.", "", rownames(dup_plot))

dup_log_plot <- dup_plot
dup_log_plot$log_odds <- log(as.numeric(dup_log_plot$odds))
dup_log_plot$std_error <- as.numeric(dup_log_plot$std_error)
dup_log_plot$retrotransposon <- rownames(dup_log_plot)

dup_log_plot$retrotransposon <- factor(dup_log_plot$retrotransposon, levels=dup_log_plot$retrotransposon)

p <- ggplot(data = dup_log_plot, aes(x=retrotransposon,
                                     y=log_odds))
p <- p + geom_bar(aes(fill = distance),
                  stat = "identity")
p <- p + scale_fill_manual(values=dupcols)
p <- p + geom_errorbar(aes(ymin=log_odds-std_error, ymax=log_odds+std_error),
                       width=.2, position=position_dodge(.9))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.55),
               axis.text=element_text(size=14), axis.title = element_text(size=14), axis.title.x=element_blank())
p <- p + ylab("Log odds")
#p <- p + xlab("Retrotransposon")
#p <- p + guides(fill=FALSE)
pdf(paste0(plot_dir, "/DE_repeat_enrichment_at_duplications.pdf"))
print(p)
dev.off()

#####

DE_rt_odds_del <- lapply(DE_inserts_split, function(y) {
  print(unique(y$repeat_id))
  res <- lapply(deletion_borders, oddsRatio, non_deletion, y)
  names(res) <- as.character(as.integer(unlist(dist)))
  # remove NAs:
  res <- na.omit.list(res)
  # remove empty elements:
  res <- res[lapply(res, length) > 0]
  # remove elements with p > 0.05:
  res <- res[unlist(lapply(res, function(a) a$pval < 0.05))]
  res <- res[unlist(lapply(res, function(a) a$odds != "Inf"))]
  
  k <<- 1
  return(res)
})
DE_rt_odds_del <- DE_rt_odds_del[lapply(DE_rt_odds_del, length) > 0]

del_plots <- lapply(as.character(as.integer(dist)), function(y) {
  del <- lapply(DE_rt_odds_del, function(x) {
    return(x[grep(paste0(y, "$"), names(x))])
  })
  del <- del[lapply(del, length) > 0]
  if ( length(del) > 0 ) {
    for ( d in 1:length(del) ) {
      if (d==1) {
        res <- data.frame(del[[d]])
      } else {
        res <- rbind(res, unlist(del[[d]]))
      }
    }
    rownames(res) <- names(del)
    colnames(res) <- c("odds", "std_error", "pval")
    return(res)
  }
})
names(del_plots) <- as.character(as.integer(dist))

######

del_plot <- do.call("rbind", del_plots)
del_plot$distance <- gsub("\\..*$", "", rownames(del_plot))
del_plot <- del_plot[c(1:3),]
rownames(del_plot) <- gsub("^.*\\.", "", rownames(del_plot))

del_log_plot <- del_plot
del_log_plot$log_odds <- log(as.numeric(del_log_plot$odds))
del_log_plot$std_error <- as.numeric(del_log_plot$std_error)
del_log_plot$retrotransposon <- rownames(del_log_plot)

del_log_plot$retrotransposon <- factor(del_log_plot$retrotransposon, levels=del_log_plot$retrotransposon)

p <- ggplot(data = del_log_plot, aes(x=retrotransposon,
                                     y=log_odds))
p <- p + geom_bar(aes(fill = distance),
                  stat = "identity")
p <- p + scale_fill_manual(values=delcols)
p <- p + geom_errorbar(aes(ymin=log_odds-std_error, ymax=log_odds+std_error),
                       width=.2, position=position_dodge(.9))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.55),
  axis.text=element_text(size=14), axis.title = element_text(size=14), axis.title.x=element_blank())
p <- p + ylab("Log odds")
#p <- p + xlab("Retrotransposon")
#p <- p + guides(fill=FALSE)
pdf(paste0(plot_dir, "/DE_repeat_enrichment_at_deletions.pdf"))
print(p)
dev.off()

######

del_rt_gtf <- rep_inserts_pooled_gtf[rep_inserts_pooled_gtf$repeat_id %in% del_log_plot$retrotransposon]
export(del_rt_gtf, paste0(table_dir, "/rt_insertions_at_deletion_breakpoints.gtf"))

dup_rt_gtf <- rep_inserts_pooled_gtf[rep_inserts_pooled_gtf$repeat_id %in% dup_log_plot$retrotransposon]
export(dup_rt_gtf, paste0(table_dir, "/rt_insertions_at_duplication_breakpoints.gtf"))

# for ( n in 1:length(dup_plots) ) {
#   if ( !is.null(dup_plots[[n]]) ) {
#     dup_log_plot <- dup_plots[[n]]
#     dup_log_plot$log_odds <- log(as.numeric(dup_log_plot$odds))
#     dup_log_plot$std_error <- as.numeric(dup_log_plot$std_error)
#     dup_log_plot$retrotransposon <- rownames(dup_plots[[n]])
#     
#     p <- ggplot(data = dup_log_plot, aes(x=retrotransposon, 
#                                          y=log_odds))
#     p <- p + geom_bar(aes(fill = retrotransposon), 
#                       stat = "identity")
#     #p <- p + scale_fill_manual(values=duprtcols)
#     p <- p + geom_errorbar(aes(ymin=log_odds-std_error, ymax=log_odds+std_error), 
#                            width=.2, position=position_dodge(.9))
#     p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.55), 
#                    axis.text=element_text(size=14), axis.title = element_text(size=14))
#     p <- p + ylab("Log odds")
#     p <- p + xlab("Retrotransposon")
#     p <- p + guides(fill=FALSE)
#     pdf(paste0(plot_dir, "/repeat_enrichment_", names(dup_plots)[n], "_from_duplications.pdf"))
#     print(p)
#     dev.off()
#   }
# }






sig_rep_odds_del <- lapply(rep_inserts_split, function(y) {
  print(unique(y$repeat_id))
  res <- lapply(deletion_borders, oddsRatio, non_deletion, y)
  names(res) <- as.character(as.integer(unlist(dist)))
  # remove NAs:
  res <- na.omit.list(res)
  # remove empty elements:
  res <- res[lapply(res, length) > 0]
  # remove elements with p > 0.05:
  res <- res[unlist(lapply(res, function(a) a$pval < 0.05))]
  res <- res[unlist(lapply(res, function(a) a$odds != "Inf"))]
  
  k <<- 1
  return(res)
})
sig_rep_odds_del <- sig_rep_odds_del[lapply(sig_rep_odds_del, length) > 0]

odds_del_1000 <- lapply(sig_rep_odds_del, function(x) return(x[[1]]))
odds_del_2000 <- lapply(sig_rep_odds_del, function(x) return(x[[2]]))
odds_del_5000 <- lapply(sig_rep_odds_del, function(x) return(x[[3]]))
odds_del_10000 <- lapply(sig_rep_odds_del, function(x) return(x[[4]]))
odds_del_100000 <- lapply(sig_rep_odds_del, function(x) return(x[[5]]))

save.image(paste0(Robject_dir, "/individual_repeat_enrichment_odds_calculated.RData"))

# p <- ggplot(data = duplication_log_odds, aes(x=distance_interval_from_breakpoint, 
#                                              y=log_odds))
# p <- p + geom_bar(aes(fill = distance_interval_from_breakpoint), 
#                   stat = "identity")
# p <- p + scale_fill_manual(values=dupcols)
# p <- p + geom_errorbar(aes(ymin=log_odds-std_error, ymax=log_odds+std_error), 
#                        width=.2, position=position_dodge(.9))
# p <- p + theme(axis.text.x = element_text(angle = 90, vjust=0.55), 
#                axis.text=element_text(size=14), axis.title = element_text(size=14))
# p <- p + ylab("Log odds")
# p <- p + xlab("Distance interval from breakpoint")
# p <- p + guides(fill=FALSE)
# pdf(paste0(plot_dir, "/repeat_enrichment_at_different_intervals_from_duplications.pdf"))
# p
# dev.off()


