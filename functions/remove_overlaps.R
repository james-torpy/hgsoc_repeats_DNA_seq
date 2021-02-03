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