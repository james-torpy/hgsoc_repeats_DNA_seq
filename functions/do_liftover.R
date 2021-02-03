do_liftover <- function(gr, original, output) {
  
  library(rtracklayer)

  output <- gsub("h", "H", output)
  chain <- import.chain(paste0(ref_dir, "liftover/", original, 
    "To", output, ".over.chain"))

  return(liftOver(gr, chain))

}