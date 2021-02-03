return_top_DE <- function(
  in_dir,
  up_filename,
  down_filename,
  max_pval = 0.05,
  min_logfc = 0.7
) {

  # load DE retrotransposons:
  up_DE <- read.table(
    paste0(in_dir, up_filename),
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )
  sig_up <- up_DE[up_DE$FDR < max_pval,]
  sig_up <- up_DE[up_DE$logFC >= min_logfc,]
  sig_up <- sig_up[order(sig_up$logFC, decreasing = T),]
  
  down_DE <- read.table(
    paste0(in_dir, down_filename),
    header = T,
    sep = "\t",
    stringsAsFactors = F
  )
  sig_down <- down_DE[down_DE$FDR < max_pval,]
  sig_down <- down_DE[down_DE$logFC <= -min_logfc,]
  sig_down <- sig_down[order(sig_down$logFC),]
  
  # take the top 10 up and 10 down retrotransposons:
  if (nrow(sig_up) > 9) {
    top_up <- sig_up[1:10,]
  } else {
    top_up <- sig_up
  }
  
  if (nrow(sig_down) > 9) {
    top_down <- sig_down[1:10,]
  } else {
    top_down <- sig_down
  }

  return(
    c(rownames(top_up), rownames(top_down))
  )

}