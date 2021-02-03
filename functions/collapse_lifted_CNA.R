collapse_lifted_CNA <- function(lifted_gr) {
  
  if (length(lifted_gr) > 0) {

    main_chr <- seqnames(lifted_gr)@values[
      which.max(seqnames(lifted_gr)@lengths)
    ]
    main_gr <- sort(lifted_gr[seqnames(lifted_gr) == main_chr])
    new_gr <- GRanges(
      seqnames = seqnames(main_gr)[1],
      ranges = IRanges(
        start = start(main_gr[1]), 
        end = end(main_gr[length(main_gr)])
      ),
      strand = Rle("*")
    )

    return(new_gr)

  } else {
    return(NULL)
  }

}