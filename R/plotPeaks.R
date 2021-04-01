#' Plot tool for called peaks
#' 
#' Plot the log2 transformed signal along a given chromosome with the peaks.
#' 
#' @param res output of \link{HMMtBroadPeak}.
#' @param seqname the chromosome to plot.
#' @return an ggplot object
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom ggplot2 ggplot aes_string geom_line
#' @importFrom S4Vectors runValue runLength
#' @export
#' 
plotPeaks <- function(res, seqname){
  stopifnot("The res must be output of HMMtBroadPeak function." =
              is.list(res) && all(c("counts", "peaks") %in% names(res)))
  stopifnot(is.character(seqname))
  stopifnot(length(seqname)==1)
  stopifnot("seqname is not shown in counts table." = 
              seqname %in% seqlevels(res$counts))
  cvg <- coverage(res$counts, weight = res$counts$log2signal)
  peak <- coverage(res$peaks)
  cvg <- cvg[[seqname]]
  peak <- peak[[seqname]]
  plotdata <- 
    rbind(data.frame(
            x=cumsum(c(1, runLength(cvg)))[-(length(runLength(cvg))+1)],
            value = runValue(cvg),
            data = "log2signal"),
          data.frame(x=cumsum(runLength(cvg)),
                     value = runValue(cvg),
                     data = "log2signal"),
          data.frame(
            x=cumsum(c(1, runLength(peak)))[-(length(runLength(peak))+1)],
            value = runValue(peak)*
              quantile(x=cvg, probs=.8, na.rm = TRUE, names = FALSE),
            data = "called peak"),
          data.frame(x=cumsum(runLength(peak)),
                     value = runValue(peak)*
                       quantile(x=cvg, probs=.8, na.rm = TRUE, names = FALSE),
                     data = "called peak"))
  ggplot(data = plotdata,
         mapping = aes_string(x="x", y="value", color = "data")) +
    geom_line()
}