#' Call broad peak by Hidden Markov Models with t emissions
#' 
#' Call very broad peaks for data such as LAD domains, NAD domains.
#' Reads will be count by each bins. Only bins with at least given reads
#' (defined by background parameter) for all samples (pool all reads for 
#' each bin) will be subsequently normalized. These bins will be first 
#' normalized to CPM (count per million) reads and then do log2 transform
#' for the ratio over control with a pseudocount. 
#' The peaks were defined by running a hidden markov model over the 
#' normalized values (using the R-package HMMt).
#' 
#' @param treatment,control Bam file of treatments and controls. Make sure the
#' index file keep same prefix name with bam file.
#' @param binSize The size of bins for count
#' @param background Only bins with at least background reads (pool all reads)
#' will be subsequently normalized.
#' @param pseudocount default 1.
#' @param gapwidth The Ranges of peaks separated by a gap less than gapwith 
#' positions will be merged.
#' @param ... parameters passed to \link[HMMt]{BaumWelchT} 
#' except m (fixed to 2).
#' @return a list with elements counts and peaks. Bothe counts and peaks 
#' are GRanges objects.
#' @export
#' @import GenomicRanges
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom SummarizedExperiment assays rowRanges
#' @importFrom S4Vectors DataFrame
#' @importFrom HMMt BaumWelchT
#' @examples 
#' treatment <- system.file("extdata", "LB1.KD.chr1_1_5000000.bam",
#'                           package = "HMMtBroadPeak",
#'                           mustWork = TRUE)
#' control   <- system.file("extdata", "LB1.WT.chr1_1_5000000.bam",
#'                           package = "HMMtBroadPeak",
#'                           mustWork = TRUE)
#' called <- HMMtBroadPeak(treatment, control)
#' called$peaks
#' plotPeaks(called, seqname="chr1")
#' 
HMMtBroadPeak <- function(treatment, control, binSize=5e3,
                          background=10, pseudocount=1, 
                          gapwidth=binSize,
                          ...){
  if(missing(treatment) || missing(control)){
    stop("treatment/control is required")
  }
  args <- list(...)
  if("m" %in% names(args)){
    stop("The number of states is fixed to 2. Please remove the m parameter.")
  }
  ## check header are identical
  seqlen <- check_identical_header(c(treatment, control))
  ## create tileGenome
  rr <- tileGenome(seqlengths = seqlen, tilewidth = binSize)
  rr <- unlist(rr)
  suppressWarnings({
    reads_treatment <- summarizeOverlaps(features = rr,
                                       reads = treatment,
                                       mode = count_by_overlap,
                                       ignore.strand=TRUE)
    reads_control <- summarizeOverlaps(features = rr,
                                       reads = control,
                                       mode = count_by_overlap,
                                       ignore.strand=TRUE)
  })
  if(!identical(rr, rowRanges(reads_treatment)) || 
     !identical(rr, rowRanges(reads_control))){
    stop("unexpect event happened.")
  }
  counts_treatment <- rowSums(assays(reads_treatment)$counts)
  counts_control <- rowSums(assays(reads_control)$counts)
  mcols(rr) <- DataFrame(treatment = counts_treatment,
                         control = counts_control)
  rr$filter <- rowSums(as.data.frame(mcols(rr)))>=background
  rr$treatment <- log2((rr$treatment*1e6/sum(rr$treatment))+pseudocount)
  rr$control <- log2((rr$control*1e6/sum(rr$control))+pseudocount)
  rr$log2ratio <- rr$treatment - rr$control
  rr$log2ratio[!rr$filter] <- 0
  rr <- split(rr, seqnames(rr))
  hmmt <- lapply(rr, function(.ele){
    BaumWelchT(.ele$log2ratio, m=2, ...)
  })
  rr <- mapply(function(.rr, .hmmt){
    .rr$ViterbiPath <- .hmmt@ViterbiPath
    .rr
  }, rr, hmmt)
  rr <- unlist(GRangesList(rr))
  out <- list(counts=rr)
  ## rr$ViterbiPath==2 is peak
  rr <- rr[rr$ViterbiPath==2]
  rr <- reduce(rr, min.gapwidth = gapwidth)
  out$peaks <- rr
  return(out)
}

#' helper function to check the bam file header
#' @param bams Bam file of treatments and controls
#' @return seqlengths
#' @importFrom Rsamtools scanBamHeader
#' @import GenomicRanges
#' @import IRanges
check_identical_header <- function(bams){
  if(missing(bams)){
    stop("bams is required.")
  }
  headers <- lapply(bams, scanBamHeader)
  headers_SQ <- lapply(headers, function(.ele){
    .ele[[1]]$targets
  })
  SQ <- headers_SQ[[1]]
  if(length(bams)>1){
    ident <- TRUE
    for(i in seq_along(headers_SQ)[-1]){
      if(!identical(SQ, headers_SQ[[i]])){
        ident <- FALSE
        warning("Assemblies in bam files are not identical.")
        break()
      }
    }
    if(!ident){
      warning("Assemblies are not identical, use the first bam file header")
      SQ_name <- lapply(headers_SQ, names)
      SQ_name <- Reduce(intersect, SQ_name)
      SQ <- SQ[SQ_name]
    }
  }
  SQ
}

count_by_overlap <- function(features, reads,  ignore.strand, inter.feature) { 
  ## perform filtering, or subsetting etc. 
  countOverlaps(features, reads, minoverlap = 1L,
                ignore.strand = ignore.strand)
}
