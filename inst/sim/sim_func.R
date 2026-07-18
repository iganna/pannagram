# All function to find similarities

# ---- PATCHED: assemble-then-threshold similarity search ----
# Union length (bp) covered by inclusive [a,b] intervals (merged).
mergeLen <- function(b, e){
  if(length(b) == 0) return(0)
  iv <- cbind(pmin(b, e), pmax(b, e))
  iv <- iv[order(iv[,1]), , drop = FALSE]
  tot <- 0; cs <- iv[1,1]; ce <- iv[1,2]
  if(nrow(iv) > 1){
    for(i in 2:nrow(iv)){
      if(iv[i,1] <= ce + 1){ ce <- max(ce, iv[i,2]) }
      else { tot <- tot + (ce - cs + 1); cs <- iv[i,1]; ce <- iv[i,2] }
    }
  }
  tot + (ce - cs + 1)
}

# Merged inclusive [a,b] intervals encoded as "a-b;a-b;..." (genome sub-segments
# of a copy; gaps = insertions/nesting, excluded from the painted bp).
mergeIvStr <- function(b, e){
  if(length(b) == 0) return("")
  iv <- cbind(pmin(b, e), pmax(b, e))
  iv <- iv[order(iv[,1]), , drop = FALSE]
  segs <- character(0); cs <- iv[1,1]; ce <- iv[1,2]
  if(nrow(iv) > 1){
    for(i in 2:nrow(iv)){
      if(iv[i,1] <= ce + 1){ ce <- max(ce, iv[i,2]) }
      else { segs <- c(segs, paste0(cs, "-", ce)); cs <- iv[i,1]; ce <- iv[i,2] }
    }
  }
  paste(c(segs, paste0(cs, "-", ce)), collapse = ";")
}

#' Find Hits in Reference (patched, assemble-then-threshold)
#'
#' Groups same-(query, chromosome, strand) BLAST HSPs into colinear copies by
#' gap-clustering along the genome, then applies the similarity / coverage
#' thresholds ONCE to each assembled copy:
#'   - consensus coverage = merged query-side bp / query length  >= coverage
#'   - identity = alignment-length-weighted mean pident          >= sim.cutoff
#'   - genome footprint    = ref span / query length             >= coverage
#'   - symmetric (optional): query length / ref span             >= coverage
#'     (when FALSE, copies whose genomic footprint is LONGER than the consensus
#'      -- i.e. carry an insertion / nested element -- are still kept.)
#'
#' This replaces the per-HSP "exact match" + fragile fragment-stitching logic so
#' that single sub-coverage HSPs get partial credit and inserted/nested copies
#' are not rejected. No pre-filtering of HSPs by pident before assembly.
findHitsInRef <- function(v, sim.cutoff, coverage = NULL, symmetric = FALSE,
                          gap.factor = 1.0, gap.abs = 20000,
                          min.new.frac = 0.5, echo = FALSE){
  if(is.null(coverage)) coverage <- sim.cutoff
  if(!('len1' %in% colnames(v))) stop('No column len1 in the data.frame')
  if(!((sim.cutoff >= 0) & (sim.cutoff <= 1)))
    stop(paste("Similarity cutoff should be between 0 and 1, now", sim.cutoff))
  if(!((coverage >= 0) & (coverage <= 1)))
    stop(paste("Coverage cutoff should be between 0 and 1, now", coverage))

  empty <- data.frame(V1=character(), V2=numeric(), V3=numeric(),
                      V4=numeric(), V5=numeric(), V6=numeric(), V7=numeric(),
                      V8=character(), len1=numeric(), strand=character(),
                      subiv=character(),
                      stringsAsFactors = FALSE)
  if(nrow(v) == 0) return(empty)

  gmin  <- pmin(v$V4, v$V5)
  gmax  <- pmax(v$V4, v$V5)
  strnd <- ifelse(v$V4 > v$V5, '-', '+')
  alen  <- v$V7   # alignment length (kept as col 7 from BLAST) -> identity weight
  pid   <- v$V6
  # PATCHED: RM-style substitution-only identity when the BLAST table carries
  # the optional mismatch/gaps columns (mism = V11, gaps = V12). RepeatMasker's
  # %div counts substitutions only and reports indels separately (del%/ins%),
  # whereas blastn `pident` charges gap columns against identity -- so for
  # insertion-rich diverged copies blastn pident sits ~5-7 pts below RM identity
  # and trips the sim gate. We reconcile by excluding indel columns:
  #   identity = identical_bases / aligned_non-gap_columns.
  has.sub <- all(c('mism', 'gaps') %in% colnames(v))
  if(has.sub){
    matches <- v$V7 - v$mism - v$gaps   # identical aligned bases
    nongap  <- v$V7 - v$gaps            # aligned columns excluding indels (= matches + mism)
  }

  key <- paste(v$V1, v$V8, strnd, sep = '\r')
  # PATCHED (speed): group HSPs by key in a single O(N) pass via split() instead
  # of `for(k in unique(key)) which(key==k)` which was O(n_keys * n_hsps) and
  # dominated runtime (~23s -> ~1s). Output is identical (final res is sorted
  # downstream; copy IDs are arbitrary labels).
  groups <- split(seq_along(key), key)
  res.list <- vector("list", length(groups))
  ng <- 0L
  for(idx in groups){
    o   <- idx[order(gmin[idx])]
    len1 <- v$len1[o[1]]
    # Genomic gap threshold for chaining: relative to the consensus length but
    # capped at an absolute size (a nested insertion is at most a few kb-tens of
    # kb; without the cap, gap.factor * len1 becomes hundreds of kb for long
    # consensi and merges unrelated copies).
    gthr <- min(max(1, gap.factor * len1), gap.abs)
    n.o <- length(o)
    if(n.o == 1){
      cl <- 1L
    } else {
      # Collinear chaining (vectorised, O(n)). A new copy starts when the
      # genomic gap exceeds gthr. Additionally, for LONG-RANGE joins (gap beyond
      # one consensus length -- the part that extends stock behaviour) the next
      # HSP must add NEW consensus: an HSP that mostly re-covers the previous
      # one's consensus span is a SEPARATE copy of the same family, and merging
      # it would create a chimera. The gate is not applied to small gaps, so
      # stock chaining (and its handling of overlapping HSPs of one copy) is
      # preserved -- the fix only ADDS safe long-range joins.
      gapd <- gmin[o][-1] - gmax[o][-n.o]
      qs <- v$V2[o]; qe <- v$V3[o]
      ov <- pmax(0, pmin(qe[-1], qe[-n.o]) - pmax(qs[-1], qs[-n.o]) + 1)
      new.frac <- 1 - ov / (qe[-1] - qs[-1] + 1)
      brk <- (gapd > gthr) | ((gapd > len1) & (new.frac < min.new.frac))
      cl <- cumsum(c(1, brk * 1))
    }
    for(ci in unique(cl)){
      g <- o[cl == ci]
      qcov  <- mergeLen(v$V2[g], v$V3[g])
      if(has.sub){
        # RM-style: substitution-only identity (indel columns excluded).
        ident <- 100 * sum(matches[g]) / sum(nongap[g])
      } else {
        # fallback: alignment-length-weighted mean of gap-inclusive blast pident.
        ident <- sum(pid[g] * alen[g]) / sum(alen[g])
      }
      gspan <- max(gmax[g]) - min(gmin[g]) + 1
      ok <- (qcov / len1 >= coverage) &
            (ident >= sim.cutoff * 100) &
            (gspan / len1 >= coverage)
      if(symmetric) ok <- ok & (len1 / gspan >= coverage)
      if(isTRUE(ok)){
        ng <- ng + 1L
        res.list[[ng]] <- data.frame(
          V1 = v$V1[g[1]], V2 = min(v$V2[g]), V3 = max(v$V3[g]),
          V4 = min(gmin[g]), V5 = max(gmax[g]), V6 = ident, V7 = qcov,
          V8 = v$V8[g[1]], len1 = len1, strand = strnd[g[1]],
          subiv = mergeIvStr(gmin[g], gmax[g]),
          stringsAsFactors = FALSE)
      }
    }
  }
  if(ng == 0) return(empty)
  do.call(rbind, res.list[seq_len(ng)])
}


#' Convert BLAST results to GFF format
#'
#' @description
#' `blastres2gff` converts a data frame containing BLAST results into a GFF formatted file.
#'
#' @param v.blast A data frame containing the BLAST results. Expected columns are V1, V4, V5, V7, V8, len1, and strand.
#' @param gff.file The file path where the GFF output will be saved.
#' @param to.sort Boolean value indicating whether the output should be sorted. If `TRUE` (default), 
#' the output is sorted by column 4 and then by column 1. If `FALSE`, the output is not sorted.
#'
#' @return This function does not return a value. It writes the GFF formatted data to a file specified by `gff.file`.
#'
#' @examples
#' # Example usage (assuming `blast_results` is your data frame with BLAST results):
#' blastres2gff(blast_results, "output.gff")
#' 
#' @author Anna A. Igolkina 
blastres2gff <- function(v.blast, gff.file, to.sort = T){
  v.gff = data.frame(col1 = v.blast$V8,
                     col2 = 'blast2gff',
                     col3 = 'query',
                     col4 = v.blast$V4,
                     col5 = v.blast$V5,
                     col6 = '.',
                     col7 = v.blast$strand,
                     col8 = '.',
                     col9 = paste0('ID=Q', 1:nrow(v.blast),
                                   ';query=',v.blast$V1,
                                   ';length=', v.blast$len1,
                                   ';similarity=', round(v.blast$V6, 1),
                                   ';coverage=', round(v.blast$V7 / v.blast$len1 * 100, 1) ))
  
  # Sorting
  if(to.sort){
    v.gff = v.gff[order(v.gff$col4),]
    v.gff = v.gff[order(v.gff$col1),]
  }
  writeGFF(v.gff, gff.file)
}

#' Find Nestedness in Data
#'
#' This function computes the nestedness of given data, considering strand direction if specified.
#' It rearranges sequence start and end positions based on strand direction, creates unique identifiers,
#' and computes coverage for each side of the data.
#'
#' @param v.res A data frame that should contain the following columns:
#'   - V1: An identifier column for sequences #1.
#'   - V2: A numeric column representing the start position of the sequences #1.
#'   - V3: A numeric column representing the end position of the sequences #1.

#'   - V8: An identifier column for sequences #2.
#'   - V4: A numeric column representing the start position of the sequences #2.
#'   - V5: A numeric column representing the end position of the sequences #2.

#' @param use.strand Logical, if TRUE, strand information is considered in the processing.
#'
#' @return Returns a data frame `v.cover` with the following structure:
#'   - C1, C8: Coverage data calculated for each side.
#'   - V1, V8:Identifiers of sequences..
#'   - dir: A column specifying the direction of the strand ('+', '-') if `use.strand` is TRUE, or '.'.
#'   - Additional columns might be included based on the implementation of `getOneSideCoverage` and other calculations within the function.
#'
#' @examples
#' # Example usage:
#' # result <- findNestedness(data, TRUE)
#'
#' @export
findNestedness <- function(v.res, use.strand = T){
  
  idx.strand = v.res$V4 > v.res$V5
  tmp = v.res$V4[idx.strand]
  v.res$V4[idx.strand] = v.res$V5[idx.strand]
  v.res$V5[idx.strand] = tmp
  
  
  # Make unique names for strands
  if(use.strand == T){
    str.strand = c('+', '-')
    dir.strand = str.strand[idx.strand * 1 + 1]
    v.res$V8 = paste(v.res$V8, dir.strand, sep = '|')
  } 
  v.res$comb = paste(v.res$V1, v.res$V8, sep = '||')
  
  v.res$V1 = v.res$comb
  v.res$V8 = v.res$comb
  cover1.info = getOneSideCoverage(v.res)
  cover1 = cover1.info$coverage
  cover8.info = getOneSideCoverage(v.res, side = 1)
  cover8 = cover8.info$coverage
  
  
  v.cover = data.frame(C1 = cover1)
  v.cover$C8 = cover8[rownames(v.cover)]
  v.cover[,c('V1', 'V8')] = stringr::str_split_fixed(rownames(v.cover), "\\|\\|", 2)
  rownames(v.cover) = NULL
  
  if(use.strand){
    v.cover$dir = substr(v.cover$V8, nchar(v.cover$V8), nchar(v.cover$V8))
    v.cover$V8 = substr(v.cover$V8, 1, (nchar(v.cover$V8) - 2))
  } else {
    v.cover$dir = '.'
  }
  
  return(v.cover)
  # v.nest = v.cover[(v.cover$p1 >= sim.cutoff) | (v.cover$p8 >= sim.cutoff),]
  
}

#' Calculate Coverage one sequence.
#'
#' @param v.rest A data frame that should contain the following columns:
#'   - V1: An identifier column for sequences #1.
#'   - V2: A numeric column representing the start position of the sequences #1.
#'   - V3: A numeric column representing the end position of the sequences #1.

#'   if `side = 1` is provided, then also:
#'   - V8: An identifier column for sequences #2.
#'   - V4: A numeric column representing the start position of the sequences #2.
#'   - V5: A numeric column representing the end position of the sequences #2.

#' @param side An integer indicating which side to consider for coverage calculation of #1 or of #2. 
#'             Defaults to 0, so coverage of sequences "#1.
#'
#' @return Returns a named vector where each name corresponds to an identifier in the V1 column,
#'         and each value is the sum of coverage values for that identifier.
#'
#' @export
getOneSideCoverage <- function(v.rest, side = 0){
  
  # print(head(v.rest))
  if(side == 1){
    
    
    idx.tmp = v.rest$V4 > v.rest$V5
    if(sum(idx.tmp) > 0){
      tmp = v.rest$V4[idx.tmp]
      v.rest$V4[idx.tmp] = v.rest$V5[idx.tmp]
      v.rest$V5[idx.tmp] = tmp
    }
    
    v.rest$V1 = v.rest$V8
    v.rest$V2 = v.rest$V4
    v.rest$V3 = v.rest$V5
    
  }
  v.rest = v.rest[, 1:3]
  # print(head(v.rest))
  
  # - - - - - - - - - - - - - - - - - - - - - - - - 
  
  v.rest = v.rest[order(-v.rest$V3),]
  v.rest = v.rest[order(v.rest$V2),]
  v.rest$V1 = as.factor(v.rest$V1)
  v.rest = v.rest[order(v.rest$V1),]
  
  # remove nestedness
  idx.nested = 1
  while(length(idx.nested) > 0){
    idx.nested = which((v.rest$V2[-1] >= v.rest$V2[-nrow(v.rest)]) & 
                         (v.rest$V3[-1] <=v.rest$V3[-nrow(v.rest)]) & 
                         (v.rest$V1[-1] == v.rest$V1[-nrow(v.rest)])) + 1
    # print(length(idx.nested))
    if(length(idx.nested) == 0) next
    v.rest = v.rest[-idx.nested,]  
  }
  
  v.rest$cover = v.rest$V3 - v.rest$V2 + 1
  v.rest$overlap = c(v.rest$V2[-1] - v.rest$V3[-nrow(v.rest)] - 1, 0)
  v.rest$overlap[v.rest$overlap > 0] = 0
  idx.diff = which(v.rest$V1[-1] != v.rest$V1[-nrow(v.rest)])
  v.rest$overlap[idx.diff] = 0
  v.rest$cover = v.rest$cover + v.rest$overlap
  
  coverage = tapply(v.rest$cover, v.rest$V1, sum)
  
  beg = tapply(v.rest$cover, v.rest$V2, min)
  end = tapply(v.rest$cover, v.rest$V3, max)
  return(list(coverage = coverage,
              beg = beg,
              end = end))
}
















