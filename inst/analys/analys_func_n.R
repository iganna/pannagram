# ============================================================================
#  Interval-format (_n) counterparts of the analys_func.R readers.
#
#  Only the functions that read /accs are duplicated here; everything else in
#  analys_func.R is format-independent and is reused as is. Names carry the `_n`
#  suffix so both versions can live in the package at the same time.
#
#  See utils/interval_func.R for the format itself.
# ============================================================================

#' Convert GFF Annotations Using MSA Data
#'
#' This function converts GFF annotations using multiple sequence alignment (MSA) data.
#' It processes data for a specified number of chromosomes and modifies annotations according to
#' the mappings in the MSA.
#'
#' @param path.cons String path to MSA files.
#' @param acc1 String representing the first accession (which has the annotation).
#' @param acc2 String representing the second accession (which annotation you want to produce).
#' @param gff1 DataFrame of GFF annotations of the first accession to be converted.
#' @param n.chr Number of chromosomes to process (default is 5).
#' @param exact.match Logical flag determining the use of exact position matching (default is TRUE).
#' @param max.chr.len Maximum chromosome length (default is 35 * 10^6).
#' @param gr.accs.e Path to access data in the HDF5 file (default is "accs/").
#' @param echo Logical flag for messages during execution (default is TRUE).
#' @param ref.acc String identifier for the accession, which was used to sort the MSA position order.
#' @param aln.type String specifying the prefix or type of alignment data being used; this might include prefixes such as 'pan' for generic MSA or 'comb_' for combined reference-based alignments.
#' @param pangenome.name String specifying a name to represent the pangenomic coordinate system if one of the accessions is named 'pangen', enabling transfers with pangenome coordinates.
#' @param s.chr String used to split the chromosomal name and extract the chromosome number, following a specific pattern such as '*_ChrX' where 'X' denotes a chromosome number.
#' 
#' @return A list with three elements: gff2.loosing (lost annotations),
#' gff2.remain (remaining annotations), idx.remain (indexes of remaining annotations).
#' 
#' @examples
#' library(rhdf5)
#' # Example of function usage:
#' gff_data <- gffgff("path/to/consensus/", "acc1", "acc2", gff1)
#' 
#' @export
gff2gff_n <- function(acc1, acc2, # if one of the accessions is called 'pangen', then transfer is with pangenome coordinate
                    gff1, 
                    path.proj=NULL,
                    aln.type = 'pan',  # please provide correct prefix. For example, in case of reference-based, it's 'comb_'
                    ref.acc='',
                    exact.match=T, 
                    s.chr = '_Chr', # in this case the pattern is "*_ChrX", where X is the number
                    echo=FALSE,
                    pangenome.name='Pangen',
                    remain=F, ...
                    ){
  
  # --- Determine Pannagram version and corresponding paths ---
  dot.args <- list(...)
  pannagram.paths = getPannagramPaths(path.proj, dot.args)
  path.cons <- pannagram.paths$path.msa
  
  # --- Variables ---
  gr.accs.e = "accs/"
  
  # Set of names of accettions, which can be used to specify pangenomes coordinates
  pangenome.names = unique(c(pangenome.name, 'Pangen', 'Pangenome', 'Pannagram'))
  colnames.full1 = colnames(gff1)
  
  if(ref.acc == ''){
    ref.suff = ref.acc
  } else {
    # ref.suff = paste0('_ref_', ref.acc)
    ref.suff = paste0('_', ref.acc)
  }
  
  if(ref.acc != ''){
    aln.type = 'ref'
  }
  
  if(acc1 == acc2) stop('Accessions provided are identical')
  
  # --- Determine the number of chromosomes ---
  
  i.chr <- 1
  while (TRUE) {
    file.msa.tmp <- paste0(path.cons, aln.type, '_', i.chr, '_', i.chr, ref.suff, '.h5')
    if (!file.exists(file.msa.tmp)) break
    i.chr <- i.chr + 1
  }
  n.chr <- i.chr - 1
  if (n.chr == 0) {
    pokaz(path.cons)
    pokaz(aln.type)
    pokazAttention('Required format:', paste0(path.cons, aln.type, '_', 'X_X', ref.suff, '.h5'))
    stop('Files in the required format do not exist.')
  }
  if(echo) pokaz("Number of chromosome is", n.chr)
 
  # Remove zero or negative positions
  idx_negative <- ((gff1$V4 <= 0) | (gff1$V5 <= 0))
  if (any(idx_negative)) {
    pokazAttention("Some positions in the input data are <= 0. They will be removed")
    gff1 <- gff1[!idx_negative, ]
  }
  
  # Indexing
  gff1$idx = 1:nrow(gff1)
  
  # Get chromosomes by format
  gff1 =  extractChrByFormat(gff1, s.chr)
  gff1 = gff1[order(gff1$chr),]

  # Fitler out blocks
  gff1 = filterBlocks_n(acc1, gff1, pangenome.names, n.chr, path.cons, aln.type, ref.suff, gr.accs.e)
  
  # Construct 2
  colnames.1.to.9 = colnames(gff1)[1:9]
  colnames(gff1)[1:9] = paste0('V', 1:9)
  # Prepare new annotation
  gff2 = gff1
  gff2$len.init = gff2$V5 - gff2$V4 + 1
  gff2$V9 = paste0(gff2$V9, ';len_init=', gff2$len.init)
  gff2$V2 = 'pannagram'
  gff2$V4 = 0
  gff2$V5 = 0
  gff2$V1 = gsub(acc1, acc2, gff2$V1)
  
  for(i.chr in 1:n.chr){
    pokaz('Chromosome', i.chr)
    
    # ---
    # If there some regions to annotate from the chromosome i.chr
    idx.chr = which(gff2$chr == i.chr)
    if(length(idx.chr) == 0) next
    # ---
    
    if(echo) pokaz('Chromosome', i.chr)
    file.msa = paste0(path.cons, aln.type, '_', i.chr, '_', i.chr, ref.suff, '.h5')
    
    if(tolower(acc1) %in% tolower(pangenome.names)){
      v = h5VecRead(file.msa, acc2)
      v = cbind(1:length(v), v)
    } else if (tolower(acc2) %in% tolower(pangenome.names)){
      v = h5VecRead(file.msa, acc1)
      v = cbind(v, 1:length(v))
    } else {  # Two different accessions
      v = cbind(h5VecRead(file.msa, acc1),
                h5VecRead(file.msa, acc2))  
    }
    
    max.chr.len = max(nrow(v), max(abs(v[!is.na(v)])))
    idx.chr = idx.chr[gff1$V5[idx.chr] <= max.chr.len]
    
    v = v[v[,1]!=0,]
    v = v[!is.na(v[,1]),]
    v = v[!is.na(v[,2]),]
    idx.v.neg = which(v[,1] < 0)
    if(length(idx.v.neg) > 0){
      v[idx.v.neg,] = v[idx.v.neg,] * (-1)
    }
    # v = v[order(v[,1]),]  # Not necessarily 
    
    # Get correspondence between two accessions
    v.corr = rep(0, max.chr.len)
    v.corr[v[,1]] = v[,2]
    
    if(echo) pokaz('Number of provided annotations:', length(idx.chr))
    
    if(exact.match){
      # Exact match of positions
      gff2$V4[idx.chr] = v.corr[gff1$V4[idx.chr]]
      gff2$V5[idx.chr] = v.corr[gff1$V5[idx.chr]]
    } else {
      w = getPrevNext(v.corr)
      w$v.next[is.na(w$v.next)] = 0
      w$v.prev[is.na(w$v.prev)] = 0
      
      gff2$V4[idx.chr] = w$v.next[gff1$V4[idx.chr]]
      gff2$V5[idx.chr] = w$v.prev[gff1$V5[idx.chr]]
      
      # save(list = ls(), file = "tmp_workspace_gff2gff.RData")
      
      # If ids > length of the alignment, it will introduce NAs
      # gff2$V4[is.na(gff2$V4)] = 0
      # gff2$V5[is.na(gff2$V5)] = 0
      
      # if(sum(is.na(gff2$V4[idx.chr]) > 0)) stop('NA in gff2$V4[idx.chr]')
      # if(sum(is.na(gff2$V5[idx.chr]) > 0)) stop('NA in gff2$V5[idx.chr]')
    }
    idx.chr = idx.chr[(gff2$V4[idx.chr] != 0) & (gff2$V5[idx.chr] != 0)]
    
    # Information about blocks
    blocks = getSequntialBlocks(v[,2])
    bl.beg = blocks[abs(gff2$V4[idx.chr])]
    bl.end = blocks[abs(gff2$V5[idx.chr])]
    
    # save(list = ls(), file = "tmp_workspace_gff2gff.RData")
    
    idx.noneq = (bl.beg != bl.end) *1
    idx.noneq[is.na(idx.noneq)] = 1
    
    if(sum(idx.noneq) > 0){
      gff2$V4[idx.chr[idx.noneq == 1]] = 0
      gff2$V5[idx.chr[idx.noneq == 1]] = 0  
    }
    
    if(sum(is.na(gff2$V4[idx.chr])) > 0) stop('NA in gff2$V4[idx.chr]')
    if(sum(is.na(gff2$V5[idx.chr])) > 0) stop('NA in gff2$V5[idx.chr]')
  }
  
  idx.wrong.blocks = (sign(gff2$V4 * gff2$V5) != 1)
  
  # gff2[idx.wrong.blocks, c('V4', 'V5')] = 0
  gff2.loosing = gff2[idx.wrong.blocks,]
  gff2.remain = gff2[!idx.wrong.blocks,]
  
  # Fix direction
  s.strand = c('+'='-', '-'='+')
  idx.neg = gff2.remain$V4 < 0
  tmp = abs(gff2.remain$V4[idx.neg])
  # print(head(gff2.remain[idx.neg,]))
  gff2.remain$V4[idx.neg] = abs(gff2.remain$V5[idx.neg])
  gff2.remain$V5[idx.neg] = tmp
  gff2.remain$V7[idx.neg] = s.strand[gff2.remain$V7[idx.neg]]
  
  idx.trans = which(gff2.remain$V4 > gff2.remain$V5)
  if(length(idx.trans) > 0){
    gff2.remain = gff2.remain[-idx.trans,]
  }
  
  # Fitler out blocks
  gff2.remain = filterBlocks_n(acc2, gff2.remain, pangenome.names, n.chr, path.cons, aln.type, ref.suff, gr.accs.e)
  
  
  # Prepare results
  gff2.remain$len.new = gff2.remain$V5 - gff2.remain$V4 + 1
  gff2.remain$V9 = paste0(gff2.remain$V9, ';len_new=', gff2.remain$len.new)  # add new length
  gff2.remain$V9 = paste(gff2.remain$V9, ';len_ratio=', 
                         round(gff2.remain$len.new / gff2.remain$len.init, 2), sep='')  # add new length
  
  colnames(gff2.loosing)[1:9] = colnames.1.to.9
  colnames(gff2.remain)[1:9] = colnames.1.to.9
  
  
  if(remain){
    return(gff2.remain[,colnames.full1])
  }else {
    return(gff2.remain[,1:9])
  }
}


#' Convert BED Annotations Using MSA Data
#'
#' This function utilizes Multiple Sequence Alignment (MSA) data to convert BED annotations
#' from one genomic accession to another. It processes annotations and modifies them based
#' on the MSA mappings.
#'
#' @param path.cons String specifying the path to the directory containing MSA files.
#' @param acc1 String identifier for the first accession, which contains the original annotations.
#' @param acc2 String identifier for the second accession, to which annotations are being transferred.
#' @param bed1 DataFrame containing the BED annotations of the first accession.
#' @param ... Additional parameters passed to the gff2gff_n function, such as:
#'   - ref.acc: String identifier for the accession used in to sort the positions in MSA.
#'   - n.chr: Integer specifying the number of chromosomes to process.
#'   - exact.match: Logical flag indicating whether exact position matching should be used.
#'   - gr.accs.e: String specifying the path within the HDF5 file to access data.
#'   - aln.type: String specifying the prefix of alignment files being used.
#'   - echo: Logical flag to enable or disable messages during the function's execution.
#'   - pangenome.name: String specifying a name to represent the pangenomic coordinate system.
#'   - s.chr: String used to split chromosomal names and extract the chromosome number.
#'
#' @return A DataFrame in BED format containing the converted annotations.
#'
#' @examples
#' library(rhdf5) # hypothetical example, adjust according to actual usage
#' bed_data <- read.table("path/to/data.bed", header = TRUE, sep = "\t")
#' converted_bed <- bed2bed_n("path/to/consensus/", "acc1", "acc2", bed_data)
#'
#' @export
bed2bed_n <- function(bed1, 
                    ... # Use '...' to capture all other arguments
) {
  pokaz('bed2bed_n')
  # Convert BED to GFF-like format
  colnames.bed1 = colnames(bed1)
  colnames(bed1) = c('chrom', 'beg', 'end', 'name', 'score', 'strand')
  gff1 = data.frame(V1 = bed1$chrom,
                    V2 = 'tmp',
                    V3 = 'type',
                    V4 = bed1$beg + 1, # counting starts from 0
                    V5 = bed1$end,
                    V6 = bed1$score,
                    V7 = bed1$strand,
                    V8 = '.',
                    V9 = bed1$name
  )
  
  # Call gff2gff_n function, passing all additional parameters through '...'
  gff2 = gff2gff_n(gff1 = gff1, 
                 ...)
  gff2$V4 = gff2$V4 - 1  # counting starts from 0
  
  # Convert the output back to BED format
  bed2 = gff2[,c(1, 4, 5, 9, 6, 7)]
  colnames(bed2) = colnames.bed1[1:6]
  
  return(bed2)
}


#' @export
pos2pos_n <- function(pos1, 
                    ... # Use '...' to capture all other arguments
) {
  
  # Convert Positions to GFF-like format
  colnames.pos1 = colnames(pos1)
  colnames(pos1) = c('chrom', 'pos', 'info')
  gff1 = data.frame(V1 = pos1$chrom,
                    V2 = 'tmp',
                    V3 = 'type',
                    V4 = pos1$pos,
                    V5 = pos1$pos,
                    V6 = 0,
                    V7 = '+',
                    V8 = '.',
                    V9 = pos1$info
  )
  
  # Call gff2gff_n function, passing all additional parameters through '...'
  gff2 = gff2gff_n(gff1 = gff1, 
                 ...)
  
  # Convert the output back to Positions format
  pos2 = gff2[,c(1, 5, 9)]
  colnames(pos2) = colnames.pos1[1:3]
  
  return(pos2)
}


plotMsaFragment_n <- function(path.cons, 
                       acc,
                       chr,
                       beg,
                       end,
                       diff.mode=F,
                       ref.acc = '0', 
                       exact.match=T,
                       gr.accs.e = "accs/",
                       aln.type = 'pan',  # please provide correct prefix. For example, in case of reference-based, it's 'comb_'
                       echo=T,
                       pangenome.name='Pangen',
                       s.chr = '_Chr' # in this case the pattern is "*_ChrX", where X is the number
                       ){

  pokazAttention('Better not to use `plotMsaFragment_n`, but use `getMxFragment_n` together with `msaplot` ')
  seq.mx = getMxFragment_n(path.cons = path.cons, 
                         acc = acc,
                         chr = chr,
                         beg = beg,
                         end = end,
                         ref.acc=ref.acc,
                         exact.match=exact.match,
                         gr.accs.e=gr.accs.e,
                         aln.type=aln.type,
                         echo=echo,
                         pangenome.name=pangenome.name,
                         s.chr=s.chr)
  
  if(diff.mode){
    p3 = msadiff(seq.mx)
    s.diff = '_diff'
  } else {
    p3 = msaplot(seq.mx)
    s.diff = ''
  }  
  

  return(p3)
}


getMxFragment_n <- function(path.cons, 
                            acc,
                            chr,
                            beg,
                            end,
                            ref.acc = '0', 
                            exact.match=T,
                            gr.accs.e = "accs/",
                            aln.type = 'pan',  # please provide correct prefix. For example, in case of reference-based, it's 'comb_'
                            echo=T,
                            pangenome.name='Pangen',
                            s.chr = '_Chr' # in this case the pattern is "*_ChrX", where X is the number
){
  
  gr.accs.e <- "accs/"
  gr.accs.b <- "/accs"
  
  i.chr = chr
  pos1 = beg
  pos2 = end
  
  if(pos1 > pos2){
    pokazAttention('Positions were sorted')
    tmp = pos1
    pos1 = pos2
    pos2 = tmp
  }
  
  # Get positions in the pangenome coordinates
  pangenome.names = unique(c(pangenome.name, 'Pangen', 'Pangenome', 'Pannagram'))
  if(tolower(acc) %in% tolower(pangenome.names)){
    pos1.acc = pos1
    pos2.acc = pos2
  } else {
    file.msa = paste0(path.cons, aln.type, '_', i.chr, '_', i.chr, '_ref_',ref.acc,'.h5')
    ivl.acc = h5IvlRead(file.msa, acc)
    # Binary search instead of a scan over the whole pangenome. The legacy code
    # compared the SIGNED stored value, i.e. it matched plus-strand positions
    # only -- keep that, so both versions accept exactly the same input.
    mapPos = function(q){
      i = ivlPanAt(ivl.acc, q)
      i = i[i > 0]
      if(length(i) == 0) return(integer(0))
      i[ivlToVecRange(ivl.acc, i[1], i[1]) == q]
    }
    pos1.acc = mapPos(pos1)
    pos2.acc = mapPos(pos2)  
  }
  
  if(length(pos1.acc) == 0) stop('First position does not exist')
  if(length(pos2.acc) == 0) stop('Second position does not exist')
  
  # Get Alignment
  file.seq.msa = paste0(path.cons, 'seq_', i.chr, '_', i.chr, '_ref_',ref.acc,'.h5')
  
  rhdf5::h5ls(file.seq.msa)
  
  # Get accession names
  groups = rhdf5::h5ls(file.seq.msa)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # Initialize vector and load MSA data for each accession
  seq.mx = matrix('-', nrow = length(accessions), ncol = pos2.acc - pos1.acc + 1)
  
  s.verbose = c('|', rep('-', length(accessions)), '|\n')
  if(echo) cat(paste0(s.verbose, collapse = ''))
  if(echo) cat('|')
  for(i.acc in 1:length(accessions)){
    # pokaz(accessions[i.acc])
    if(echo) cat('.')
    # Hyperslab instead of "read the whole accession, then subset in R". With
    # seq_*.h5 written in 64K chunks this touches only the chunks the window falls
    # into; the legacy line inflated the entire accession for every call.
    seq.mx[i.acc,] = rhdf5::h5read(file.seq.msa, paste0(gr.accs.e, accessions[i.acc]),
                                   index = list(pos1.acc:pos2.acc))
  }
  rownames(seq.mx) = accessions
  if(echo) cat('|\n')
  
  return(seq.mx)
}


filterBlocks_n <- function(acc, gff, pangenome.names, n.chr, path.cons, aln.type, ref.suff, gr.accs.e) {
  # Convert pangenome.names to lowercase once for efficiency
  pangenome.names <- tolower(pangenome.names)
  
  # Initialize an empty vector for indices to remove
  idx.remove <- c()
  if (!(tolower(acc) %in% pangenome.names)) {
    # Loop over each chromosome
    for (i.chr in 1:n.chr) {
      # Find indices for the current chromosome
      idx.chr <- which(gff$chr == i.chr)
      gff.chr <- gff[idx.chr, ]
      
      # Generate the MSA file name
      file.msa <- paste0(path.cons, aln.type, '_', i.chr, '_', i.chr, ref.suff, '.h5')
      
      # Check if the MSA file exists before reading
      if (file.exists(file.msa)) {
        # The block id is a column of the interval table, so nothing is
        # materialised here: ivlBlockAt() answers by binary search what the
        # legacy code read out of a vector as long as the accession chromosome.
        ivl.acc <- h5IvlRead(file.msa, acc)
        
        bl.beg <- ivlBlockAt(ivl.acc, gff.chr$beg)
        bl.end <- ivlBlockAt(ivl.acc, gff.chr$end)
        
        # defineBlocks() returned a vector of length max(abs(v)). Inside that range
        # but outside every block span the value was 0; only PAST the end did the
        # indexing yield NA. ivlBlockAt() returns 0 for both cases, so the
        # distinction has to be restored: the filter below drops a row on NA but
        # keeps it when both ends are 0, i.e. a feature lying entirely between two
        # blocks survives in the legacy code.
        len.acc.max <- if (nrow(ivl.acc) == 0) 0 else
          max(pmax(ivl.acc$acc.beg,
                   ivl.acc$acc.beg + ivl.acc$delta * (ivl.acc$len - 1L)))
        bl.beg[gff.chr$beg > len.acc.max] <- NA
        bl.end[gff.chr$end > len.acc.max] <- NA
        
        # Append indices where block starts and ends do not match
        idx.remove <- c(idx.remove, idx.chr[which(bl.beg != bl.end)])
        # Append indices where block start or end is NA
        idx.remove <- c(idx.remove, idx.chr[which(is.na(bl.beg) | is.na(bl.end))])
      } else {
        message(paste("File not found:", file.msa))
      }
    }
  }
  
  # Remove rows from gff based on idx.remove if there are any indices
  if (length(idx.remove) > 0) {
    idx.remove <- unique(idx.remove)
    gff <- gff[-idx.remove, ]
  }
  
  # Return the filtered gff
  return(gff)
}
