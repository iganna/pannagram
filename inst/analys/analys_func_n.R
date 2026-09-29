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
  gff1$idx = seq_len(nrow(gff1))
  
  # Get chromosomes by format
  gff1 =  extractChrByFormat(gff1, s.chr)
  gff1 = gff1[order(gff1$chr),]

  # Fitler out blocks
  gff1 = filterBlocks_n(acc1, gff1, pangenome.names, n.chr, path.cons, aln.type, ref.suff, gr.accs.e,
                        exact.match = exact.match)

  if(nrow(gff1) == 0){
    pokazAttention('No annotations to convert')
    if(remain){
      return(gff1[,colnames.full1])
    }else {
      return(gff1[,1:9])
    }
  }
  
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
    
    v = v[v[,1]!=0,,drop=F]
    v = v[!is.na(v[,1]),,drop=F]
    v = v[!is.na(v[,2]),,drop=F]
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
  s.strand = c('+'='-', '-'='+', '.'='.')
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
  if(nrow(gff2.remain) > 0){
    gff2.remain$V9 = paste0(gff2.remain$V9, ';len_new=', gff2.remain$len.new)  # add new length
    gff2.remain$V9 = paste(gff2.remain$V9, ';len_ratio=', 
                           round(gff2.remain$len.new / gff2.remain$len.init, 2), sep='')  # add new length
  } else {
    pokazAttention('No annotations were converted')
  }
  
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
  for(i.acc in seq_along(accessions)){
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


filterBlocks_n <- function(acc, gff, pangenome.names, n.chr, path.cons, aln.type, ref.suff, gr.accs.e, exact.match = T) {
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
        
        bl.beg <- ivlBlockAt(ivl.acc, gff.chr$V4)
        bl.end <- ivlBlockAt(ivl.acc, gff.chr$V5)

        # Not exact match: an end out of blocks is moved inwards, to the nearest block
        if(!exact.match){
          bl.spans <- ivlBlockSpans(ivl.acc)
          bl.beg <- c(bl.spans$block, NA)[findInterval(gff.chr$V4 - 1, bl.spans$hi) + 1]
          bl.end <- c(NA, bl.spans$block)[findInterval(gff.chr$V5, bl.spans$lo) + 1]
        }

        # Append indices where block starts and ends do not match
        idx.remove <- c(idx.remove, idx.chr[which(bl.beg != bl.end)])
        # Append indices where block start or end is NA
        idx.remove <- c(idx.remove, idx.chr[which(is.na(bl.beg) | is.na(bl.end))])
        # Append indices where block start or end is outside of all blocks:
        # ivlBlockAt() returns 0 both between the blocks and past the end
        idx.remove <- c(idx.remove, idx.chr[which((bl.beg == 0) | (bl.end == 0))])
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


#' Get SNPs from the Pangenome Alignment
#'
#' Library version of the features step `-snp`: finds SNPs on all chromosomes
#' of the alignment and saves them into VCF files, one per chromosome.
#' A position is a SNP if at least one accession differs from
#' the consensus there (gaps are ignored).
#' Uses the files produced by features: `seq_<comb>.h5` and
#' `seq_cons_<comb>.fasta` (consensus folder) and the alignment files `<aln.type>_<comb>.h5`.
#'
#' @param path.proj Path to the project folder.
#' @param acc Accession in whose coordinates SNPs are returned,
#'   or 'pangenome' (also 'pangen', 'pannagram') for pangenome coordinates.
#' @param aln.type Prefix of the alignment files (default 'pan').
#' @param ref.acc Reference accession, only if the alignment is reference-based.
#' @param positions Which positions to return: 'snp' (default) - only variable positions,
#'   'invariant' - only non-variable positions, 'all' - all positions of the alignment.
#'   'invariant' and 'all' are large for whole genomes!
#' @param cores Number of cores to read accessions in parallel (default 1).
#' @param singer If TRUE, the VCF is prepared as input for SINGER (only with positions = 'snp'):
#'   only biallelic positions where all accessions have a nucleotide (no gaps/N),
#'   haploid genotypes `0`/`1`. Run SINGER with `-ploidy 1`.
#'   Files get the prefix `singer_`.
#'
#' @note INTERVAL FORMAT (_n). Reads alignments written by `pannagram_n`, where
#'   /accs/<acc> is an interval table (see utils/interval_func.R). The seq_*.h5 and
#'   consensus files are the same as for `getSNPs()`.
#'
#' @return Invisibly, paths to the saved VCF files. Files are saved into
#'   `<path.proj>/features/snp/` with names
#'   `snps_`/`invariant_`/`all_` + `<comb><ref.suff><aln.suff>_<pangen|acc>.vcf`,
#'   where `<aln.suff>` is empty for the default alignment ('pan') and `_<aln.type>` for the others.
#'   For `acc` = pangenome, REF is the consensus nucleotide; otherwise it is the nucleotide of `acc`.
#'   Positions are processed by chunks to limit memory.
#'
#' @examples
#' \dontrun{
#' getSNPs_n(path.proj = 'path_to_alignment_project/')
#' getSNPs_n(path.proj = 'path_to_alignment_project/', acc = 'name_genome1')
#' getSNPs_n(path.proj = 'path_to_alignment_project/', aln.type = 'synteny', singer = TRUE)
#' }
#'
#' @export
getSNPs_n <- function(path.proj = NULL,
                      acc = 'pangenome',
                      aln.type = 'pan',
                      ref.acc = '',
                      positions = c('snp', 'invariant', 'all'),
                      cores = 1,
                      singer = F){

  positions <- match.arg(positions)
  singer <- isTRUE(as.logical(singer))
  if(singer && positions != 'snp') stop("singer = TRUE works only with positions = 'snp'")

  # --- Paths ---
  if(is.null(path.proj)) stop("Path to the project folder must be provided!")
  path.proj = paste0(sub('/+$', '', path.proj), '/')
  path.msa = paste0(path.proj, 'features/alignments/')
  path.seq = paste0(path.proj, 'features/consensus/')
  path.out = paste0(path.proj, 'features/snp/')

  # --- Variables ---
  gr.accs.e = "accs/"
  gr.accs.b = "/accs"
  pangenome.names = c('pangen', 'pangenome', 'pannagram')
  acc.pangen = tolower(acc) %in% pangenome.names

  ref.suff = if (ref.acc == '') '' else paste0('_', ref.acc)
  aln.pref = paste0(aln.type, '_')
  # The default alignment ('pan') keeps the old names: seq_1_1.h5; the rest get seq_1_1_<aln.type>.h5
  seq.suff = if (aln.type == 'pan') '' else paste0('_', aln.type)

  # --- Combinations (chromosomes) from the alignment files ---
  s.combinations <- list.files(path = path.msa, pattern = paste0("^", aln.pref, ".*\\.h5$"))
  s.combinations = sub(aln.pref, "", s.combinations)
  s.combinations = sub("\\.h5$", "", s.combinations)
  if(ref.suff != ''){
    s.combinations = s.combinations[endsWith(s.combinations, ref.suff)]
    s.combinations = substr(s.combinations, 1, nchar(s.combinations) - nchar(ref.suff))
  }
  s.combinations = s.combinations[grepl("^[0-9]+_[0-9]+$", s.combinations)]
  if(length(s.combinations) == 0){
    pokazAttention('Required format:', paste0(path.msa, aln.pref, 'X_X', ref.suff, '.h5'))
    stop('Files in the required format do not exist.')
  }
  s.combinations = s.combinations[order(sapply(s.combinations, comb2ref))]

  # --- Check that sequence files from features (-seq) exist ---
  files.seq.cons = paste0(path.seq, "seq_cons_", s.combinations, ref.suff, seq.suff, ".fasta")
  files.seq      = paste0(path.seq, "seq_",      s.combinations, ref.suff, seq.suff, ".h5")
  files.missing  = c(files.seq.cons, files.seq)[!file.exists(c(files.seq.cons, files.seq))]
  if(length(files.missing) > 0){
    pokazAttention('Missing files:', files.missing)
    stop('Sequence files do not exist. Run the features step -seq first.')
  }

  # --- Parallel ---
  # The cluster is created anew for each round to free the memory of workers
  if(cores > 1) `%dopar%` <- foreach::`%dopar%`

  # --- Output files ---
  if(!dir.exists(path.out)) dir.create(path.out, recursive = TRUE)
  if(!dir.exists(path.out)) stop(paste("The output folder was not created:", path.out))
  file.pref = c(snp = 'snps_', invariant = 'invariant_', all = 'all_')[positions]
  if(singer) file.pref = 'singer_'
  file.suff = if(acc.pangen) '_pangen' else paste0('_', acc)
  min.alleles = if(positions == 'snp') 2 else 1
  n.chunk = 10^6  # max number of positions in memory at once

  files.vcf = c()
  for (s.comb in s.combinations) {
    pokaz("Combination", s.comb)
    i.chr = comb2ref(s.comb)

    # Get Consensus
    file.seq.cons = paste0(path.seq, "seq_cons_", s.comb, ref.suff, seq.suff, ".fasta")
    s.pangen = seq2nt(readFastaMy(file.seq.cons))

    # Get accessions
    file.seq = paste0(path.seq, "seq_", s.comb, ref.suff, seq.suff, ".h5")
    groups = rhdf5::h5ls(file.seq)
    accessions = groups$name[groups$group == gr.accs.b]
    if(!acc.pangen && !(acc %in% accessions)) stop(sprintf("Accession '%s' is not in the alignment", acc))

    # Round 1: positions of differences
    if(positions == 'all'){
      pos = seq_along(s.pangen)
    } else {
      if(cores == 1){
        pos.diff.list = lapply(accessions, function(acc.tmp) {
          v = rhdf5::h5read(file.seq, paste0(gr.accs.e, acc.tmp))
          which((v != s.pangen) & (v != "-"))
        })
      } else {
        myCluster <- parallel::makeCluster(cores, type = "PSOCK")
        doParallel::registerDoParallel(myCluster)
        pos.diff.list = tryCatch(
          foreach::foreach(acc.tmp = accessions, .errorhandling = "stop") %dopar% {
            v = rhdf5::h5read(file.seq, paste0(gr.accs.e, acc.tmp))
            which((v != s.pangen) & (v != "-"))
          },
          finally = parallel::stopCluster(myCluster))
        rm(myCluster)
      }
      pos = sort(unique(unlist(pos.diff.list, use.names = FALSE)))
      rm(pos.diff.list)
      if(positions == 'invariant') pos = setdiff(seq_along(s.pangen), pos)
    }

    # Positions in the output coordinates
    if(acc.pangen){
      chr.name = paste0("PanGen_Chr", i.chr)
      pos.out = pos
    } else {
      # Convert to the coordinates of the accession
      chr.name = paste0(acc, "_Chr", i.chr)
      file.comb = paste0(path.msa, aln.pref, s.comb, ref.suff, ".h5")
      pos.acc = h5VecRead(file.comb, acc)[pos]
      pos = pos[pos.acc != 0]
      pos.out = abs(pos.acc[pos.acc != 0])
      rm(pos.acc)

      ord = order(pos.out)
      pos = pos[ord]
      pos.out = pos.out[ord]
      rm(ord)
    }

    if (length(pos) == 0) {
      pokaz("No positions were found..")
      next
    }

    # Round 2: nucleotides in positions, by chunks
    file.vcf = paste0(path.out, file.pref, s.comb, ref.suff, seq.suff, file.suff, ".vcf")
    idx.chunks = split(seq_along(pos), ceiling(seq_along(pos) / n.chunk))
    for(i.chunk in seq_along(idx.chunks)){
      pos.chunk = pos[idx.chunks[[i.chunk]]]
      # Only the range of the chunk is read from h5 (start/count: much faster than index = for strings)
      p.beg = min(pos.chunk)
      p.end = max(pos.chunk)
      pos.chunk = pos.chunk - p.beg + 1

      if(cores == 1){
        val.list = lapply(accessions, function(acc.tmp) {
          rhdf5::h5read(file.seq, paste0(gr.accs.e, acc.tmp), start = p.beg, count = p.end - p.beg + 1)[pos.chunk]
        })
      } else {
        myCluster <- parallel::makeCluster(cores, type = "PSOCK")
        doParallel::registerDoParallel(myCluster)
        val.list = tryCatch(
          foreach::foreach(acc.tmp = accessions, .errorhandling = "stop") %dopar% {
            rhdf5::h5read(file.seq, paste0(gr.accs.e, acc.tmp), start = p.beg, count = p.end - p.beg + 1)[pos.chunk]
          },
          finally = parallel::stopCluster(myCluster))
        rm(myCluster)
      }
      snp.val = do.call(cbind, val.list)
      colnames(snp.val) = accessions
      rm(val.list)

      snp.ref = if(acc.pangen) s.pangen[pos.chunk + p.beg - 1] else snp.val[, acc]
      snp.pos = pos.out[idx.chunks[[i.chunk]]]

      saveVCF2(snp.val, snp.pos, chr.name = chr.name, file.vcf = file.vcf,
               append = (i.chunk > 1), snp.ref = snp.ref, min.alleles = min.alleles,
               singer = singer)
      rm(snp.val, snp.ref, snp.pos, pos.chunk)
      gc()
    }
    pokaz('Saved', file.vcf)
    files.vcf = c(files.vcf, file.vcf)

    rm(s.pangen, pos, pos.out, idx.chunks)
    rhdf5::h5closeAll()
    gc()
  }

  return(invisible(files.vcf))
}
