#' Extract chromosome number by a provided format
#' 
#' This function matches a specified chromosome format 
#' in the chromosome names of a GFF file,
#' extracts and saves the corresponding chromosome numbers.
#' 
#' @param gff A data frame containing the GFF file data.
#' @param s.chr The chromosome format to match. For example, "chr" for formats like "chr1", "chr2", etc.
#' @return A modified data frame with chromosome information extracted.
#' 
#' @examples
#' # Example usage
#' gff_data <- read.table("data.gff", header = FALSE)
#' matchChromosomeFormat(gff_data, "chr") # Matches chromosome format starting with "chr"
#' 
#' @export
extractChrByFormat <- function(gff, s.chr){
  
  n.gff = nrow(gff)
  matched.chr.format <- grep(paste0(".*",s.chr,"\\d+"), gff$V1, value = TRUE)
  if(length(matched.chr.format) < 0.7 * n.gff) stop('Check the chromosome formats (#1)')
  
  # Cheromosome ID
  idx.chr.format <- grep(paste0(".*",s.chr,"\\d+"), gff$V1, value = F)
  gff = gff[idx.chr.format,]
  if(nrow(gff) < 0.7 * n.gff) stop('Check the chromosome formats (#2)')
  
  gff$chr <- sub(paste(".*",s.chr,"(\\d+)", sep = ''), "\\1", gff$V1)
  if(sum(is.na(gff$chr)) > 0) stop('Check the chromosome formats (#3)')
  
  if(nrow(gff) != n.gff){
    pokazAttention('Not all the gff records were recognized, only',nrow(gff), 'out of', n.gff)
  }
  
  return(gff)
}


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
#' @param aln.type String specifying the prefix or type of alignment data being used; this might include prefixes such as 'msa_' for generic MSA or 'comb_' for combined reference-based alignments.
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
gff2gff <- function(path.cons, 
                    acc1, acc2, # if one of the accessions is called 'pangen', then transfer is with pangenome coordinate
                    gff1, 
                    ref.acc='0', 
                    n.chr = 5,
                    exact.match=T, 
                    gr.accs.e = "accs/",
                    aln.type = 'msa_',  # please provide correct prefix. For example, in case of reference-based, it's 'comb_'
                    echo=T,
                    pangenome.name='Pangen',
                    s.chr = '_Chr' # in this case the pattern is "*_ChrX", where X is the number
                    ){
  
  # Set of names of accettions, which can be used to specify pangenomes coordinates
  pangenome.names = unique(c(pangenome.name, 'Pangen', 'Pangenome', 'Pannagram'))
  
  if(acc1 == acc2) stop('Accessions provided are identical')
  
  gff1$idx = 1:nrow(gff1)
  # Get chromosome number from the first column of the gff file
  # if(!('chr' %in% colnames(gff1))){
    gff1 =  extractChrByFormat(gff1, s.chr)
    gff1 = gff1[order(gff1$chr),]
  # }
  
  colnames.1.to.9 = colnames(gff1)[1:9]
  colnames(gff1)[1:9] = paste0('V', 1:9)
  # Prepare new annotation
  gff2 = gff1
  gff2$len.init = gff2$V5 - gff2$V4 + 1
  gff2$V9 = paste0(gff2$V9, ';len_init=', gff2$len.init)
  gff2$V2 = 'pannagram'
  gff2$V4 = -1
  gff2$V5 = -1
  gff2$V1 = gsub(acc1, acc2, gff2$V1)
  
  
  for(i.chr in 1:n.chr){
    
    # ---
    # If there some regions to annotate from the chromosome i.chr
    idx.chr = which(gff2$chr == i.chr)
    if(length(idx.chr) == 0) next
    # ---
    
    if(echo) pokaz('Chromosome', i.chr)
    file.msa = paste0(path.cons, aln.type, i.chr, '_', i.chr, '_ref_', ref.acc,'.h5')
    
    if(tolower(acc1) %in% tolower(pangenome.names)){
      v = h5read(file.msa, paste0(gr.accs.e, acc2))
      v = cbind(1:length(v), v)
    } else if (tolower(acc2) %in% tolower(pangenome.names)){
      v = h5read(file.msa, paste0(gr.accs.e, acc1))
      v = cbind(v, 1:length(v))
    } else {  # Two different accessions
      v = cbind(h5read(file.msa, paste0(gr.accs.e, acc1)),
                h5read(file.msa, paste0(gr.accs.e, acc2)))  
    }
    
    max.chr.len = max(nrow(v), max(abs(v[!is.na(v)])))
    
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
    
    
    if(echo) pokaz('Number of annotations:', length(idx.chr))
    
    if(exact.match){
      # Exact match of positions
      gff2$V4[idx.chr] = v.corr[gff1$V4[idx.chr]]
      gff2$V5[idx.chr] = v.corr[gff1$V5[idx.chr]]
    } else {
      w = getPrevNext(v.corr)
      gff2$V4[idx.chr] = w$v.next[gff1$V4[idx.chr]]
      gff2$V5[idx.chr] = w$v.prev[gff1$V5[idx.chr]]
    }
    
    # TODO: add information about blocks
    
  }
  
  idx.wrong.blocks = sign(gff2$V4 * gff2$V5) != 1
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
  
  gff2.remain$len.new = gff2.remain$V5 - gff2.remain$V4 + 1
  gff2.remain$V9 = paste0(gff2.remain$V9, ';len_new=', gff2.remain$len.new)  # add new length
  gff2.remain$V9 = paste(gff2.remain$V9, ';len_ratio=', 
                         round(gff2.remain$len.new / gff2.remain$len.init, 2), sep='')  # add new length
  
  colnames(gff2.loosing)[1:9] = colnames.1.to.9
  colnames(gff2.remain)[1:9] = colnames.1.to.9
  return(gff2.remain = gff2.remain[,1:9])
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
#' @param ... Additional parameters passed to the gff2gff function, such as:
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
#' converted_bed <- bed2bed("path/to/consensus/", "acc1", "acc2", bed_data)
#'
#' @export
bed2bed <- function(path.cons, 
                    acc1, acc2,
                    bed1, 
                    ... # Use '...' to capture all other arguments
) {
  
  # Convert BED to GFF-like format
  colnames.bed1 = colnames(bed1)
  colnames(bed1) = c('chrom', 'beg', 'end', 'name', 'score', 'strand')
  gff1 = data.frame(V1 = bed1$chrom,
                    V2 = 'tmp',
                    V3 = 'type',
                    V4 = bed1$beg,
                    V5 = bed1$end,
                    V6 = bed1$score,
                    V7 = bed1$strand,
                    V8 = '.',
                    V9 = bed1$name
  )
  
  # Call gff2gff function, passing all additional parameters through '...'
  gff2 = gff2gff(path.cons = path.cons, 
                 acc1 = acc1, 
                 acc2 = acc2,
                 gff1 = gff1, 
                 ...)
  
  # Convert the output back to BED format
  bed2 = gff2[,c(1, 4, 5, 9, 6, 7)]
  colnames(bed2) = colnames.bed1[1:6]
  
  return(bed2)
}


pos2pos <- function(path.cons, 
                    acc1, acc2,
                    pos1, 
                    ... # Use '...' to capture all other arguments
) {
  
  # Convert Positions to GFF-like format
  colnames.pos1 = colnames(pos1)
  colnames(pos1) = c('chrom', 'beg', 'end')
  gff1 = data.frame(V1 = pos1$chrom,
                    V2 = 'tmp',
                    V3 = 'type',
                    V4 = pos1$beg,
                    V5 = pos1$end,
                    V6 = 0,
                    V7 = '+',
                    V8 = '.',
                    V9 = 'pos'
  )
  
  # Call gff2gff function, passing all additional parameters through '...'
  gff2 = gff2gff(path.cons = path.cons, 
                 acc1 = acc1, 
                 acc2 = acc2,
                 gff1 = gff1, 
                 ...)
  
  # Convert the output back to Positions format
  pos2 = gff2[,c(1, 4, 5)]
  colnames(pos2) = colnames.pos1[1:3]
  
  return(pos2)
}


plotMsaFragment <- function(path.cons, 
                       acc,
                       chr,
                       beg,
                       end,
                       diff.mode=F,
                       ref.acc = '0', 
                       exact.match=T,
                       gr.accs.e = "accs/",
                       aln.type = 'msa_',  # please provide correct prefix. For example, in case of reference-based, it's 'comb_'
                       echo=T,
                       pangenome.name='Pangen',
                       s.chr = '_Chr' # in this case the pattern is "*_ChrX", where X is the number
                       ){

  pokazAttention('Better not to use `plotMsaFragment`, but use `getMxFragment` together with `msaplot` ')
  seq.mx = getMxFragment(path.cons = path.cons, 
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

getMxFragment <- function(path.cons, 
                            acc,
                            chr,
                            beg,
                            end,
                            ref.acc = '0', 
                            exact.match=T,
                            gr.accs.e = "accs/",
                            aln.type = 'msa_',  # please provide correct prefix. For example, in case of reference-based, it's 'comb_'
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
    file.msa = paste0(path.cons, aln.type, i.chr, '_', i.chr, '_ref_',ref.acc,'.h5')
    v.acc = h5read(file.msa, paste0(gr.accs.e, acc))
    pos1.acc = which(v.acc == pos1)
    pos2.acc = which(v.acc == pos2)  
  }
  
  if(length(pos1.acc) == 0) stop('First position doesn’t exist')
  if(length(pos2.acc) == 0) stop('Second position doesn’t exist')
  
  # Get Alignment
  file.seq.msa = paste0(path.cons, 'seq_', i.chr, '_', i.chr, '_ref_',ref.acc,'.h5')
  
  h5ls(file.seq.msa)
  
  # Get accession names
  groups = h5ls(file.seq.msa)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # Initialize vector and load MSA data for each accession
  seq.mx = matrix('-', nrow = length(accessions), ncol = pos2.acc - pos1.acc + 1)
  
  s.verbose = c('|', rep('-', length(accessions)), '|\n')
  if(echo) cat(paste0(s.verbose, collapse = ''))
  if(echo) cat('|')
  for(i.acc in 1:length(accessions)){
    # pokaz(accessions[i.acc])
    if(echo) cat('.')
    s.acc = h5read(file.seq.msa, paste0(gr.accs.e, accessions[i.acc]))
    
    # if(length(s.acc) != length(v.acc))  stop('MSA and seq do not match')
    seq.mx[i.acc,] = s.acc[pos1.acc:pos2.acc]
  }
  rownames(seq.mx) = accessions
  if(echo) cat('|\n')
  
  return(seq.mx)
}



#' Fill a vector with 1 corresponding to Begin and End positions in GFF Annotations
#' 
#' @param len.acc Integer, length of the vector.
#' @param gff Data frame with `beg` and `end` columns for GFF annotations.
#' 
#' @return Numeric vector with marks for GFF annotation regions.
#' 
#' @examples
#' fillBegEnd(100, data.frame(beg = c(10, 30), end = c(20, 40)))
#' @export
fillBegEnd <- function(len.acc, gff){
  g = rep(0, len.acc)
  
  gff[gff$end <= len.acc,] # remove genes which are over the aligned regions
  
  # Check for the overlap
  diff = gff$beg[-1] < gff$end[-nrow(gff)]
  if(sum(diff) > 0) {
    pokazAttention(gff[min(which(diff)) + (0:1),])
    stop(paste('Overlap in GFF, accession', acc)) 
  }
  
  if(!('idx' %in% colnames(gff))){
    gff$idx = 1:nrow(idx)
  }
  g[gff$beg] = gff$idx
  g[gff$end] = gff$idx * (-1)
  g = cumsum(g)
  g[gff$end] = gff$idx
  return(g)
}



#' Identify Included and Overlapping Annotations in GFF Data
#'
#' This function processes a GFF data frame to find and distinguish between 'include' relationships 
#' (where one annotation is completely within another of the same type) and overlaps (where annotations 
#' of the same type intersect but neither is completely within the other). It iteratively removes identified 
#' includes and overlaps from the dataset, returning a list with indices of includes, overlaps, and the 
#' remaining GFF data.
#'
#' @param gff.cut A data frame representing GFF annotations, expected to contain `beg`, `end`, and `V1` 
#' columns, where `beg` and `end` represent the start and end of the annotations, respectively, and 
#' `V1` indicates the annotation type. The data frame should also contain an `idx` column, which uniquely 
#' identifies each annotation.
#'
#' @return A list containing three elements: `idx.include`, a data frame with `child` and `parent` indices 
#' indicating include relationships; `idx.overlap`, a vector of indices indicating overlaps; and `gff.cut`, 
#' the modified input data frame with includes and overlaps removed.
#'
#' @export

findIncludeAndOverlap <- function(gff, echo = F){
  
  gff.cut = gff
  
  # Find "Include"
  idx.include = c()
  while(T){
    n.gff.cut = nrow(gff.cut)
    idx = which((gff.cut$beg[-1] < gff.cut$end[-n.gff.cut]) &
                  (gff.cut$end[-1] <= gff.cut$end[-n.gff.cut]) &
                  (gff.cut$V1[-1] == gff.cut$V1[-n.gff.cut]))
    if(length(idx) == 0) break
    idx = idx + 1
    idx = setdiff(idx, idx + 1)
    if (echo) pokaz(length(idx))
    idx.include = rbind(idx.include, cbind(gff.cut$idx[idx], gff.cut$idx[idx-1]))
    gff.cut = gff.cut[-idx,]
  }
  idx.include = data.frame(child = idx.include[,1], parent = idx.include[,2])
  idx.include = idx.include[order(idx.include$parent),]
  
  # Found Overlap
  idx.overlap = c()
  while(T){
    n.gff.cut = nrow(gff.cut)
    idx = which((gff.cut$beg[-1] < gff.cut$end[-n.gff.cut]) & 
                  (gff.cut$V1[-1] == gff.cut$V1[-n.gff.cut]))
    if(length(idx) == 0) break
    idx = idx + 1
    idx = setdiff(idx, idx + 1)
    if (echo) pokaz(length(idx))
    idx.overlap = c(idx.overlap, gff.cut$idx[idx])
    gff.cut = gff.cut[-idx,]
  }
  
  pokaz('Remained types:', unique(gff.cut$type))
  pokaz('Inlcude types:', unique(gff[gff$idx %in% idx.include$child,]$type))
  pokaz('Overlap types:', unique(gff[gff$idx %in% idx.overlap,]$type))
  
  return(list(idx.include = idx.include,
              idx.overlap = idx.overlap,
              gff.cut = gff.cut))
}


#' Check for Translocations in Numeric Vector
#'
#' Evaluates a numeric vector for translocations by identifying repeated 
#' non-zero elements after removing consecutive duplicates.
#'
#' @param g Numeric vector to be checked for translocations.
#' 
#' @return NULL
#' 
#' @examples
#' g <- c(0, 1, 2, 2, 3, 1, 4) # This will trigger a translocation error
#' checkTranslocations(g)
#' 
#' @export
checkTranslocations <- function(g, raise.error = F){
  g = g[g != 0]
  idx.diff = c(1, which(diff(g) != 0)+1)
  g = g[idx.diff]
  
  if(length(g) != length(unique(g))) {
    if(raise.error) stop('Translocation')
    g.trans = uniqueDuplicates(g)
    pokazAttention('Translocation was found,', length(g.trans), 'gene affected')
    return(g.trans)
  }
  return(NULL)
}




#' 
getGeneBlocks <- function(g.tmp, len.pan, v.acc){
  
  idx.fill = (v.acc != 0)
  g = rep(0, len.pan)
  g[abs(v.acc[idx.fill])] = g.tmp[idx.fill]
  g.trans = checkTranslocations(g)
  if(length(g.trans) > 0){
    g[g %in% g.trans] = 0
  }
  
  g = cbind(g, 1:length(g))
  g = g[g[,1] != 0,]
  
  idx.beg = c(1, which(diff(g[,1]) != 0)+1)
  idx.end = c(which(diff(g[,1]) != 0), nrow(g))
  
  if(sum(g[idx.beg,1] != g[idx.end,1])) stop("GENE IDs DONT MATCH")
  
  return(data.frame(beg = g[idx.beg,2],
                    end = g[idx.end,2],
                    idx = g[idx.beg,1]))
}


saveVCF <- function(snp.val, snp.pos, chr.name, file.vcf, append=F) {
  
  if(length(snp.pos) != nrow(snp.val)) stop('Dimentions of the SNP matrix and the vectop of positions should match')
  
  if(!append){
    # Open a connection to the output file
    file.vcf.conn <- file(file.vcf, "w")
    # Print the VCF header
    cat("##fileformat=VCFv4.2\n", file = file.vcf.conn)
    sample_names <- colnames(snp.val) # Exclude the first column (reference)
    cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", paste(sample_names, collapse="\t"), "\n", sep="\t", file = file.vcf.conn)
  } else {
    file.vcf.conn <- file(file.vcf, "a")
  }
  
  for (i in 1:nrow(snp.val)) {
    ref <- snp.val[i, 1] # Reference nucleotide
    if(ref == '-'){
      alts.cnt = table(snp.val[i,snp.val[i,] != '-'])
      ref = names(alts.cnt[alts.cnt == max(alts.cnt)])[1]
    }
    alts <- unique(snp.val[i,]) # Unique alternative nucleotides, excluding reference
    alts <- alts[alts != "-" & alts != ref] # Exclude dashes and reference nucleotides
    alt <- paste(alts, collapse = ",") # Concatenate alternative nucleotides
    
    # Skip if there are no alternative alleles
    if (length(alts) == 0) next
    
    # Prepare genotypes for each accession
    genotypes <- sapply(1:ncol(snp.val), function(j) {
      allele <- snp.val[i, j]
      if (allele == ref) {
        return("0/0") # Homozygous reference
      } else if (allele == "-") {
        return("./.") # Unknown genotype
      } else {
        alt_index <- match(allele, alts)
        return(paste(alt_index, alt_index, sep="/")) # Homozygous alternative
      }
    })
    
    # Output VCF line with genotypes
    cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
                chr.name,  # Chromosome
                snp.pos[i],  # Position
                '.',    # SNP ID, placeholder used here
                ref,    # Reference nucleotide
                alt,    # Concatenate alternative nucleotides
                '.',    # Quality, placeholder used here
                'PASS', # Filter, placeholder used here
                '.',    # Additional information, placeholder used here, 
                "GT",   # Genotype format, using GT to denote genotype
                paste(genotypes, collapse="\t")),  # In accessions
        file = file.vcf.conn)
  }
  
  # Close the file connection
  close(file.vcf.conn)
}





