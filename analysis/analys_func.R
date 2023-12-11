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
#' @param flag.exact Logical flag determining the use of exact position matching (default is TRUE).
#' @param max.chr.len Maximum chromosome length (default is 35 * 10^6).
#' @param gr.accs.e Path to access data in the HDF5 file (default is "accs/").
#' @param echo Logical flag for messages during execution (default is TRUE).
#' 
#' @return A list with three elements: gff2.loosing (lost annotations),
#' gff2.remain (remaining annotations), idx.remain (indexes of remaining annotations).
#' 
#' @examples
#' library(rhdf5)
#' # Example of function usage:
#' gff_data <- gffgff("path/to/consensus/", "acc1", "acc2", gff1)
#' 
gffgff <- function(path.cons, acc1, acc2, ref.acc, gff1, 
                   n.chr = 5,
                   flag.exact=T, 
                   max.chr.len = 35 * 10^6,
                   gr.accs.e = "accs/",
                   echo=T){
  
  # Prepare new annotation
  gff2 = gff1
  gff2$len.init = gff2$V5 - gff2$V4 + 1
  gff2$V9 = paste(gff2$V9, ';len_init=', gff2$len.init, sep='')
  gff2$V2 = 'panConvertor'
  gff2$V4 = -1
  gff2$V5 = -1
  
  for(i.chr in 1:n.chr){
    if(echo) pokaz('Chromosome', i.chr)
    file.msa = paste(path.cons, 'msa_', i.chr, '_', i.chr, '_ref_', ref.acc,'.h5', sep = '')
    
    v = cbind(h5read(file.msa, paste(gr.accs.e, acc1, sep = '')),
              h5read(file.msa, paste(gr.accs.e, acc2, sep = '')))
    v = v[v[,1]!=0,]
    idx.v.neg = which(v[,1] < 0)
    if(length(idx.v.neg) > 0){
      v[idx.v.neg,] = v[idx.v.neg,] * (-1)
    }
    # v = v[order(v[,1]),]  # Not necessarily 
    
    # Get correspondence between two accessions
    v.corr = rep(0, max.chr.len)
    v.corr[v[,1]] = v[,2]
    
    idx.chr = which(gff2$chr == i.chr)
    if(echo) pokaz('Number of annotations:', length(idx.chr))
    
    if(flag.exact){
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
  gff2.remain$V4[idx.neg] = abs(gff2.remain$V5[idx.neg])
  gff2.remain$V5[idx.neg] = tmp
  gff2.remain$V7[idx.neg] = s.strand[gff2.remain$V7[idx.neg]]
  
  gff2.remain$len.new = gff2.remain$V5 - gff2.remain$V4 + 1
  gff2.remain$V9 = paste(gff2.remain$V9, ';len_new=', gff2.remain$len.new, sep='')  # add new length
  gff2.remain$V9 = paste(gff2.remain$V9, ';len_ratio=', 
                         round(gff2.remain$len.new / gff2.remain$len.init, 2), sep='')  # add new length
  
  return(list(gff2.loosing = gff2.loosing, gff2.remain = gff2.remain, idx.remain = which(!idx.wrong.blocks)))
}