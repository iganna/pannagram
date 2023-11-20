#' Find Open Reading Frames (ORFs) in DNA Sequences
#'
#' This function identifies all possible open reading frames (ORFs) in a given set of DNA sequences in both strands.
#' This function scans the input sequences,translates them into amino acid sequences, 
#' and then identifies the ORFs based on the specified minimum and maximum lengths.
#'
#' @param seqs A character vector of DNA sequences.
#' @param min.len The minimum length of ORF to be considered (in amino acids). Default: `30` (no limit).
#' @param max.len The maximum length of ORF to be considered (in amino acids). Default: `Inf` (no limit).
#'
#' @return A character vector of amino acid sequences corresponding to the identified ORFs.
#' Each sequence name contains the original sequence name, the peptide identifier, start and end positions,
#' strand, and length.
#'
#' @examples
#' # Example DNA sequence
#' seqs <- c("CTAAAAAGCTTCGTTCTCTCCCAGACACGGTAATCAAGACATCTTCTGTAATCATGCTTACCGGCGCCGGTGAGACTCATGTTCGCCGGCGTCGGGACCTGCAAAATCTTCTTAGTCTGTTTCAGTCTTTCTTAAGCTTCTGGTCTTGGTTTTTCTCTCCGTTAAATTCTTCTCTCTCGCAAAATCCAAAACTTTCCTTTCCTTTGGTCTTCTTCTTTGAGGGTCTTCTTTGCGACGGAAGAGTGAATGTTTACTCTACAAGCCCTGCCCCGAATCATAACTCAGACGCAGATCCGCCTTTGAACCCGACAAGAGAAGAATCGCACTTCTCCATCATCGCCCTTTCCAAATTTCTAGAGCAGTTTTTCATCGACACATCCTTTGGCTCCAGATCTGCCCTATTAGCAGTAACGCCTTTAGTCCCTTCTCACCGGATTTCGTTGATGGCAAGAGATTCTCTTGGACCATCAGTTGTGCCTTTTCCGGTGAGAACGTCGACGACACACGACGAGCTTCATTGCTGGCAGGTGTTTTTCCCAGGTGCTTCACCGTCGACCATCACTTCCGATCTTTACGACGGCGACGCTCTCCCTTCCCTTCTCTCAACCGTAACAACTCTAGGGTTTGGTGGGATGGGCCTCAATTTTCTTGACGGCCCAGTTGCAACAATCCATGCTTTTGGTAGTTTGGGCCTCTCTCAATTCGGCGGCCCAAAAACTCTAAACCTAATTTGGCCCATATCCAGGTTTGCTTTACTCAGAGTTTTGTCGACTCTCTCTCTGAACCATTGCTTCACGATAACGGATGGTCGTGGACTCCATTCAATCGTTGGAACCTTAGTCTTAGAGGATCCTTGCAATCACAATACTAGAAACCTGTTACATAGAGTAGATTGGTTGCCAACAGATGTAAGTTTAGTCTTGAAGCTTAGTGGATTAGTTTTGATCTCGAACCATAGGCTATGTGATGCGGCTCACGCCTTTGTTTTTATATGTGATCATATAATGAGTCTGAACCTTCCCACTATCATGCTTACTCATCCCGTCCCATGTCGTTTGTTGGACGATGCATCACAACCCTCTTGGAATTGTTGGCTAAAACTGCTTCGGCCTGCTCTATCAACCCTCTTCTCAAGCCCCATGAATCCGGTTAATCTCTGTCATTTGATTAGAAGATTATCTGTGCTTTGGACTTGTAATCTTTGTTGCCGATTGCATTATGTTTACTCCTCTAAGGAGATACCTTCCCTGAAACCTTACTCCTCTAAGGAGATATCCTCTTTGAAACTTGAACCCGGAGAAAAGTCTCACACCGCAGCCTCCTCTTTGGAGCATTTCAACAAAAGCTTTTATGCGTCGATCGCACGAGCTTACTATGATCACAAGTTGATCCGAGATTTCTCGAAGCTAGTCTTCAAGCTTGAATCTCTGCAGGTCCCAATAGATCTACCATCCCTTGCAACCTTCTCGAAGTTATGCATTTTGGAAGACTCATGGAATCTTTATCTCTATCTTTCATGTAATGCTTTAATTTTACTTTGTATGAACCCACCTCTCTTCAGGAGTTAAATTAATGAAAGCTTCCTTTGCCCAAAAAAAAAAAAAAAAAAAAAACAAAAATATAT")
#' orfs <- orfFinder(seqs, 50, echo=F)
#' print(orfs)
#'
#'
orfFinder <- function(seqs, min.len = 30, max.len = Inf, echo=F){
  
  codon.table <- c(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  
  if(is.null(names(seqs))){
    if(echo) pokazAttention('New names will be defined')
    names(seqs) = paste('Seq_', 1:length(seqs), '|', nchar(seqs), sep = '')
  }
  
  seqs.aa = c()
  if(echo) pokazAttention('minimum length of AA sequence is', min.len, ', maximum', max.len)
  for(i.s in 1:length(seqs)){
    if (i.s %% 100 == 0) print(i.s)
    s = seq2nt(seqs[i.s])
    
    for(i.orf in 0:2){
      
      # Remove first nts
      s.tmp = s
      if(i.orf != 0){
        s.tmp = s.tmp[-(1:(i.orf))]
      }
      
      # Length should fit mod 3
      remove.nt <- length(s.tmp) %% 3
      s.use <- s.tmp[1:(length(s.tmp) - remove.nt)]
      
      s.strand = '+'
      n = length(s.use)
      for(i.tmp in 1:2){
        z = matrix(s.use, nrow = 3)
        z = apply(z, 2, function(x) paste0(x, collapse = ''))
        a = codon.table[z]
        names(a) = NULL
        
        # Find potisions of Stop-codons
        idx = which(a == '*')
        
        # Lengths of seqeunces between stop-codons
        len.seq = length(a)
        if(!(len.seq %in% idx)) {
          idx = c(idx, len.seq + 1)
        }
        idx.len = c(idx, len.seq + 1) -  c(0, idx) - 1
        
        # Get ORFs with the appropriate lengths
        idx.remain = which((idx.len >= min.len) & (idx.len <= max.len))
        if(length(idx.remain) != 0){
          # Get sequences
          idx =  c(0, idx)
          for(i in idx.remain){
            pos1 = (idx[i] + 1)
            pos2 = (idx[i+1]-1)
            s.tmp  = nt2seq(a[pos1:pos2])
            if(s.strand == '+'){
              s.tmp.name = paste(names(seqs)[i.s], 'pept', (pos1 - 1) * 3 + i.orf + 1, (pos2 + 1) * 3 + i.orf, i.orf, s.strand, pos2 - pos1 + 1, sep = '|')
            } else {
              s.tmp.name = paste(names(seqs)[i.s], 'pept', n - (pos2 + 1) * 3 + 1 + i.orf, n - (pos1) * 3 + 1 + i.orf + 2, i.orf, s.strand, pos2 - pos1 + 1, sep = '|')
            }
            # print(s.tmp.name)
            seqs.aa[s.tmp.name] = s.tmp
          }
        }
        
        # Reverse complement
        s.strand = '-'
        s.use = revCompl(s.use)
      }
      
      
    }
    
  }
  
  return(seqs.aa)
  
}