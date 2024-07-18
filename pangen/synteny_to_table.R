# GAP, SNP, DIV


source("utils/utils.R")
source('visualisation/dotplot.R')



path.aln = ''
pref.comb = ''

file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')


file.aln.pre = '/Volumes/Samsung_T5/vienn/pannagram_test/symA_test2/alignments_ml2/ml4_1_1_full.rds'

x = readRDS(file.aln.pre)

head(x)

# Dataframe to produce
# name.q beg.q end.q len.q name.r beg.r end.r len.r
x.tab <- data.frame(
  name.q = character(),
  beg.q = numeric(),
  end.q = numeric(),
  len.q = numeric(),
  name.r = character(),
  beg.r = numeric(),
  end.r = numeric(),
  len.r = numeric(),
  stringsAsFactors = FALSE
)


diff.tab <- data.frame(
  beg.q = numeric(),
  end.q = numeric(),
  beg.r = numeric(),
  end.r = numeric(),
  type = character(),
  dir = numeric(),
  stringsAsFactors = FALSE
)

# If it is GAP mode, the positions beg and end will be flanking positions (which actually correspond to each other in the alignment).
# If it is SNP mode, the positions will and correspond to the SNP position.
# If it is HDR mode, then again, as in the case with GAP, there will be flanking aligned positions.

for(irow in 1:nrow(x)){
  
  # Sequences as vectors of nucleotides
  s.q = seq2nt(x$V8[irow])
  s.r = seq2nt(x$V9[irow])
  aln.len = length(s.q)
  
  # Positions in alignment
  p.q = rep(0, aln.len)
  p.r = rep(0, aln.len)
  
  p.q[s.q != '-'] = x$V2[irow]:x$V3[irow]
  p.r[s.r != '-'] = x$V4[irow]:x$V5[irow]
  
  # Gaps
  idx.gap = (s.q == '-') | (s.r == '-')
  pos.gap = findOnes(idx.gap * 1)
  # Dataframe to save
  if(nrow(pos.gap) != 0) {
    # stop('gap')
    diff.tmp = data.frame(beg.q = p.q[pos.gap$beg-1], 
                          end.q = p.q[pos.gap$end+1],
                          beg.r = p.r[pos.gap$beg-1], 
                          end.r = p.r[pos.gap$end+1],
                          type = 'GAP',
                          dir = x$dir[irow])
    
    diff.tab = rbind(diff.tab, diff.tmp)
  }
  
  # SNPs in the regious without gaps
  idx.snp = (s.q != s.r) & !(idx.gap)
  pos.snp = findOnes(idx.snp * 1)
  # Dataframe to save
  if(nrow(pos.snp) != 0) {
    # stop('snp')
    diff.tmp = data.frame(beg.q = p.q[pos.snp$beg], 
                          end.q = p.q[pos.snp$end],
                          beg.r = p.r[pos.snp$beg], 
                          end.r = p.r[pos.snp$end],
                          type = 'SNP',
                          dir = x$dir[irow])
    diff.tab = rbind(diff.tab, diff.tmp)
  }
  rm(diff.tmp)
}

# Between synteny blocks
n.row = nrow(x)
idx.diff = which( ( ((abs(x$V2[-1] - x$V3[-n.row])) != 1) | 
                    ((abs(x$V4[-1] - x$V5[-n.row])) != 1) )  &
                  ( x$dir[-1] == x$dir[-n.row] ) )
diff.tmp = data.frame(beg.q = x$V3[idx.diff], 
                      end.q = x$V2[idx.diff+1],
                      beg.r = x$V5[idx.diff], 
                      end.r = x$V4[idx.diff+1],
                      type = 'HDR',
                      dir = x$dir[idx.diff])

diff.tab = rbind(diff.tab, diff.tmp)
rm(diff.tmp)

# Flanking for inversions


# Inversions needed?


# Sorting
diff.tab = diff.tab[order(diff.tab$beg.q),]
diff.tab$len.q = abs(diff.tab$end.q - diff.tab$beg.q) - 1
diff.tab$len.r = abs(diff.tab$end.r - diff.tab$beg.r) - 1
diff.tab$len.q[diff.tab$type == 'SNP'] = 0
diff.tab$len.r[diff.tab$type == 'SNP'] = 0

diff.tab$type[(diff.tab$type == 'HDR') &(diff.tab$len.q*diff.tab$len.r == 0) & (diff.tab$len.q+diff.tab$len.r != 0)] = 'GAP+'


# ***********************************************************************
# ---- Manual testing ----

if(F){
  library(ggplot2)
  source("../utils/utils.R")
  source('../visualisation/dotplot.R')
  source('../visualisation/orfplot.R')
}


dotplot(s1, s2, 15, 13)


