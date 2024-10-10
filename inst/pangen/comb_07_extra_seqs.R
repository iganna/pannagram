# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  library(muscle)  #BiocManager::install("muscle")
  # library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
# source("synteny_funcs.R")

# pokazStage('Step 10. Prepare sequences for MAFFT')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--path.cons"), type="character", default=NULL, 
              help="path to consensus directory", metavar="character"),
  make_option(c("--path.chromosomes"), type="character", default=NULL, 
              help="path to directory with chromosomes", metavar="character"),
  make_option(c("--path.extra"), type="character", default=NULL, 
              help="path to directory, where to combine fasta files for mafft runs", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer"),
  make_option(c("--path.log"), type = "character", default = NULL,
              help = "Path for log files", metavar = "character"),
  make_option(c("--log.level"), type = "character", default = NULL,
              help = "Level of log to be shown on the screen", metavar = "character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files

# ***********************************************************************

#TODO: SHOULD BE PARAMATERS ?

len.short = 50
# len.large = 40000
n.flank = 30

aln.type.in = aln.type.pre
aln.type.out = aln.type.extra

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Whong number of cores: NULL')

pokaz('Number of cores', num.cores)
if(num.cores > 1){
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
}
# num.cores.max = 10
# num.cores <- min(num.cores.max, ifelse(!is.null(opt$cores), opt$cores, num.cores.max))
# if(num.cores > 1){
#   myCluster <- makeCluster(num.cores, type = "PSOCK") 
#   registerDoParallel(myCluster)
# }

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder doesn???t exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# Path to mafft input
if (!is.null(opt$path.extra)) path.extra <- opt$path.extra

# **************************************************************************
# ---- Combinations of chromosomes query-base to create the alignments ----

s.pattern <- paste0("^", aln.type.in, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.type.in, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

if(length(pref.combinations) == 0) {
  stop('No files with the ref-based alignments are found')
}

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)

# ***********************************************************************
# ---- MAIN program body ----

echo = F
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  q.chr = strsplit(s.comb, '_')[[1]][1]
  
  file.comb = paste0(path.cons, aln.type.in, s.comb,'.h5')
  
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  
  # Get breaks
  breaks.init = c()
  for(acc in accessions){
    pokaz(acc)
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb, s.acc)
    breaks.acc = findBreaks(v)
    breaks.init = rbind(breaks.init, breaks.acc)
  }
  
  # Fix the borders
  breaks.init$idx.beg = breaks.init$idx.beg
  breaks.init$idx.end = breaks.init$idx.end
  
  breaks.init = breaks.init[order(-breaks.init$idx.end),]
  breaks.init = breaks.init[order(breaks.init$idx.beg),]
  breaks.init$id = 1:nrow(breaks.init)
  
  # Breaks from the previous data
  
  breaks <- mergeOverlaps(breaks.init)
  n.digits <- nchar(as.character(nrow(breaks)))
  format.digits <- paste0("%0", n.digits, "d")
  
  # Get the consensus sequences ans save then into files
  
  seqs.add = c()
  for(acc in accessions){
    
    # Read the chromosome
    file.chromosome = paste(path.chromosomes, 
                            acc, 
                            '_chr', q.chr, '.fasta', sep = '')
    genome = readFasta(file.chromosome)
    genome = seq2nt(genome)
    
    # Read the alignment
    s.acc = paste0(gr.accs.e, acc)
    v = h5read(file.comb, s.acc)
    
    for(i.b in 1:nrow(breaks)){
      
      pos.b <- breaks$idx.beg[i.b] + 1
      pos.e <- breaks$idx.end[i.b] - 1
      
      v.b = v[pos.b:pos.e]
      s.b = rep('-', length(v.b))
      s.b[v.b != 0] = genome[abs(v.b)]
      idx.rc = which(v.b < 0)
      if(length(idx.rc) > 0){
        s.b[idx.rc] = justCompl(s.b[idx.rc])
      }
      s.b = nt2seq(s.b)
      s.b = setNames(s.b, 
                     paste0('acc_', acc, '_br_', i.b))
      
      file.br.fasta = paste0(path.extra, 'break_', sprintf(format.digits, i.b), '.fasta')
      if(!file.exists(file.br.fasta)){
        writeFasta(s.b, file.br.fasta)  
      } else {
        writeFasta(s.b, file.br.fasta, append=T)
      }
      
    }
    
    for(i.b in which(breaks.init$acc == acc)){
      
      pos.b <- breaks.init$val.beg[i.b] + 1
      pos.e <- breaks.init$val.end[i.b] - 1
      
      poses = pos.b:pos.e
      s.b = genome[abs(poses)]
      idx.rc = which(poses < 0)
      if(length(idx.rc) > 0){
        s.b[idx.rc] = justCompl(s.b[idx.rc])
      }
      s.b = nt2seq(s.b)
      s.b.name = paste0('break_id_', breaks.init$id[i.b])
      s.b = setNames(s.b, s.b.name)
      
      seqs.add = c(seqs.add, s.b)
    }
    
  }
  
  if(length(seqs.add) != nrow(breaks.init)) stop('Not all of the sequences were created')

  # ---- Additional alignments ----
  
  for (i.b in 1:nrow(breaks)) {

    pos.b <- breaks$idx.beg[i.b]
    pos.e <- breaks$idx.end[i.b]

    # Select initial breaks that fit within the current merged break
    breaks.tmp <- breaks.init[(breaks.init$idx.beg >= pos.b) & (breaks.init$idx.end <= pos.e),]

    # Sort by length (starting from the largest)
    breaks.tmp <- breaks.tmp[order(breaks.tmp$len.comb),]

    file.br.fasta = paste0(path.extra, 'break_', sprintf(format.digits, i.b), '.fasta')
    aln = readFasta(file.br.fasta)
    aln = aln2mx(aln)

    mx.cons = mx2cons(aln, amount = 3)
    mx.cons = cbind('-', mx.cons)
    mx.cons = cbind(mx.cons, '-')
    idx.cons = matrix((pos.b):(pos.e), nrow = 1)

    if(length(idx.cons) != ncol(mx.cons)) stop('Something is wrong with lengths')

    for(irow in 1:nrow(breaks.tmp)){
      pokaz(irow)
      s.b.name = paste0('break_id_', breaks.tmp$id[irow])
      s.b = seqs.add[s.b.name]
      
      pos.acc = (breaks.tmp$val.beg[irow]+1):(breaks.tmp$val.end[irow]-1)
      
      if(nchar(s.b) != length(pos.acc)) stop('Mismatch sequence and positions')
      
      
      i.beg = which(idx.cons[1,] == breaks.tmp$idx.beg[irow])
      i.end = which(idx.cons[1,] == breaks.tmp$idx.end[irow])
      i.n = i.end - i.beg - 1
      if(i.n == 0){
        idx1 = 1:i.beg
        idx2 = i.end:ncol(mx.cons)
        n.insert = breaks.tmp$len.acc[irow]
        
        mx.add = matrix('-', nrow = nrow(mx.cons), ncol = n.insert)
        mx.add[1,] = seq2nt(s.b)
        mx.cons <- cbind(mx.cons[, idx1], 
                         mx.add, 
                         mx.cons[, idx2])
        
        idx.cons = rbind(idx.cons, 0)
        idx.add = matrix(0, nrow = nrow(idx.cons), ncol = n.insert)
        idx.add[nrow(idx.add), ] = pos.acc
        
        idx.cons = cbind(idx.cons[, idx1], 
                         idx.add, 
                         idx.cons[, idx2])
        
      } else {
        i.add = (i.beg+1):(i.end - 1)
        mx.add = mx.cons[,i.add,drop=F]
        mx.add = mx.add[rowSums(mx.add != '-') > 0, , drop=F]
        
        seq1 = mx2aln(mx.add)
        seq2 = s.b
        
        path.work = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/extra_regions/tmp/'
        mafft.res = mafftAdd(seq1, seq2, path.work)
        
        result = mafft.res$result
        df = mafft.res$df
        pos.mx = mafft.res$pos.mx
        
        ## ---- Check all the synteny chunks ----
        mafft.mx = matrix('-', nrow = 2, ncol = ncol(pos.mx))
        mafft.mx[1, pos.mx[1,] != 0] = toupper(seq2nt(seq1[1]))
        mafft.mx[2, pos.mx[2,] != 0] = toupper(seq2nt(seq2))
        
        # msaplot(mafft.mx)
        
        irow_support = c()
        cnunk.min.len = 25
        sim.score = 0.9
        for(irow in 1:nrow(result)){
          if(result$len[irow] < min(c(nchar(seq1[1])/3, nchar(seq2)/3, cnunk.min.len))) next
          idx = result$beg[irow]:result$end[irow]
          score = sum(mafft.mx[1,idx] == mafft.mx[2,idx]) / result$len[irow]
          if(score >= sim.score){
            irow_support = c(irow_support, irow)
          }
        }
        
        ## ---- BLAST----
        x = data.frame(tmp=numeric())   
        if(length(irow_support) != nrow(result)){
          x = blastTwoSeqs(seq1[1], seq2, path.work)  
        }
        # If some BLAST
        if(nrow(x) > 0){
          
          x$idx = 1:nrow(x)
          x$dir = (x$V4 > x$V5) * 1
          x = x[x$dir == 0,]   # TODO: keep inversions
          x = glueZero(x)
          
          # test df lines
          for (irow in setdiff(1:nrow(df), irow_support)) {
            
            start_df <- df$V2[irow]
            end_df <- df$V3[irow]
            
            length_df <- end_df - start_df
            
            overlap_start <- pmax(start_df, x$V2) 
            overlap_end <- pmin(end_df, x$V3)     
            
            overlap_length <- overlap_end - overlap_start
            
            valid_overlap <- overlap_length >= 0.9 * length_df & overlap_length > 0
            
            if(sum(valid_overlap) == 0) next
            
            df.tmp = x[valid_overlap,]
            
            df.tmp[,c('V4', 'V5')] = df.tmp[,c('V4', 'V5')] - df$V4[irow] + df$V2[irow]
            
            valid_overlap2 = pmin(df.tmp$V3, df.tmp$V5) - pmax(df.tmp$V2, df.tmp$V4)
            valid_overlap2 = valid_overlap2 / df.tmp$V7
            
            if(sum(valid_overlap2 > 0.9) > 0) {
              irow_support = c(irow_support, irow)
            }
            
          }
        }
        #
        if(length(irow_support) != 0){
          irow_support = sort(irow_support)
          
          n.sup = length(irow_support)
          
          result.sup = result[irow_support,,drop=F]
          result.sup$type = 0
          result.sup$id = 1:n.sup
          result.sup$len = result.sup$end - result.sup$beg + 1
          
          df.sup = df[irow_support,,drop=F]
          
          # stop()
          
          # ---- Add the gaps from the first sequence ----
          len1 = sum(pos.mx[1,] != 0)
          result.sup.tmp = data.frame(beg = c(1, df.sup$V3+1),
                                      end = c(df.sup$V2-1, len1),
                                      # id = (0:(n.sup)) + 0.3),
                                      id = (0:(n.sup)),
                                      type=1)
          
          # Lengths of all blocks
          result.sup.tmp$len = result.sup.tmp$end - result.sup.tmp$beg + 1
          result.sup.tmp = result.sup.tmp[result.sup.tmp$len > 0,,drop = F]
          
          # define ID by the order in the initial alignment
          if(nrow(result.sup.tmp) > 0){
            len.mafft.aln = ncol(pos.mx)
            for(irow in 1:nrow(result.sup.tmp)){
              tmp = which((result.sup.tmp$beg[irow] <= pos.mx[1,]) & (pos.mx[1,] <= result.sup.tmp$end[irow]))
              result.sup.tmp$id[irow] = result.sup.tmp$id[irow] + (mean(tmp)) / len.mafft.aln
            }
            result.sup = rbind(result.sup, result.sup.tmp)        
          }
          
          
          # ---- Add the gaps from the second sequence ----
          len2 = sum(pos.mx[2,] != 0)
          result.sup.tmp =  data.frame(beg = c(1, df.sup$V5+1),
                                       end = c(df.sup$V4-1, len2),
                                       # id = (0:(n.sup)) + 0.5),
                                       id = (0:(n.sup)),
                                       type=2)
          
          # Lengths of all blocks
          result.sup.tmp$len = result.sup.tmp$end - result.sup.tmp$beg + 1
          result.sup.tmp = result.sup.tmp[result.sup.tmp$len > 0,,drop = F]
          
          # define ID by the order in the initial alignment
          if(nrow(result.sup.tmp) > 0){
            len.mafft.aln = ncol(pos.mx)
            for(irow in 1:nrow(result.sup.tmp)){
              tmp = which((result.sup.tmp$beg[irow] <= pos.mx[2,]) & (pos.mx[2,] <= result.sup.tmp$end[irow]))
              result.sup.tmp$id[irow] = result.sup.tmp$id[irow] + mean(tmp) / len.mafft.aln
            }
            result.sup = rbind(result.sup, result.sup.tmp)
          }
          
          
          # ---- Order ----
          result.sup = result.sup[order(result.sup$id),]
          
          result.sup$a.beg = cumsum(c(0, result.sup$len[-nrow(result.sup)])) + 1
          result.sup$a.end = result.sup$a.beg + result.sup$len - 1
          
          if(sum((result.sup$a.end - result.sup$a.beg + 1) != result.sup$len) > 0) {
            stop("SOMETHING IS WROTG")
          }
          
          # ---- Consensus positions ----
          
          mx.comb = matrix(0, nrow = 2, ncol = result.sup$a.end[nrow(result.sup)])
          for(irow in 1:nrow(result.sup)){
            if(result.sup$type[irow] == 0){
              mx.comb[,result.sup$a.beg[irow]:result.sup$a.end[irow]] = 
                pos.mx[,result.sup$beg[irow]:result.sup$end[irow]]
            } else {
              
              mx.comb[result.sup$type[irow],result.sup$a.beg[irow]:result.sup$a.end[irow]] = 
                result.sup$beg[irow]:result.sup$end[irow]
            }
          }
          
          rowSums(pos.mx != 0)
          
          if(sum(rowSums(mx.comb != 0) !=  rowSums(pos.mx != 0)) > 0){
            stop('MISSED POSITIONS')
            setdiff(pos.mx[1,], mx.comb[1,])
            setdiff(pos.mx[2,], mx.comb[2,])
          } 
          
          # stop()
          
          for(irow in 1:2){
            pp = mx.comb[irow,]
            pp = pp[pp != 0]
            if(is.unsorted(pp)){
              stop('WRONG SORTING')
            }
            
            if(length(pp) != length(unique(pp))){
              stop('WRONG UNIQUE')
            }
          }
          
        } else {
          
          n1 = nchar(s1)
          n2 = nchar(s2)
          mx.comb = matrix(0, nrow = 2, ncol = n1 + n2)
          mx.comb[1, 1:n1] = 1:n1
          mx.comb[2, n1 + (1:n2)] = 1:n2
        }
        
        
        
        non.zero.indices.1 <- mx.comb[1,] != 0
        non.zero.indices.2 <- mx.comb[2,] != 0
        
        mx.cons.new = matrix('-',
                             nrow = (nrow(mx.cons) + 1),
                             ncol = ncol(mx.comb))
        mx.cons.new[1:(nrow(mx.cons)), non.zero.indices.1] = mx.cons[,i.add,drop=F]
        mx.cons.new[nrow(mx.cons.new), non.zero.indices.2] = seq2nt(seq2)
        
        
        # Combine
        idx1 = 1:i.beg
        idx2 = i.end:ncol(mx.cons)
        mx.cons = rbind(mx.cons, '-')
        mx.cons = cbind(mx.cons[,idx1],
                            mx.cons.new,
                        mx.cons[,idx2])
        
        s.cons = mx2cons(mx.cons)
        mx.cons[1,] = s.cons
        mx.cons = mx.cons[1:nrow(mx.cons)-1,]
        
        idx.new = matrix(0, nrow = (nrow(idx.cons) + 1), ncol = ncol(mx.cons.new))
        idx.new[1:nrow(idx.cons), non.zero.indices.1] = idx.cons[,i.add]
        idx.new[nrow(idx.new),non.zero.indices.2] = pos.acc
        
        idx.cons = rbind(idx.cons, 0)
        idx.cons = cbind(idx.cons[, idx1], 
                         idx.new, 
                         idx.cons[, idx2])
        
        if(ncol(idx.cons) != ncol(mx.cons)) stop('Chrckpoint lengths')
        
        
      }
    }


    # msaplot(aln)
    # msaplot(mx.cons)



  }

 
  H5close()
  gc()
}


if(num.cores > 1){
  stopCluster(myCluster)
}

warnings()


# ***********************************************************************
# ---- Manual testing ----

if(F){

  library(rhdf5)
  path.cons = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/consensus/'
  path.mafft.in = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/mafft_in/'
  path.chromosomes = '/Volumes/Samsung_T5/vienn/test/symA_test_0/intermediate/chromosomes/'
  
}




