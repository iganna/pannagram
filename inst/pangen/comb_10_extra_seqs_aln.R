# Get positiona for an extra alignment

suppressMessages({
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
  library(rhdf5)
  # library(muscle)  #BiocManager::install("muscle")
  library(pannagram)
  library(igraph)
  library(R.utils)
  # library(Biostrings)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
source(system.file("pangen/comb_func_extra2.R", package = "pannagram"))
source(system.file("pangen/synteny_func.R", package = "pannagram"))
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

aln.type.in = aln.type.add

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores = opt$cores
if(is.null(num.cores)) stop('Whong number of cores: NULL')

pokaz('Number of cores', num.cores)

# Path with the consensus output
if (!is.null(opt$path.cons)) path.cons <- opt$path.cons
if(!dir.exists(path.cons)) stop('Consensus folder does not exist')

# Path to chromosomes
if (!is.null(opt$path.chromosomes)) path.chromosomes <- opt$path.chromosomes

# Path to mafft input
if (!is.null(opt$path.extra)) path.extra <- opt$path.extra

path.work = paste0(path.extra, 'tmp/')
if (!dir.exists(path.work)) {
  dir.create(path.work)
}

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

echo = T
for(s.comb in pref.combinations){
  
  if(echo) pokaz('* Combination', s.comb)
  
  # Load breaks
  file.breaks.info = paste0(path.extra, "breaks_info_",s.comb,".RData")
  load(file.breaks.info)  # "breaks.init", "breaks"
  
  breaks$len.true = breaks$idx.end - breaks$idx.beg - 1
  
  # Accessions
  file.comb = paste0(path.cons, aln.type.in, s.comb,'.h5')
  groups = h5ls(file.comb)
  accessions = groups$name[groups$group == gr.accs.b]
  
  # ---- Additional alignments ----
  
  
  loop.function <- function(i.b, breaks, echo.loop=T){
  # for (i.b in 1:nrow(breaks)) {
  # for (i.b in 1:nrow(breaks)) {
    
    # ---- Logs ----
    
    # Log files
    file.log.loop = paste0(path.log, 'loop_file_', 
                           i.b,
                           '.log')
    if(!file.exists(file.log.loop)){
      invisible(file.create(file.log.loop))
    }
    
    if(checkDone(file.log.loop)){
      return()
    }
    
    # ---- Main code ----
    
    file.br.group = paste0(path.extra, breaks$id.s[i.b], '_group.fasta')
    file.br.idx = paste0(path.extra, breaks$id.s[i.b], '_group.txt')
    file.br.add = paste0(path.extra, breaks$id.s[i.b], '_add.fasta')
    file.br.out = paste0(path.extra, breaks$id.s[i.b], '_out.RData')
    file.br.len = paste0(path.extra, breaks$id.s[i.b], '_len.RData')
    
    pos.b <- breaks$idx.beg[i.b]
    pos.e <- breaks$idx.end[i.b]
    
    # Select initial breaks that fit within the current merged break
    breaks.tmp <- breaks.init[(breaks.init$idx.beg >= pos.b) & (breaks.init$idx.end <= pos.e),]
    
    # Sort by length (starting from the largest)
    breaks.tmp <- breaks.tmp[order(-breaks.tmp$len.acc),]
    idx.small = which((breaks.tmp$len.acc < 100)&(breaks.tmp$len.comb < 100))
    idx.reordered = c(idx.small, setdiff(1:nrow(breaks.tmp), idx.small))
    breaks.tmp = breaks.tmp[idx.reordered,]
    rownames(breaks.tmp) = NULL
    
    head(breaks.tmp[,-ncol(breaks.tmp),drop=F])
    
    # Read the group idxs
    # pokaz(file.br.idx)
    if(file.exists(file.br.idx)){
      indexes = read.table(file.br.idx, row.names = 1, header = F)
      
      # Read the group alignment
      aln = readFasta(file.br.group)
      aln = aln2mx(aln)
      
      aln.non.gap = colSums(aln != '-')
      if(min(aln.non.gap) == 0) stop('Alignment is not full')
      
      mx.cons = mx2cons(aln, amount = 3)
      mx.cons = cbind('-', mx.cons)
      mx.cons = cbind(mx.cons, '-')
      mx.cons = mx.cons[rowSums(mx.cons != '-') > 0,,drop=F]
      idx.cons = matrix((pos.b):(pos.e), nrow = 1)
      
      # save(list = ls(), file ="tmp_workspace_extra.RData")
      
      if(length(idx.cons) != ncol(mx.cons)){
        save(list = ls(), file ="tmp_workspace_extra.RData")
        stop('Something is wrong with lengths')
      } 
    } else {
      mx.cons = c()
      mx.cons = cbind('-', mx.cons)
      mx.cons = cbind(mx.cons, '-')
      # mx.cons = mx.cons[rowSums(mx.cons != '-') > 0,,drop=F]
      idx.cons = matrix((pos.b):(pos.e), nrow = 1)
    }
    
    
    for(irow in 1:nrow(breaks.tmp)){
      pokaz("Dim mx.cons", dim(mx.cons))
      # for(irow in 54:nrow(breaks.tmp)){
      pokaz('irow', irow, '/', nrow(breaks.tmp))
      s.b = breaks.tmp$seq[irow]
      
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
        mx.cons <- cbind(mx.cons[, idx1, drop=F], 
                         mx.add, 
                         mx.cons[, idx2, drop=F])
        
        idx.cons = rbind(idx.cons, 0)
        idx.add = matrix(0, nrow = nrow(idx.cons), ncol = n.insert)
        idx.add[nrow(idx.add), ] = pos.acc
        
        idx.cons = cbind(idx.cons[, idx1, drop=F], 
                         idx.add, 
                         idx.cons[, idx2, drop=F])
        
      } else {
        i.add = (i.beg+1):(i.end - 1)
        mx.add = mx.cons[,i.add,drop=F]
        mx.add = mx.add[rowSums(mx.add != '-') > 0, , drop=F]
        
        seq1 = mx2aln(mx.add)
        seq2 = s.b
        
        # pokaz('point0')
        
        ## ---- Mafft ----
        pokaz('Length1', nchar(seq1), 'Length2', nchar(seq2), ". Mafft is running")
        
        if(max(c(nchar(seq1), nchar(seq2))) < 15000){
          
          path.work.tmp.pref = tempfile(pattern = paste0("mafft_", breaks$id.s[i.b]), tmpdir = path.work, fileext = "_")
          if(max(c(nchar(seq1), nchar(seq2))) < 25){
            mafft.res = mafftAdd(seq1[1], seq2, path.work.tmp.pref, n.flank = 0)  
          } else {
            mafft.res = mafftAdd(seq1[1], seq2, path.work.tmp.pref, n.flank = 20)    
          }
          
          pokaz("Mafft done")
          
          result = mafft.res$result
          df = mafft.res$df
          pos.mx = mafft.res$pos.mx
          
          # pokaz('point1')
          
          ## ---- Check all the synteny chunks ----
          mafft.mx = matrix('-', nrow = 2, ncol = ncol(pos.mx))
          mafft.mx[1, pos.mx[1,] != 0] = toupper(seq2nt(seq1[1]))
          mafft.mx[2, pos.mx[2,] != 0] = toupper(seq2nt(seq2))
          
          # msaplot(mafft.mx)
          
          
          if(nrow(result) > 0){
            result$score = 0
            
            for(jrow in 1:nrow(result)){
              idx = result$beg[jrow]:result$end[jrow]
              result$score[jrow] = sum(mafft.mx[1,idx] == mafft.mx[2,idx]) / result$len[jrow]
            }
            cnunk.min.len = 25
            len.cutoff = min(c(nchar(seq1[1])/3, nchar(seq2)/3, cnunk.min.len))
            sim.score = 0.9
            irow_support = which((result$score >= sim.score) & (result$len >= len.cutoff)) 
            
            result = result[,-which(colnames(result) == 'score'), drop=F]
            
          } else {
            irow_support = c()
          }
          
          
          # pokaz('point2')
          
          # save(list = ls(), file ="tmp_workspace_point2.RData")
          
          ## ---- BLAST----
          x = data.frame(tmp=numeric())   
          if((length(irow_support) != nrow(result)) & (nrow(result) > 0) & (length(irow_support) > 0)){
            path.work.tmp.pref = tempfile(pattern = paste0("blast_", breaks$id.s[i.b]), tmpdir = path.work, fileext = "_")
            x = blastTwoSeqs(seq1[1], seq2, path.work.tmp.pref)  
          }
          # If some BLAST
          if(nrow(x) > 0){
            
            x$idx = 1:nrow(x)
            x$dir = (x$V4 > x$V5) * 1
            x = x[x$dir == 0,]   # TODO: keep inversions
            # x = glueZero(x)
            
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
              stop("SOMETHING IS WRONG")
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
                save(list = ls(), file ="tmp_workspace_ws1.RData")
                stop('WRONG SORTING')
              }
              
              if(length(pp) != length(unique(pp))){
                stop('WRONG UNIQUE')
              }
            }
            
          } else {
            
            # save(list = ls(), file ="tmp_workspace_s1.RData")
            
            n1 = nchar(seq1)[1]
            n2 = nchar(seq2)[1]
            
            # pokaz(n1, n2)
            
            mx.comb = matrix(0, nrow = 2, ncol = n1 + n2)
            mx.comb[1, 1:n1] = 1:n1
            mx.comb[2, n1 + (1:n2)] = 1:n2
          }
          
        }else {
          
          ## ---- LONG: BLAST only ----
          
          path.work.tmp.pref = tempfile(pattern = paste0("blast_", breaks$id.s[i.b]), tmpdir = path.work, fileext = "_")
          x = blastTwoSeqs(seq1[1], seq2, path.work.tmp.pref) 
          x = x[x$V4 < x$V5 , , drop = F]
          x = x[x$V2 < x$V3 , , drop = F]  
          
          if(nrow(x) != 0){
            
            b = x[,2:5, drop = F]
            x.min = min(c(b$V2, b$V3))
            x.max = max(c(b$V2, b$V3))
            y.min = min(c(b$V4, b$V5))
            y.max = max(c(b$V4, b$V5))
            
            b = rbind(c(x.min, x.min, y.min, y.min),
                      b,
                      c(x.max, x.max, y.max, y.max))
            
            # plotSynDot(b)
            
            d = computeDistBwBlastHits(b)
            
            # Convert the matrix into a graph
            graph <- graph_from_adjacency_matrix(d$dist, mode = "directed", weighted = TRUE)  # from igraph package
            
            # Shortest path
            shortest_path <- shortest_paths(graph, from = 1, to = nrow(b), weights = E(graph)$weight)
            path.idx = shortest_path$vpath[[1]] 
            path.idx = d$order[path.idx]
            path.idx = setdiff(path.idx, c(1, nrow(b)))
            path.idx = path.idx - 1
            
            # plotSynDot(x)
            # plotSynDot(x[path.idx,,drop=F])
            df.sup = x[path.idx,,drop=F]
            
            if(nrow(df.sup) > 1){
              df.sup$dir = 0
              df.sup = cleanOverlaps(df.sup)  
              
              if(sum(nchar(df.sup$V8) != nchar(df.sup$V9)) > 0) stop('Wrong length of sequences')
              
              # Checkup
              for(irow in 1:nrow(df.sup)){
                
                s1 = seq2nt(df.sup$V8[irow])
                s2 = seq2nt(df.sup$V9[irow])
                
                if(sum(s1 != '-') != (df.sup$V3[irow] - df.sup$V2[irow] + 1)){
                  save(list = ls(), file ="tmp_workspace_wl1.RData")
                  stop('Wrong len 1')
                } 
                if(sum(s2 != '-') != (df.sup$V5[irow] - df.sup$V4[irow] + 1)) stop('Wrong len 2')
              }
            }
            
            
            
            # ---- mx.comb ----
            n.sup = nrow(df.sup)
            result.sup = c()
            # ---- Add the gaps from the first sequence ----
            len1 = nchar(seq1)[1]
            result.sup.tmp = data.frame(beg = c(1, df.sup$V3+1),
                                        end = c(df.sup$V2-1, len1),
                                        id = (0:(n.sup)) + 0.3,
                                        # id = (0:(n.sup)),
                                        type=1)
            
            # Lengths of all blocks
            result.sup.tmp$len = result.sup.tmp$end - result.sup.tmp$beg + 1
            result.sup.tmp = result.sup.tmp[result.sup.tmp$len > 0,,drop = F]
            
            # Save
            result.sup = rbind(result.sup, result.sup.tmp)
            
            # ---- Add the gaps from the second sequence ----
            len2 = nchar(seq2)[1]
            result.sup.tmp =  data.frame(beg = c(1, df.sup$V5+1),
                                         end = c(df.sup$V4-1, len2),
                                         id = (0:(n.sup)) + 0.5,
                                         # id = (0:(n.sup)),
                                         type=2)
            
            # Lengths of all blocks
            result.sup.tmp$len = result.sup.tmp$end - result.sup.tmp$beg + 1
            result.sup.tmp = result.sup.tmp[result.sup.tmp$len > 0,,drop = F]
            
            # Save
            result.sup = rbind(result.sup, result.sup.tmp)
            
            # --- Add fake ----
            
            result.sup = rbind(result.sup,
                               data.frame(beg = -1, end = -1, id = 1:n.sup, type = 0, len = nchar(df.sup$V8)))
            
            # ---- Order ----
            result.sup = result.sup[order(result.sup$id),]
            
            result.sup$a.beg = cumsum(c(0, result.sup$len[-nrow(result.sup)])) + 1
            result.sup$a.end = result.sup$a.beg + result.sup$len - 1
            
            if(sum((result.sup$a.end - result.sup$a.beg + 1) != result.sup$len) > 0) {
              stop("SOMETHING IS WRONG")
            }
            
            mx.comb = matrix(0, nrow = 2, ncol = max(result.sup$a.end))
            for(irow in 1:nrow(result.sup)){
              if(result.sup$type[irow] == 0){
                
                df.irow = df.sup[result.sup$id[irow],]
                
                pos = result.sup$a.beg[irow]:result.sup$a.end[irow]
                pos1 = df.irow$V2:df.irow$V3
                pos2 = df.irow$V4:df.irow$V5
                
                s1 = seq2nt(df.irow$V8)
                s2 = seq2nt(df.irow$V9)
                
                mx.comb[1,pos[s1 != '-']] = pos1
                mx.comb[2,pos[s2 != '-']] = pos2
              } else {
                pos = result.sup$a.beg[irow]:result.sup$a.end[irow]
                posi = result.sup$beg[irow]:result.sup$end[irow]
                mx.comb[result.sup$type[irow],pos] = posi
                
              }
            }
            
            
            for(irow in 1:2){
              pp = mx.comb[irow,]
              pp = pp[pp != 0]
              if(is.unsorted(pp)){
                save(list = ls(), file ="tmp_workspace_ws2.RData")
                stop('WRONG SORTING')
              }
              
              if(length(pp) != length(unique(pp))){
                stop('WRONG UNIQUE')
              }
            }
            
          } else {
            n1 = nchar(seq1)[1]
            n2 = nchar(seq2)[1]
            
            # pokaz(n1, n2)
            
            mx.comb = matrix(0, nrow = 2, ncol = n1 + n2)
            mx.comb[1, 1:n1] = 1:n1
            mx.comb[2, n1 + (1:n2)] = 1:n2
          }
          
        }
        
        # pokaz('point3')
        
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
        mx.cons = cbind(mx.cons[,idx1,drop=F],
                        mx.cons.new,
                        mx.cons[,idx2,drop=F])
        
        s.cons = mx2cons(mx.cons)
        mx.cons[1,] = s.cons
        mx.cons = mx.cons[1:nrow(mx.cons)-1,,drop=F]
        
        idx.new = matrix(0, nrow = (nrow(idx.cons) + 1), ncol = ncol(mx.cons.new))
        idx.new[1:nrow(idx.cons), non.zero.indices.1] = idx.cons[,i.add]
        idx.new[nrow(idx.new),non.zero.indices.2] = pos.acc
        
        idx.cons = rbind(idx.cons, 0)
        idx.cons = cbind(idx.cons[, idx1,drop=F], 
                         idx.new, 
                         idx.cons[, idx2,drop=F])
        
        if(ncol(idx.cons) != ncol(mx.cons)) stop('Checkpoint lengths')
        
        
      }
    }
    
    if(nrow(idx.cons) != (nrow(breaks.tmp) + 1)) stop('Wrong number of rows')
    
    
    # ---- Combine the alignments ---- 
    ## ---- Initialisation ----
    msa.new = matrix('-', nrow = length(accessions), ncol = ncol(mx.cons),
                     dimnames = list(accessions, NULL))
    
    idx.new = matrix(0, nrow = length(accessions), ncol = ncol(mx.cons),
                     dimnames = list(accessions, NULL))
    
    ## ---- Add previous from initial alignment ----
    # Sequences
    
    pokaz(file.br.idx)
    if(file.exists(file.br.idx)){
      idx.aln = which(idx.cons[1,] != 0) 
      idx.aln = idx.aln[-c(1, length(idx.aln))]
      accs = sapply(rownames(aln), function(s) strsplit(s, '_break_')[[1]][1])
      accs <- sub("^acc_", "", accs)
      msa.new[accs, idx.aln] = aln
      # Positions
      idx.new[accs, idx.aln] <- as.matrix(indexes)
    }

    ## ---- Add previous from new alignment ----
    for(irow in 1:nrow(breaks.tmp)){
      if(sum(msa.new[breaks.tmp$acc[irow], idx.cons[irow + 1,] != 0] != '-') > 0) stop('wrong with accessions msa')
      if(sum(idx.new[breaks.tmp$acc[irow], idx.cons[irow + 1,] != 0] != 0) > 0) stop('wrong with accessions idx')
      
      msa.new[breaks.tmp$acc[irow], idx.cons[irow + 1,] != 0] = seq2nt(breaks.tmp$seq[irow])
      pos.irow = breaks.tmp$val.beg[irow]:breaks.tmp$val.end[irow]
      pos.irow = pos.irow[-1]
      pos.irow = pos.irow[-length(pos.irow)]
      idx.new[breaks.tmp$acc[irow], idx.cons[irow + 1,] != 0] = pos.irow
    }
    
    # msaplot(aln)
    # msaplot(msa.new)
    # 
    # msadiff(msa.new)
    
    # Empty columns
    idx.zero = which(colSums(idx.new != 0) == 0)
    idx.gap = which(colSums(msa.new != '-') == 0)
    if(length(idx.zero) > 0){

      if(!setequal(idx.zero, idx.gap)) stop('Wrong zeros')
      idx.new = idx.new[,-idx.zero, drop = F]
      msa.new = msa.new[,-idx.gap, drop = F]
    }
    
    # Length
    len.idx = ncol(idx.new)
    len.msa = ncol(msa.new)
    if(len.idx != len.msa) stop('Lwngth are different')
    len.aln = len.idx
    
    
    if(len.aln < breaks$len.true[i.b]) stop(paste('The length has decreased', i.b))
    
    # Save
    save(list = c("len.aln"), file = file.br.len)
    save(list = c("idx.new", "msa.new"), file = file.br.out)
    
    pokaz('Done.', file=file.log.loop, echo=echo.loop)
    return(NULL)
  }
  
  max.tile.loop = 600
  if(num.cores == 1){
    for(i.b in 1:nrow(breaks)){
      tryCatch({
        # Set a timeout of 600 seconds (10 minutes) for the function execution
        withTimeout({
          loop.function(i.b, breaks, echo.loop = echo.loop)
        }, timeout = max.tile.loop)
      }, TimeoutException = function(ex) {
        pokazAttention("Timeout reached for i.b = %d. Skipping to the next iteration.", i.b)
      })
    }
  } else {
    # Set the number of cores for parallel processing
    myCluster <- makeCluster(num.cores, type = "PSOCK") 
    registerDoParallel(myCluster) 
    
    tmp = foreach(i.b = 1:nrow(breaks), 
                  .packages=c('rhdf5', 'crayon', 'igraph', 'pannagram', 'R.utils'))  %dopar% { 
                    tryCatch({
                      # Set a timeout of 600 seconds (10 minutes) for each parallel task
                      withTimeout({
                        loop.function(i.b, breaks, echo.loop = echo.loop)
                      }, timeout = max.tile.loop)
                    }, TimeoutException = function(ex) {
                      pokazAttention("Timeout reached for i.b = %d. Skipping to the next iteration.", i.b)
                      return(NULL) # Return NULL for skipped tasks
                    })
                  }
    stopCluster(myCluster)
  }
  
  
  
  H5close()
  gc()
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

