suppressMessages({
  library(Biostrings)
  library(seqinr)
  #library(spgs)  # reverseComplement("actg")
  library(foreach)
  library(doParallel)
  library(stringr)
  library(optparse)
  library(crayon)
})


source("utils/utils.R")
source("pangen/synteny_funcs.R")

pokazStage('Step 6. Alignment-2. Fill the gaps between synteny blocks')

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-q", "--path.query"), type="character", default=NULL, 
              help="path to query chomosome fasta files", metavar="character"),  
  make_option(c("-o", "--path.aln"), type="character", default=NULL, 
              help="path to the output directory with alignments", metavar="character"),
  make_option(c("-p", "--pref"), type="character", default=NULL, 
              help="prefix of the reference file", metavar="character"),
  make_option(c("-g", "--path.gaps"), type="character", default=NULL, 
              help="prefix of the directory with gaps", metavar="character"),
  make_option(c("-r", "--path.ref"), type="character", default=NULL, 
              help="path to the reference file", metavar="character"),
  make_option(c("-n", "--n.chr.ref"), type="character", default=NULL, 
              help="number of chromosomes in the reference genome", metavar="character"),
  make_option(c("-m", "--n.chr.acc"), type="character", default=NULL, 
              help="number of chromosomes in the accessions", metavar="character"),
  make_option(c("-a", "--all.vs.all"), type="character", default=NULL, 
              help="alignment of all chromosomes vs all or not: T/F", metavar="character"),
  make_option(c("-c", "--cores"), type = "integer", default = 1, 
              help = "number of cores to use for parallel processing", metavar = "integer")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);

# print(opt)

# ***********************************************************************
# ---- Values of parameters ----

# Number of cores
num.cores <- ifelse(!is.null(opt$cores), opt$cores, 30)


if (!is.null(opt$path.query)) path.query <- opt$path.query
if (!is.null(opt$path.aln)) path.aln <- opt$path.aln
if (!is.null(opt$path.ref)) path.base <- opt$path.ref
if (!is.null(opt$pref)) base.acc <- opt$pref
if (!is.null(opt$path.gaps)) path.gaps <- opt$path.gaps

if (!is.null(opt$n.chr.ref)) n.chr.ref <- opt$n.chr.ref
if (!is.null(opt$n.chr.acc)) n.chr.acc <- opt$n.chr.acc
if (!is.null(opt$all.vs.all)) all.vs.all <- as.logical(opt$all.vs.all)


if(!dir.exists(path.aln)) dir.create(path.aln)
if(!dir.exists(path.gaps)) dir.create(path.gaps)


# ***********************************************************************
# ---- Preparation ----

files.maj <- list.files(path.aln, pattern = "\\maj.rds$")
pokaz('Number of alignment files:', length(files.maj))
# if(length(files.maj) == 0) stop('No alignment files provided')



# ***********************************************************************
# ---- MAIN program body ----


loop.function <- function(f.maj, echo = T){
  
  # Remove extensions
  pref.comb <- sub("\\maj.rds$", "", f.maj)
  
  # Parser for BLAST-resulr file
  parts <- strsplit(pref.comb, "_")[[1]]
  
  base.chr <- parts[length(parts) - 1]
  query.chr <- parts[length(parts)]
  parts = parts[-c(length(parts) - 1, length(parts))]
  acc <- paste0(parts, collapse = '_')
  
  pokaz(acc, query.chr, base.chr)
  
  
  pref.comb = paste0(acc, '_', query.chr, '_', base.chr, collapse = '')
  
  # If the previous file was not created - next
  file.aln.pre <- paste(path.aln, paste0(pref.comb, '_maj.rds', collapse = ''), sep = '')
  if(!file.exists(file.aln.pre)) {
    return(NULL)
  }
  
  # If the full file has beed already created - next
  file.aln.full <- paste(path.aln, paste0(pref.comb,  '_full.rds', collapse = ''), sep = '')
  if(file.exists(file.aln.full)) {
    return(NULL)
  }
  
  pokaz('Alignment:', acc, query.chr, base.chr)
  
  # # ---- Read genomes ----
  # 
  # # Read reference sequences
  # base.file = paste0(base.acc, '_chr', base.chr , '.fasta', collapse = '')
  # pokaz('Base:', base.file)
  # base.fas.fw = readFastaMy(paste(path.base, base.file, sep = ''))
  # base.fas.fw = seq2nt(base.fas.fw)
  # base.fas.bw = revCompl(base.fas.fw)
  # base.len = length(base.fas.bw)
  # 
  # # Read query sequences
  # query.file = paste(acc, '_chr',query.chr, '.fasta', sep = '')
  # pokaz('Query:', query.file)
  # 
  # query.fas.chr = readFastaMy(paste(path.query, query.file, sep = ''))
  # query.fas.chr = seq2nt(query.fas.chr)
  # query.len = length(query.fas.chr)
  
  # ---- Read Major ----
  
  pokaz('Read the skeleton alignment..')
  x.sk = readRDS(file.aln.pre)
  
  # ---- Read Gaps bwteen normal blocks ----
  
  complexity.threshold = 200  # Max number of blast hits between two synteny blocks
  
  pos.beg.info = 2  # Position in sequence's name, where the genome's begin is started
  
  file.gaps.out = paste0(path.gaps,
                         'acc_', acc, 
                         '_qchr_', query.chr, '_bchr_', base.chr, '_out.txt', collapse = '')
  
  pokaz('gap file', file.gaps.out)
  
  if.table <- T
  if(file.exists(file.gaps.out)){
    
    pokaz('Read blast of good gaps..')
    
    tryCatch({
      x.gap = read.table(file.gaps.out, stringsAsFactors = F)
    }, error = function(e) {
      if.table <<- F
    })
  } else {
    if.table <- F
  }
  
  if(if.table){
    pokaz('Number of gaps', nrow(x.gap))
    
    x.gap$pref1 = sapply(x.gap$V1, function(s) strsplit(s, '_query')[[1]][1])
    x.gap$pref2 = sapply(x.gap$V10, function(s) strsplit(s, '_base')[[1]][1])
    x.gap = x.gap[x.gap$pref1 == x.gap$pref2,]
    
    x.gap$q.beg = as.numeric(sapply(x.gap$V1, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    x.gap$b.beg = as.numeric(sapply(x.gap$V10, function(s) strsplit(s, '\\|')[[1]][pos.beg.info])) - 1
    x.gap$V2 = x.gap$V2 + x.gap$q.beg
    x.gap$V3 = x.gap$V3 + x.gap$q.beg
    x.gap$V4 = x.gap$V4 + x.gap$b.beg
    x.gap$V5 = x.gap$V5 + x.gap$b.beg
    x.gap$idx = 1:nrow(x.gap)
    x.gap$dir = (x.gap$V4 > x.gap$V5) * 1
    
    cnt = table(x.gap$pref1)
    
    idx.good = c()
    for(i in 1:length(cnt)){
      if(i %% 100 == 0) pokaz('Number of analysed gaps', i)
      
      # name of the node
      s = names(cnt)[i]
      
      # If only one BLAST-hit - get it as it is.
      if(cnt[s] == 1) {
        idx.good = c(idx.good, x.gap$idx[x.gap$pref1 == s])
        next
      }
      
      # Add two extra "nodes", for the begin and end
      x.tmp = x.gap[x.gap$pref1 == s,c('V2', 'V3', 'V4', 'V5', 'idx')]
      
      if(nrow(x.tmp) > complexity.threshold) next
      
      q.beg = min(x.tmp$V2) - 1  # change
      q.end = max(x.tmp$V3) + 1  # change
      b.beg = min(c(x.tmp$V4, x.tmp$V5)) - 1  # change
      b.end = max(c(x.tmp$V4, x.tmp$V5))+ 1  # change
      
      x.tmp = rbind(c(q.beg, q.beg, b.beg, b.beg, 0), x.tmp)
      x.tmp = rbind(x.tmp, c(q.end, q.end, b.end, b.end, 0))
      x.tmp = x.tmp[order(x.tmp$V2),]
      x.tmp$w = abs(x.tmp$V3 - x.tmp$V2) + abs(x.tmp$V4 - x.tmp$V5)
      
      # Set up the direction
      x.tmp0 = x.tmp
      idx.dir = which(x.tmp$V5 < x.tmp$V4)
      tmp = x.tmp$V4[idx.dir]
      x.tmp$V4[idx.dir] = x.tmp$V5[idx.dir]
      x.tmp$V5[idx.dir] = tmp
      
      visit.info = initVisitInfo(nrow(x.tmp))
      
      visit.info2 = graphTraverseWnd(x.tmp, 1, x.tmp$V3[1], x.tmp$V5[1], 0, 0, visit.info)
      
      # Reconstruct the optimal path
      idx.opt = reconstructTraverse(visit.info2$v.prev)
      idx.opt = idx.opt[-c(1, length(idx.opt))]
      
      idx.good = c(idx.good, x.tmp$idx[idx.opt])
      
    }
    
    x.res = x.gap[idx.good,]
    rmSafe(x.gap)
    
    # Clean overlaps from both (base and query) sides
    x.res = cleanOverlaps(x.res)
    
  } else {
    # Create empty
    x.res <- data.frame(matrix(NA, nrow = 0, ncol = length(colnames(x.sk)), 
                               dimnames = list(NULL, colnames(x.sk))))
  }
  
  # ---- Read additional alignments ----
  
  file.gaps.out = paste0(path.gaps,
                         'acc_', acc, 
                         '_qchr_', query.chr, '_bchr_', base.chr, '_residual_out.txt', collapse = '')
  
  if(file.exists(file.gaps.out)) {
    # Fill up and down
    ## ----  Prepare ----
    x.sk = x.sk[order(x.sk$V2),]  # because all alignment will be according to the sorted "query"
    
    # Find positions of blocks, so that only at block edges somethings is changing, so we need not-trivial gaps
    x.sk$bl.end <- c(diff(x.sk$block.id) != 0, 1)
    x.sk$bl.beg <- c(1, diff(x.sk$block.id) != 0)
    x.sk$idx = 1:nrow(x.sk)
    
    y.tops = x.sk$p.beg  # top positions
    y.bots = x.sk$p.end  # bottom positions
    
    # Read blast results on "between blocks"
    pokaz('Read blast of "bad" gaps..', file.gaps.out)
    pokaz('File:', file.gaps.out)
    x.gap = read.table(file.gaps.out, stringsAsFactors = F)
    
    # Get real positions of fragments
    x.gap$q.beg = as.numeric(sapply(x.gap$V1, function(s) strsplit(s, '\\|')[[1]][2])) - 1
    x.gap$b.beg = as.numeric(sapply(x.gap$V10, function(s) strsplit(s, '\\|')[[1]][2])) - 1
    x.gap$V2 = x.gap$V2 + x.gap$q.beg
    x.gap$V3 = x.gap$V3 + x.gap$q.beg
    x.gap$V4 = x.gap$V4 + x.gap$b.beg
    x.gap$V5 = x.gap$V5 + x.gap$b.beg
    x.gap$idx = -(1:nrow(x.gap))
    x.gap$dir = (x.gap$V4 > x.gap$V5) * 1  # BUT DON'T change the order of V4 and V5
    x.gap$bl.beg = -1
    x.gap$bl.end = -1
    
    
    # Combine new core skeleton and new gaps
    comb.names = intersect(colnames(x.sk), colnames(x.gap))
    x.comb = rbind(x.sk[, comb.names], 
                   x.gap[, comb.names])
    x.comb = x.comb[order(x.comb$V2),]  # This sorting is relevant, because x.sk was initially also sorted!!
    
    idx.end = which(x.comb$bl.end == 1)
    idx.beg = which(x.comb$bl.beg == 1)
    
    ## ---- Remain in "up" direction ----
    idx.remain = c()
    for(i.cur in idx.end){
      # pokaz(i.cur)
      # if(x.comb$dir[i.cur] == 1) next
      n.gaps = countValueStretch(x.comb$bl.end, i.cur)
      if(n.gaps == 0) next
      
      # Get all possible BLAST-hits, which can be "after" the block.
      x.tmp = x.comb[i.cur + (0:(n.gaps)),]
      
      if(nrow(x.tmp) > complexity.threshold) next
      
      x.tmp$w = abs(x.tmp$V3 - x.tmp$V2) + abs(x.tmp$V4 - x.tmp$V5)
      x.tmp$w[1] = 0
      
      # # Remain the proper direction
      # x.tmp = x.tmp[x.tmp$dir == x.tmp$dir[1],]
      
      if(x.comb$dir[i.cur] == 0){
        # Find next "top" base, so you should not consider BLAST-hits higher than this
        y.top = min(c(Inf, y.tops[y.tops > x.comb$V5[i.cur]]))
        y.bot = min(x.comb$V4[i.cur], x.comb$V5[i.cur])  # min - потому что хочу взять и саму затравку
        
        x.tmp = x.tmp[(x.tmp$V4 >= y.bot) & (x.tmp$V5 >= y.bot),]
        x.tmp = x.tmp[(x.tmp$V4 <= y.top) & (x.tmp$V5 <= y.top),]
        if(nrow(x.tmp) <= 1) next
        
        # Run path search
        idx.visit = pathUpPlus(x.tmp)
      } else {
        
        # Find next "bottom" base, so you should not consider BLAST-hits higher than this
        y.bot = max(c(-Inf, y.bots[y.bots < x.comb$V5[i.cur]]))
        y.top = max(x.comb$V4[i.cur], x.comb$V5[i.cur])
        
        x.tmp = x.tmp[(x.tmp$V4 >= y.bot) & (x.tmp$V5 >= y.bot),]
        x.tmp = x.tmp[(x.tmp$V4 <= y.top) & (x.tmp$V5 <= y.top),]
        
        if(nrow(x.tmp) <= 1) next
        
        idx.visit = pathUpMinus(x.tmp)
      }
      
      # If found something - add
      if(!is.null(idx.visit)){
        idx.remain = c(idx.remain, x.tmp$idx[idx.visit][-1])
      }
    }
    
    ## ---- Remain in "Downwards" direction ----
    for(i.cur in idx.beg){
      
      # if(x.comb$dir[i.cur] == 1) next
      n.gaps = countValueStretchBW(x.comb$bl.end, i.cur)
      if(n.gaps == 0) next
      
      
      # Get all possible BLAST-hits, which can be "after" the block.
      x.tmp = x.comb[i.cur + ((-n.gaps):0),]
      
      if(nrow(x.tmp) > complexity.threshold) next
      
      x.tmp$w = abs(x.tmp$V3 - x.tmp$V2) + abs(x.tmp$V4 - x.tmp$V5)
      x.tmp$w[1] = 0
      
      # # Remain the proper direction - NOT NEEDED, because I've changed to "top" and "bottom"
      # x.tmp = x.tmp[x.tmp$dir == x.comb$dir[i.cur],]
      
      if(x.comb$dir[i.cur] == 0){
        # Find next "bottom" base, so you should not consider BLAST-hits higher than this
        y.bot = max(c(-Inf, y.bots[y.bots < x.comb$V5[i.cur]]))
        y.top = max(x.comb$V4[i.cur], x.comb$V5[i.cur])
        
        x.tmp = x.tmp[(x.tmp$V4 >= y.bot) & (x.tmp$V5 >= y.bot),]
        x.tmp = x.tmp[(x.tmp$V4 <= y.top) & (x.tmp$V5 <= y.top),]
        
        if(nrow(x.tmp) <= 1) next
        
        # Run path search
        idx.visit = pathDownPlus(x.tmp)
        
      } else {
        # Find next "top" base, so you should not consider BLAST-hits higher than this
        y.top = min(c(Inf, y.tops[y.tops > x.comb$V5[i.cur]]))
        y.bot = min(x.comb$V4[i.cur], x.comb$V5[i.cur])  # min - потому что хочу взять и саму затравку
        
        x.tmp = x.tmp[(x.tmp$V4 >= y.bot) & (x.tmp$V5 >= y.bot),]
        x.tmp = x.tmp[(x.tmp$V4 <= y.top) & (x.tmp$V5 <= y.top),]
        
        if(nrow(x.tmp) <= 1) next
        # stop('DownMinus direction found')
        
        # Run path search
        idx.visit = pathDownMinus(x.tmp)
        
      }
      
      # If found something - add
      if(!is.null(idx.visit)){
        idx.remain = c(idx.remain, x.tmp$idx[idx.visit])
      }
    }
    
    if(length(idx.remain) != 0){
      x.bw = x.gap[abs(idx.remain),]
      # Clean overlaps from both base and query sides
      x.bw = cleanOverlaps(x.bw)
    } else {
      x.bw = x.gap[F,]
    }
    
    rmSafe(x.gap)
    
  } else {
    
    # Create empty
    x.bw <- data.frame(matrix(NA, nrow = 0, ncol = length(colnames(x.sk)), 
                              dimnames = list(NULL, colnames(x.sk))))
    
  }
  
  
  # ---- Combine all together ----
  
  comb.names = intersect(intersect(colnames(x.sk), colnames(x.bw)), colnames(x.res))
  x.comb = rbind(x.res[,comb.names], x.sk[,comb.names])
  x.comb = rbind(x.comb, x.bw[,comb.names])
  x.comb = x.comb[order(x.comb$V2),]
  rownames(x.comb) = NULL
  
  # ---- Check uniqueness ---- 
  pos.q.occup = rep(0, 40000000)
  x.sk1 = x.comb
  # x.sk1 = x.res
  # x.sk1 = x.bw
  for(irow in 1:nrow(x.sk1)){
    pp = x.sk1$V4[irow]:x.sk1$V5[irow]
    if(sum(pos.q.occup[pp]) > 0) stop('non-unique')
    pos.q.occup[pp] = pos.q.occup[pp] + 1
  }
  sum(pos.q.occup > 1)
  sum(pos.q.occup)
  
  
  # ---- Check genomes ---- 
  # x.dir = setDir(x.comb, base.len = base.len)
  # checkCorrespToGenome(x.dir, query.fas = query.fas.chr,
  #                      base.fas.fw = base.fas.fw,
  #                      base.fas.bw = base.fas.bw)
  
  saveRDS(object = x.comb, file = file.aln.full) 
  
  suppressWarnings(rm(list = c("x.sk", "x.res", "x.bw", "base.fas.bw", "base.fas.fw", "query.fas.chr")))
  gc()
}
  

# ***********************************************************************
# ---- Loop  ----


if(num.cores == 1){
  for(f.maj in files.maj){
    loop.function(f.maj)
  }
} else {
  # Set the number of cores for parallel processing
  myCluster <- makeCluster(num.cores, type = "PSOCK") 
  registerDoParallel(myCluster) 
  
  tmp = foreach(f.maj = files.maj, .packages=c('crayon','stringr','Biostrings', 'seqinr'))  %dopar% { 
                              loop.function(f.maj)
                            }
  stopCluster(myCluster)
}



