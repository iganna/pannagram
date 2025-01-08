#' @export
cutAln <- function(path.msa, i.chr, p.beg, p.end, 
                   acc = NULL,
                   aln.type=NULL, ref.acc='',
                   mode = 'seq'){
  
  s.pangenome = c('pangen', 'pannagram', 'pangenome')
  
  gr.accs.e = "accs/"
  gr.accs.b <- "/accs"
  
  if(ref.acc == ''){
    ref.suff = ref.acc
  } else {
    # ref.suff = paste0('_ref_', ref.acc)
    ref.suff = paste0('_', ref.acc)
  }
  
  
  if(is.null(aln.type)) stop('Provide the alignment type')
  file.msa = paste0(path.msa, aln.type, i.chr, '_', i.chr, ref.suff, '.h5')
  
  if(mode == 'seq'){
    file.mode = paste0(path.msa, 'seq/', 'seq_', i.chr, '_', i.chr, ref.suff, '.h5')
  } else if(mode == 'pos') {
    file.mode = paste0(path.msa, aln.type, i.chr, '_', i.chr, ref.suff, '.h5')
  } else {
    stop("Mode could be either 'seq' or 'pos'")
  }
  
  if(!file.exists(file.msa)) stop(paste('File', file.msa, 'does not exist'))
  
  if(!is.null(acc)){
    pokaz('Define new pos based on the accession', acc)
    
    if(tolower(acc) %in% s.pangenome){
      info = h5ls(file.msa)
      info = info[info$group == '/accs',]
      v = 1:as.numeric(info$dim[1])
    } else {
      v = h5read(file.msa, paste0(gr.accs.e, acc))  
    }
    
    
    p.beg.acc = which(v == p.beg)
    p.end.acc = which(v == p.end)
    
    if(length(p.beg.acc) == 0) stop(paste('Position', p.beg, 'in not found in the alignment of the accession', acc))
    if(length(p.end.acc) == 0) stop(paste('Position', p.end, 'in not found in the alignment of the accession', acc))
    
    v = v[p.beg.acc:p.end.acc]
    v = v[v != 0]
    if(is.unsorted(v)) stop('The region is not in one synteny block')
    
    p.beg = p.beg.acc
    p.end = p.end.acc
  }
  
  groups = h5ls(file.mode)
  accessions <-  groups$name[groups$group == gr.accs.b]
  aln.mx = c()
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.mode, paste0(gr.accs.e, acc))
    aln.mx = rbind(aln.mx, v[p.beg:p.end])
  }
  rownames(aln.mx) <- accessions
  return(aln.mx)
}