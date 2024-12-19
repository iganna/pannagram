#' @export
cutAln <- function(path.msa, i.chr, p.beg, p.end, 
                   acc = NULL,
                   aln.type=NULL, ref.acc='',
                   mode = 'seq'){
  
  gr.accs.e = "accs/"
  gr.accs.b <- "/accs"
  
  if(ref.acc == ''){
    ref.suff = ref.acc
  } else {
    # ref.suff = paste0('_ref_', ref.acc)
    ref.suff = paste0('_', ref.acc)
  }
  
  if(mode == 'seq'){
    if(is.null(aln.type)) aln.type = 'seq_'
    file.msa = paste0(path.msa, 'seq/', aln.type, i.chr, '_', i.chr, ref.suff, '.h5')
  } else if(mode == 'pos') {
    if(is.null(aln.type)) aln.type = 'msa_'
    file.msa = paste0(path.msa, aln.type, i.chr, '_', i.chr, ref.suff, '.h5')
  } else {
    stop("Mode could be either 'seq' or 'pos'")
  }
  
  if(!file.exists(file.msa)) stop(paste('File', file.msa, 'does not exist'))
  
  if(!is.null(acc)){
    pokaz('Define new pos based on the accession', acc)
    
    v = h5read(file.msa, paste0(gr.accs.e, acc))
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
  
  groups = h5ls(file.msa)
  accessions <-  groups$name[groups$group == gr.accs.b]
  aln.mx = c()
  for(acc in accessions){
    pokaz('Sequence of accession', acc)
    v = h5read(file.msa, paste0(gr.accs.e, acc))
    aln.mx = rbind(aln.mx, v[p.beg:p.end])
  }
  rownames(aln.mx) <- accessions
  return(aln.mx)
}