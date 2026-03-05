# ---- Prefix of the alignment type ----

aln.pref = paste0(aln.type, '_')

s.pattern <- paste0("^", aln.pref, ".*\\.h5$")
s.combinations <- list.files(path = path.features.msa, pattern = s.pattern, full.names = FALSE)
if (!length(s.combinations)) stop(paste("No .h5 files matching", aln.type, "prefix found."))
s.combinations = sub(aln.pref, "", s.combinations)
s.combinations = sub("\\.h5$", "", s.combinations)

# ---- Suffix of the reference ----

if(ref.name != ""){
  ref.suff = paste0('_', ref.name)
  
  pokaz('Reference:', ref.name)
  s.combinations <- s.combinations[grep(ref.suff, s.combinations)]
  
  if (!length(s.combinations)) stop(paste("No .h5 files matching", ref.suff, "suffix found."))
  
  s.combinations = gsub(ref.suff, "", s.combinations)
} else {
  ref.suff = ''
  s.combinations <- s.combinations[grepl("^[0-9]+_[0-9]+$", s.combinations)]
}

if(length(s.combinations) == 0){
  stop('No Combinations found.')
} else {
  if(!checkCombinations(s.combinations)){
    pokaz('Pattern', s.pattern)
    pokaz('Combinations:', s.combinations)
    stop("Wrong combination format.\nPossible hint: check that you have provided the name of the reference genome -ref.")
  }
  pokaz('Combinations:', s.combinations, file=file.log.main, echo=echo.main)
}

