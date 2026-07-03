# Combine all alignments together into the final one

suppressMessages({
library(rhdf5)
library('foreach')
library(doParallel)
library("optparse")
source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func_mafft_refine2.R", package = "pannagram"))
})

# ***********************************************************************
# ---- Command line arguments ----

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--path.cons",      type = "character", default = NULL, help = "Path to directory with the consensus"),
  make_option("--cores",          type = "integer",   default = 1,    help = "Number of cores to use for parallel processing"),
  make_option("--path.log",       type = "character", default = NULL, help = "Path for log files"),
  make_option("--log.level",      type = "character", default = NULL, help = "Level of log to be shown on the screen")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, args = args);



# ***********************************************************************
# ---- Logging ----

source(system.file("utils/chunk_logging.R", package = "pannagram")) # a common code for all R logging

# ---- HDF5 ----

source(system.file("utils/chunk_hdf5.R", package = "pannagram")) # a common code for variables in hdf5-files


# ***********************************************************************
# ---- Values of parameters ----

# Number of cores for parallel processing
num.cores <- opt$cores

if (!is.null(opt$path.cons)) path.cons <- opt$path.cons

# ***********************************************************************

# ---- Combinations of chromosomes query-base to create the alignments ----


aln.pref.msa = paste0(aln.type.msa, '_')
aln.pref.clean = paste0(aln.type.clean, '_')
aln.pref.add = paste0(aln.type.add, '_')

s.pattern <- paste0("^", aln.pref.msa, ".*")
files <- list.files(path = path.cons, pattern = s.pattern, full.names = FALSE)
pref.combinations = gsub(aln.pref.msa, "", files)
pref.combinations <- sub(".h5", "", pref.combinations)

pattern <- "^[0-9]+_[0-9]+$"
pref.combinations <- pref.combinations[grep(pattern, pref.combinations)]

pokaz('Combinations', pref.combinations, file=file.log.main, echo=echo.main)


# ***********************************************************************
# ---- MAIN program body ----

stat.comb <- data.frame(comb = character(),
                        coverage = numeric(),
                        stringsAsFactors = FALSE)


# Rebuild the set of completed items from all worker logs (any core count)
done.set <- getDoneSet(path.log)
if(length(done.set) > 0){
  pokaz('Skip already done:', length(done.set), file=file.log.main, echo=echo.main)
}
assign('.worker.id', 1, envir = .GlobalEnv)   # sequential outer loop -> single core_1.log

for(s.comb in pref.combinations){

  # ---- Checkpoint: skip already completed combinations ----
  s.comb.id <- s.comb
  if(s.comb.id %in% done.set) next

  # One log file per worker (bounded number of files); also the checkpoint ledger
  file.log.loop = initLoopLog(path.log)

  pokaz('* Combination', s.comb, file=file.log.main, echo=echo.main)

  # Get accessions
  file.cln = paste0(path.cons, aln.pref.clean, s.comb,'.h5')
  file.msa = paste0(path.cons, aln.pref.msa, s.comb,'.h5')

  file.add = paste0(path.cons, aln.pref.add, s.comb,'.h5')
  # Keep an existing h5 to resume; create it only if missing.
  if(!file.exists(file.add)){
    suppressMessages({
      h5createFile(file.add)
      h5createGroup(file.add, gr.blocks)
      h5createGroup(file.add, gr.accs.e)
    })
  }

  groups = h5ls(file.cln)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  base.len = length(h5read(file.cln, paste0(gr.accs.e, accessions[1])))
  
  ref = accessions[1]
  s.ref = paste0(gr.accs.e, ref)
  v.ref0 = h5read(file.cln, s.ref)
  v.ref1 = h5read(file.msa, s.ref)
  
  # Save the reference genome (idempotent)
  suppressMessages({
    try(h5delete(file.add, s.ref), silent = TRUE)
    h5write(v.ref1, file.add, s.ref)
  })

  for(acc in setdiff(accessions, ref)){

    # ---- Checkpoint: skip accessions already written for this combination ----
    acc.id <- paste0(s.comb.id, '_', acc)
    if(acc.id %in% done.set) next

    s.acc = paste0(gr.accs.e, acc)
    v0 =  h5read(file.cln, s.acc)
    v1 =  h5read(file.msa, s.acc)
    v1[is.na(v1)] = 0
    
    idx.add = which(!(v0 %in% v1) & (v0 != 0))
    v.add = data.frame(ref.val = v.ref0[idx.add], acc.val = v0[idx.add])
    
    pokaz(acc, 'extra positions:', nrow(v.add))
    
    v.add = v.add[v.add$ref.val != 0,]
    v.add = v.add[v.add$ref.val %in% v.ref1, ]
    
    pokaz(acc, 'after filtration:', nrow(v.add))
    
    if(nrow(v.add) == 0){
      suppressMessages({
        try(h5delete(file.add, s.acc), silent = TRUE)
        h5write(v1, file.add, s.acc)
      })
      markDone(acc.id, file=file.log.loop, echo=echo.loop)
      next
    }
    
    v.add$idx1 = match(v.add[,1], v.ref1)
    v.add$acc.val1 = v1[v.add[,3]]
    
    if(max(abs(v.add$acc.val1)) != 0) {
      # save(list = ls(), file = "tmp_workspace_insert.RData")
      pokaz('Something is wrong: max(abs(v.add$acc.val1)) = ', max(abs(v.add$acc.val1)) ) 
      stop()
    }
    
    v1[v.add$idx1] = v.add$acc.val
    
    if(length(unique(v1)) != (length(v1) - sum(v1 == 0) + 1)) stop('Non-unique positions are found')

    suppressMessages({
      try(h5delete(file.add, s.acc), silent = TRUE)
      h5write(v1, file.add, s.acc)
    })
    markDone(acc.id, file=file.log.loop, echo=echo.loop)

  }

  # ---- Checkpoint marker: combination fully processed ----
  markDone(s.comb.id, file=file.log.loop, echo=echo.loop)
}


