library(rhdf5)
library(doParallel)
library(foreach)

chunk_len <- 16384

n.chr = 5

ncores <- n.chr

# aln.pref = 'pan'
# storage_mode = 'integer'

aln.pref = 'seq'
storage_mode = 'character'


cl <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(cl)

foreach(i.chr = 2:n.chr, .packages = "rhdf5") %dopar% {
  s.comb <- paste0(i.chr, "_", i.chr)
  
  infile  <- paste0(aln.pref, "_", s.comb, ".h5")
  outfile <- paste0(aln.pref, "_", s.comb, "_chunked.h5")
  
  accs <- h5ls(infile, recursive = TRUE)
  accs <- accs[accs$group == "/accs", ]
  ids <- accs$name
  
  if (file.exists(outfile)) file.remove(outfile)
  
  h5createFile(outfile)
  h5createGroup(outfile, "accs")
  
  root_items <- h5ls(infile, recursive = FALSE)
  
  if ("len" %in% root_items$name) {
    h5write(h5read(infile, "/len"), outfile, "/len")
  }
  
  if ("ref" %in% root_items$name) {
    h5write(h5read(infile, "/ref"), outfile, "/ref")
  }
  
  for (id in ids) {
    old_path <- paste0("/accs/", id)
    new_path <- old_path
    
    cat("Processing", s.comb, id, "...\n")
    
    x <- h5read(infile, old_path)
    x_len <- length(x)
    
    if (storage_mode == "character") {
      
      max_len <- 1
      
      h5createDataset(
        file = outfile,
        dataset = new_path,
        dims = x_len,
        chunk = chunk_len,
        storage.mode = "character",
        size = max_len,
        level = 6,
        filter = "gzip",
        shuffle = TRUE
      )
      
    } else if (storage_mode == "integer") {
      h5createDataset(
        file = outfile,
        dataset = new_path,
        dims = x_len,
        chunk = chunk_size,
        storage.mode = "integer",
        level = 6,
        filter = "gzip",
        shuffle = TRUE
      )
      
    } else {
      stop("Unsupported storage_mode: ", storage_mode)
    }
    
    h5write(x, outfile, new_path)
    
    rm(x)
    gc()
  }
}

parallel::stopCluster(cl)
