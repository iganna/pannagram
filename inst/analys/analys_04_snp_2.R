suppressMessages({
  library(Biostrings)
  library(rhdf5)
  library(foreach)
  library(doParallel)
  library(optparse)
  library(crayon)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("analys/analys_func.R", package = "pannagram"))
source(system.file("utils/chunk_hdf5.R", package = "pannagram")) 

args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option("--path.features.msa", type = "character", default = NULL, help = "Path to msa dir (features)"),
  make_option("--path.snp",          type = "character", default = NULL, help = "Path to snp dir"),
  make_option("--path.seq",          type = "character", default = NULL, help = "Path to seq dir"),
  make_option("--cores",             type = "integer", default = 1,       help = "number of cores to use for parallel processing"),
  make_option("--aln.type",          type="character", default="default", help="type of alignment ('msa_', 'comb_', 'v_', etc)"),
  make_option("--ref",               type="character", default=NULL,      help="prefix of the reference file"),
  make_option("--acc",               type="character", default='',      help="For which accession to generate the VCF file")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser, args = args)

source(system.file("utils/chunk_logging.R", package = "pannagram"))

num.cores <- opt$cores

path.features.msa <- if (!is.null(opt$path.features.msa)) opt$path.features.msa else stop("Error: 'path.features.msa' is NULL. Please provide a valid path.")

path.seq <- opt$path.seq
if (!dir.exists(path.seq)) stop("Folder `path.seq` does not exist")

path.snp <- opt$path.snp
if (!dir.exists(path.snp)) {
  dir.create(path.snp)
}
if (!dir.exists(path.snp)) {
  stop("The output folder was not created")
}

if (!dir.exists(path.features.msa)) {
  stop(paste("The consensus folder does not exist:", path.features.msa))
}

if (!is.null(opt$aln.type)) {
  aln.type = opt$aln.type
} else {
  aln.type = aln.type.msa
}

ref.name <- opt$ref
if (ref.name == "NULL" || is.null(ref.name)) ref.name <- ""

source(system.file("utils/chunk_combinations.R", package = "pannagram")) 

acc.vcf <- opt$acc
if(acc.vcf == ''){
  pokaz('Empty acc')
}


# --------------------------------------------------
# main loop by s.comb, parallel inside by acc
# --------------------------------------------------
for (s.comb in s.combinations) {
  
  pokaz("Combination", s.comb)
  
  # Get Consensus
  i.chr = comb2ref(s.comb)
  file.seq.cons = paste0(path.seq, "seq_cons_", s.comb, ref.suff, ".fasta")
  if(!file.exists(file.seq.cons)){
    stop("Consensus fasta does not exist")
  }
  s.pangen = readFastaMy(file.seq.cons)
  s.pangen.name = names(s.pangen)[1]
  s.pangen = seq2nt(s.pangen)
  
  # Get accessions
  file.seq = paste0(path.seq, "seq_", s.comb, ref.suff, ".h5")
  
  groups = h5ls(file.seq)
  accessions = groups$name[groups$group == gr.accs.b]
  n.acc = length(accessions)
  
  pokaz("Round 1: get positions of differences..")
  
  # -----------------------------
  # ROUND 1: parallel by acc
  # -----------------------------
  if (num.cores == 1) {
    pos.diff.list <- lapply(accessions, function(acc) {
      pokaz("Difference in accession", acc)
      v = h5read(file.seq, paste0(gr.accs.e, acc))
      pos = which((v != s.pangen) & (v != "-"))
      rmSafe(v)
      gc()
      pos
    })
  } else {
    myCluster <- makeCluster(num.cores, type = "PSOCK")
    registerDoParallel(myCluster)
    
    pos.diff.list <- foreach(
      acc = accessions,
      .packages = c("rhdf5", "crayon", "pannagram"),
      .errorhandling = "stop"
    ) %dopar% {
      pokaz("Difference in accession", acc)
      v = h5read(file.seq, paste0(gr.accs.e, acc))
      pos = which((v != s.pangen) & (v != "-"))
      rmSafe(v)
      gc()
      pos
    }
    
    stopCluster(myCluster)
  }
  
  pos.diff = sort(unique(unlist(pos.diff.list, use.names = FALSE)))
  
  # Checkup
  # lens <- lengths(pos.diff.list)
  # pokaz("Min length pos.diff.list:", min(lens), "\n")
  # pokaz("Max length pos.diff.list:", max(lens), "\n")
  # pokaz("Length pos.diff:", length(pos.diff), "\n")
  
  if (length(pos.diff) == 0) {
    pokaz("No SNPs were found..")
    next
  }
  
  pokaz("Round 2: get diffs in common positions..")
  pos = pos.diff
  
  # -----------------------------
  # ROUND 2: parallel by acc
  # -----------------------------
  if (num.cores == 1) {
    res.list <- lapply(accessions, function(acc) {
      pokaz("Sequence of accession", acc)
      
      v = h5read(file.seq, paste0(gr.accs.e, acc))
      
      val = v[pos]
      
      rmSafe(v)
      gc()
      
      list(acc = acc, val = val)
    })
  } else {
    myCluster <- makeCluster(num.cores, type = "PSOCK")
    registerDoParallel(myCluster)
    
    res.list <- foreach(
      acc = accessions,
      .packages = c("rhdf5", "crayon", "pannagram"),
      .errorhandling = "stop"
    ) %dopar% {
      pokaz("Sequence of accession", acc)
      
      v = h5read(file.seq, paste0(gr.accs.e, acc))
      
      val = v[pos]
      
      rmSafe(v)
      gc()
      
      list(acc = acc, val = val)
    }
    
    stopCluster(myCluster)
  }
  
  # Build matrices
  pokaz("Build matrices...")
  acc.names = vapply(res.list, `[[`, character(1), "acc")
  snp.ref = s.pangen[pos]
  snp.val = do.call(cbind, lapply(res.list, `[[`, "val"))
  colnames(snp.val) = acc.names
  
  # Clean up the memory
  rm(res.list)
  
  pokaz("Save VCF-file...")
  file.vcf = paste0(path.snp, "snps_", s.comb, ref.suff, "_pangen.vcf")
  saveVCF2(snp.val, pos, chr.name = paste0("PanGen_Chr", i.chr), file.vcf = file.vcf,
           snp.ref = snp.ref)
  
  # Create the VCF-file for the reference accession
  file.comb = paste0(path.features.msa, aln.pref, s.comb, ref.suff, ".h5")
  
  acc = ''
  if (ref.name != "") {
    acc = ref.name
  } 
  
  if(acc.vcf != ""){
    acc = acc.vcf
  }
  
  if(acc == ''){
    quit(save = "no")
  }
  
  pokaz('Generating VCF-file for accession', acc)
  
  pos.acc = h5read(file.comb, paste0(gr.accs.e, acc))
  pos.acc = pos.acc[pos]
  snp.val.acc = snp.val[pos.acc != 0, , drop = FALSE]
  snp.ref.acc = snp.val.acc[, acc]
  pos.acc = abs(pos.acc[pos.acc != 0])
  
  # Sort positions
  ord = order(pos.acc)
  pos.acc = pos.acc[ord]
  snp.val.acc = snp.val.acc[ord, , drop = FALSE]
  snp.ref.acc = snp.ref.acc[ord]
  
  pokaz("Save VCF-file for the accession", acc, "...")
  file.vcf.acc = paste0(path.snp, "snps_", s.comb, ref.suff, "_", acc, ".vcf")
  saveVCF2(snp.val.acc, pos.acc, chr.name = paste0(acc, "_Chr", i.chr), file.vcf = file.vcf.acc,
           snp.ref = snp.ref.acc)
}

warnings()