#' Convert alignment HDF5 files between the legacy vector format and the interval one
#'
#' Legacy: /accs/<acc> is a vector along the pangenome coordinate.
#' Interval (_n): /accs/<acc> is an integer table pan.beg | acc.beg | len | delta | block,
#'                and the length of the pangenome axis is the dataset attribute `len.pan`.
#'
#' Usage:
#'   Rscript convert_h5_n.R --file.in X.h5 --file.out Y.h5 [--direction to.ivl|to.vec]
#'   Rscript convert_h5_n.R --path.in DIR  --path.out DIR  [--direction ...] [--pattern '\\.h5$']
#'
#' Everything outside /accs (/len, /ref, /trust, ...) is copied unchanged. The legacy
#' /blocks group is dropped when converting to the interval format: the block id
#' becomes a column of every interval table.

suppressMessages({
  library(optparse)
  library(rhdf5)
})

source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("utils/chunk_hdf5.R", package = "pannagram"))
source(system.file("utils/interval_func.R", package = "pannagram"))

option_list <- list(
  make_option("--file.in",   type = "character", default = NULL, help = "Input .h5 file"),
  make_option("--file.out",  type = "character", default = NULL, help = "Output .h5 file"),
  make_option("--path.in",   type = "character", default = NULL, help = "Input directory"),
  make_option("--path.out",  type = "character", default = NULL, help = "Output directory"),
  make_option("--pattern",   type = "character", default = "\\.h5$", help = "File pattern for --path.in"),
  make_option("--direction", type = "character", default = "to.ivl", help = "'to.ivl' or 'to.vec'")
)
opt <- parse_args(OptionParser(option_list = option_list))

if(!(opt$direction %in% c("to.ivl", "to.vec"))) stop("--direction must be 'to.ivl' or 'to.vec'")

convertOne <- function(f.in, f.out, direction){
  if(file.exists(f.out)) file.remove(f.out)
  h5createFile(f.out)
  suppressMessages(h5createGroup(f.out, gr.accs.e))

  root <- h5ls(f.in, recursive = FALSE)
  for(nm in setdiff(root$name, c(sub("/", "", gr.accs.b), sub("/$", "", gr.blocks)))){
    suppressMessages(h5write(h5read(f.in, nm), f.out, nm))
  }

  accs <- h5AccNames(f.in)
  n.ivl <- 0
  for(acc in accs){
    if(direction == "to.ivl"){
      v <- as.numeric(h5read(f.in, paste0(gr.accs.e, acc)))
      v[is.na(v)] <- 0
      h5VecWrite(f.out, acc, v)
      n.ivl <- n.ivl + nrow(ivlFromVec(v, blocks = FALSE))
    } else {
      v <- h5VecRead(f.in, acc)
      suppressMessages(h5write(v, f.out, paste0(gr.accs.e, acc)))
    }
  }

  # The pangenome length is implicit in the legacy dim, explicit in the interval format
  if(direction == "to.ivl" && !("len" %in% root$name) && length(accs) > 0){
    h5PanLenSet(f.out, h5AccLen(f.out, accs[1]))
  }
  H5close()
  list(n.acc = length(accs), n.ivl = n.ivl,
       size.in = file.size(f.in), size.out = file.size(f.out))
}

files <- NULL
if(!is.null(opt$file.in)){
  if(is.null(opt$file.out)) stop("--file.out is required with --file.in")
  files <- data.frame(i = opt$file.in, o = opt$file.out, stringsAsFactors = FALSE)
} else if(!is.null(opt$path.in)){
  if(is.null(opt$path.out)) stop("--path.out is required with --path.in")
  if(!dir.exists(opt$path.out)) dir.create(opt$path.out, recursive = TRUE)
  nm <- list.files(opt$path.in, pattern = opt$pattern)
  if(length(nm) == 0) stop(paste("No files matching", opt$pattern, "in", opt$path.in))
  files <- data.frame(i = file.path(opt$path.in, nm), o = file.path(opt$path.out, nm),
                      stringsAsFactors = FALSE)
} else stop("Provide either --file.in/--file.out or --path.in/--path.out")

for(k in seq_len(nrow(files))){
  r <- convertOne(files$i[k], files$o[k], opt$direction)
  pokaz(basename(files$i[k]), '->', basename(files$o[k]),
        ': accs', r$n.acc,
        if(opt$direction == "to.ivl") paste0(', intervals ', r$n.ivl) else '',
        ',', round(r$size.in/1e6, 2), 'MB ->', round(r$size.out/1e6, 2), 'MB')
}
