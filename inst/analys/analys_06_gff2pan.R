# Convert gff files into pangenome coordinates


# ***********************************************************************
# ---- Libraries and dependencies ----
library(crayon)
library(rhdf5)
source(system.file("utils/utils.R", package = "pannagram"))
source(system.file("pangen/comb_func.R", package = "pannagram"))
source(system.file("analys_func.R", package = "pannagram"))

# ***********************************************************************
# ---- Setup ----

# 
# # Set the working path for annotations
path.msa = '/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/pushkin/genomes_GR3013/data/msa/'
path.annot = '/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/pushkin/genomes_GR3013/data/gff/'
path.res = '/Users/annaigolkina/Library/CloudStorage/OneDrive-Personal/pushkin/genomes_GR3013/data/gff_common/'
if (!dir.exists(path.res)) {
  dir.create(path.res, recursive = TRUE)
}

ref.pref = '0'
aln.type = 'v_'
aln.type = 'msa_'

# path.msa = '/Volumes/Samsung_T5/vienn/msa/'
# path.annot = '/Volumes/Samsung_T5/vienn/annotation/'
# path.res = '/Volumes/Samsung_T5/vienn/annotation_common/'
if (!dir.exists(path.res)) {
  dir.create(path.res, recursive = TRUE)
} 


# ---- Chromosome name format to extract the chromosome number ----
# s.chr = '_Chr' # in this case the pattern is "*_ChrX", where X is the number
s.chr = '_' # in this case the pattern is "*_X", where X is the number
n.chr = 2
ref.acc = 'GR3013_prokka'
# ---- Accessions ----
files.gff <- list.files(path = path.annot, pattern = "\\.gff$", full.names = FALSE)

accessions <- sub("\\.gff$", "", files.gff)

pokaz('Accessions, (amount', length(accessions), ') :')
pokaz('  ', accessions)


# ---- Variables ----
gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'

gr.blocks = 'blocks/'

s.pannagram = 'Pannagram'


# ***********************************************************************

# Loop through each accession in the list of accessions to work with
for(acc in accessions){
  
  pokaz('Accession', acc)
  
  # Read the GFF file 
  gff.acc = read.table(paste0(path.annot, acc,'.gff'), stringsAsFactors = F)
  gff.pan = gff2gff(path.msa, acc1 = acc, acc2 = 'PanGen', pangenome.name='PanGen',
                    gff1 = gff.acc, ref.acc = ref.acc,
                    s.chr = s.chr,n.chr = n.chr,exact.match = F)
  gff.pan$V1 = gsub('contig_', 'PanGen_Chr',gff.pan$V1)
  writeGFF(gff.pan, paste0(path.res, acc,'.gff'))
  
}


