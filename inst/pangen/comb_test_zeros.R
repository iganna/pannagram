library(rhdf5)
library(pannagram)
library(crayon)
gr.accs.b <- "/accs"
gr.accs.e <- "accs/"

path.cons = '/Volumes/Samsung_T5/vienn/test/a27/intermediate/consensus3/'
# path.cons = '/Volumes/Samsung_T5/vienn/test/manuals/symA_out/intermediate/consensus/'
i.chr = 2
aln.type.in = 'msa_'

s.comb = paste0(i.chr, '_', i.chr)

file.comb = paste0(path.cons, aln.type.in, s.comb,'.h5')

groups = h5ls(file.comb)
accessions = groups$name[groups$group == gr.accs.b]
n.acc = length(accessions)
pokaz('Number of accessions', n.acc)

# Get breaks
idx = 0
v.mx = c()
pos.beg = 0
pos.end = 0
for(acc in accessions){
  pokaz(acc)
  s.acc = paste0(gr.accs.e, acc)
  v = h5read(file.comb, s.acc)
  v[is.na(v)] = 0
  idx = idx + abs(v)

  if(pos.beg == 0){
    pos.beg = which(v == 5469954)
    pos.end = which(v == 5506931)
    # stop()
  }
  v.mx = cbind(v.mx, v[pos.beg:pos.end])
}

min(idx)

which(idx == 0)