




gr.accs.e <- "accs/"
gr.accs.b <- "/accs"
gr.break.e = 'break/'
gr.break.b = '/break'


path.cons = './'
s.comb = '1_1'
ref.pref = '0'

file.msa = paste(path.cons, 'msa_',s.comb,'_ref_',ref.pref,'.h5', sep = '')

groups = h5ls(file.comb)
accessions = groups$name[groups$group == gr.accs.b]
n.acc = length(accessions)

