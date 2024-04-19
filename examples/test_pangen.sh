./pangen.sh -path_out ../pan_test/rhiz1/ -refs 'ref_1021,A18' -nchr_query 1 -nchr_ref 1 -path_in ../rhizobia/ -sort_chr_len T -cores 30 -s 4

./pangen.sh -path_out '../pan_test/anna_50' -refs 'NT1_50,TE11_50' -nchr_ref 8 -path_in '../lyrata_updated/' -nchr_query 8 -cores 30 -one2one
./pangen.sh -path_out '../pan_test/anna_norm' -refs 'NT1_220222,TE11_final' -nchr_ref 8 -path_in '../lyrata_updated/' -nchr_query 8 -cores 30 -one2one

./analys.sh -path_out '../pan_test/tom2'  -ref 0 -blocks
./sv.sh -path_out '../pan_test/tom2'  -ref 0 -gff
 ./sv.sh -path_out '../pan_test/tom2'  -ref 0 -te -te_file ../new_genes/new_te.fasta
 ./sv.sh -path_out '../pan_test/tom2'  -ref 0 -graph


./sim_search.sh -in ../new_genes/new_genes.fasta -genome ../pb_chromosomes/0_chr1.fasta -out ../tmp.txt

./pipeline_test.sh -path_out '../pan_test/tom2' -refs '0,6046-v1.1,6191-v1.1' -nchr_ref 5 -path_in '../pb_updated/' -nchr_query 5 -cores 30 -acc_anal acc_tom.txt -one2one

./pangen.sh -path_out '../pan_test/p27/' -refs '0,10024' -nchr_ref 5 -path_in '../pb_27/' -nchr_query 5 -cores 30  -one2one
./analys.sh -path_out '../pan_test/p27/'  -ref 0 -blocks -seq -cores 30
./sv.sh -path_out '../pan_test/p27/'  -ref 0  -gff -te -te_file ../new_genes/new_te.fasta -graph 