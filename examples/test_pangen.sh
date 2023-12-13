./pangen_consensus.sh -pref_global ../pan_test/rhiz/ -ref_set 'ref_1021,A18' -n_chr_query 1 -n_chr_ref 1 -path_in ../rhizobia/ -sort_chr_len T -cores 30 -s 4

./pangen_consensus.sh -pref_global '../pan_test/anna_50' -ref_set 'NT1_50,TE11_50' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 30 -all_cmp F
./pangen_consensus.sh -pref_global '../pan_test/anna_norm' -ref_set 'NT1_220222,TE11_final' -n_chr_ref 8 -path_in '../lyrata_updated/' -n_chr_query 8 -cores 30 -all_cmp F

./analys.sh -pref_global '../pan_test/tom2'  -ref_pref 0 -blocks
./sv.sh -pref_global '../pan_test/tom2'  -ref_pref 0 -gff
 ./sv.sh -pref_global '../pan_test/tom2'  -ref_pref 0 -te -te_file ../new_genes/new_te.fasta
 ./sv.sh -pref_global '../pan_test/tom2'  -ref_pref 0 -graph


./sim_search.sh -in ../new_genes/new_genes.fasta -genome ../pb_chromosomes/0_chr1.fasta -out ../tmp.txt

./pipeline_test.sh -pref_global '../pan_test/tom2' -ref_set '0,6046-v1.1,6191-v1.1' -n_chr_ref 5 -path_in '../pb_updated/' -n_chr_query 5 -cores 30 -acc_anal acc_tom.txt -all_cmp F

./pipeline_test.sh -pref_global '../pan_test/p27/' -ref_set '0,10024' -n_chr_ref 5 -path_in '../pb_27/' -n_chr_query 5 -cores 30  -all_cmp F