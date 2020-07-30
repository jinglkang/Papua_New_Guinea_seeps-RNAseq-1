# papue_new_guinea-RNA-seq

# 1. after the running of OrthoFinder, select Orthogroup that contain at least one reads per species;
# and blastp to swiss-prot database, Orthogroups/Orthogroups.GeneCount.tsv is OrthoFinder result file;
# the qualified orthogroup will be in ./orthogroup_needed_sequences;
# the blastp result will be in blastp_result
perl blastp_uni.pl Orthogroups/Orthogroups.GeneCount.tsv

# 2. select the most representative transcript for each species in each orthologous group
perl get_best_blast_orthogroup-5.pl -blast_result=blastp_result -reads_matrix=six_fishes_reads_num
